using JuMP
using HiGHS
using SparseArrays
using Statistics

#Nb(m): number of thermal units in subsystem m
#Recent_History_Inflows(-t+2,:) inflows for time t (t=1,0,-1,-2,...)
#Tendancy_Inflows(t,m) is ct(m) for t=1,...,TabS(m),m=1,...,NS.
#sigmaInflows(t,m): s.d. of \xi_t(m)
#Probabilities(t,j): probability of realization Inflow_Noises(t,:,j) for period t+1,t=1,...,T-1,scenario j, sum(Probabilities(t,:))=1.

function sddp_hydro_interstage_dependent(Nb,talpha,T,S,NS,N,gamma,M,Inflow_Noises,Probabilities,Demand,Thermal_Costs,vminT,vmaxT,vminH,vmaxH,xmax,x0,lambda,evap,phi,TabS,TabP,tol,Tendancy_Inflows,Recent_History_Inflows,sigmaInflows,Itermax)

    zsup = Inf
    zinf = -Inf
    Iter = 1

    Zsups = []
    Zinfs = []

    Alphas = [[] for _ = 1:T-1]
    Thetas = [[] for _ = 1:T-1]

    FixedSubi = Int[]
    FixedSubj = Int[]
    FixedValij = Float64[]

    Nb_Thermal = sum(Nb)

    for i=1:NS
        FixedSubi = [FixedSubi; i; i; i]
        FixedSubj = [FixedSubj; Nb_Thermal+i; Nb_Thermal+NS+i; Nb_Thermal+2*NS+i]
        FixedValij = [FixedValij; 1; 1; 1]
    end
    
    for i=1:NS
        FixedSubi = [FixedSubi; (NS+i)*ones(Int,Nb[i]+1)]
        FixedSubj = [FixedSubj; Nb_Thermal+NS+i; collect(sum(Nb[1:i-1])+1:1:sum(Nb[1:i]))]
        FixedValij = [FixedValij; ones(Int,Nb[i]+1)]
    end

    DynamicSubi = [Int[] for _ = 1:T-1]
    DynamicSubj = [Int[] for _ = 1:T-1]
    DynamicValij = [Float64[] for _ = 1:T-1]

    taille = zeros(Int,T)
    StateOrder = zeros(Int,T,NS)
    periode = zeros(Int,T,NS)
    for t=1:T
        taille[T-t+1] = NS
        for m=1:NS
            periode[T-t+1,m] = (T-t+1) % TabS[m]
            if periode[T-t+1,m] == 0
                periode[T-t+1,m] = TabS[m]
            end
            if t==1
                StateOrder[T,m] = TabP[m][periode[T,m]]
            else
                StateOrder[T-t+1,m] = max(TabP[m][periode[T-t+1,m]],StateOrder[T-t+2,m]-1)
            end
            taille[T-t+1] = taille[T-t+1]+StateOrder[T-t+1,m]
        end
    end

    Cum_Probas = zeros(T-1,M+1)
    for t=1:T-1
        Cum_Probas[t,:] = [0;cumsum(Probabilities[t,:])]
    end

    Iter = 1
    End_Algo = true
    Costs = zeros(N)
    
    Time = @elapsed begin
        while End_Algo
            print(Iter, ",")
            #Iter
            #Forward pass
            Trial_States = zeros(T-1,NS,S)
            Sampled_Inflows = zeros(T,NS,S)
            for s=1:S
                ct=0
                for t=1:T
                    if t == 1
                        Sampled_Inflows[1,:,s] = Recent_History_Inflows[1,:]
                    else
                        Alea_Uniform = rand()
                        Index = sum(Alea_Uniform .>= Cum_Probas[t-1,:])
                        if Alea_Uniform == 1
                            Sampled_Noise = Inflow_Noises[t-1,:,M]
                        else
                            Sampled_Noise = Inflow_Noises[t-1,:,Index]
                        end
                        for m=1:NS
                            Sampled_Inflows[t,m,s] = Tendancy_Inflows[periode[t,m],m]
                            for j=1:TabP[m][periode[t,m]]
                                if t-j <= 1
                                    Sampled_Inflows[t,m,s] += phi[m][periode[t,m]][j] * Recent_History_Inflows[-t+j+2,m]
                                else
                                    Sampled_Inflows[t,m,s] += phi[m][periode[t,m]][j] * Sampled_Inflows[t-j,m,s]
                                end
                            end
                            Sampled_Inflows[t,m,s] = Sampled_Inflows[t,m,s] + sigmaInflows[periode[t,m],m] * Sampled_Noise[m]
                            if Sampled_Inflows[t,m,s] < 0
                                print("Negative inflows")
                                readline()
                            end
                        end
                    end

                    blx = []
                    bux = []
                    for m=1:NS
                        blx = [blx; vminT[t,m]]
                        bux = [bux; vmaxT[t,m]]
                    end
                    blx = [blx; lambda[t,:].*xmax[t,:]; vminH[t,:]; zeros(NS)]
                    bux = [bux; xmax[t,:]; vmaxH[t,:]; Inf*ones(NS)]
                    if t < T && Iter > 1
                        blx = [blx; -Inf]
                        bux = [bux; Inf]
                    end

                    if t == T || Iter == 1
                        nbvar = Nb_Thermal + 3 * NS
                        c = zeros(nbvar)
                        blc = zeros(2*NS)
                        buc = Inf*ones(2*NS)
                    else
                        nbvar = 1 + Nb_Thermal + 3 * NS
                        c = zeros(nbvar)
                        c[nbvar] = 1
                        blc = zeros(2*NS+S*(Iter-1))
                        buc = Inf * ones(2*NS+S*(Iter-1))
                    end

                    if t < T && Iter > 1
                        subi = [FixedSubi; DynamicSubi[t]]
                        subj = [FixedSubj; DynamicSubj[t]]
                        valij = [FixedValij; DynamicValij[t]]
                    else
                        subi = FixedSubi
                        subj = FixedSubj
                        valij = FixedValij
                    end

                    for m=1:NS
                        c[sum(Nb[1:(m-1)])+1:sum(Nb[1:m])] = Thermal_Costs[t,m]
                        if t > 1
                            blc[m] = gamma[t,m] * Sampled_Inflows[t,m,s] + Trial_States[t-1,m,s] - evap[t,m]
                            buc[m] = gamma[t,m] * Sampled_Inflows[t,m,s] + Trial_States[t-1,m,s] - evap[t,m]
                        else
                            blc[m] = gamma[t,m] * Sampled_Inflows[t,m,s] + x0[m] - evap[t,m]
                            buc[m] = gamma[t,m] * Sampled_Inflows[t,m,s] + x0[m] - evap[t,m]
                        end
                        blc[m+NS] = -(1-gamma[t,m]) * Sampled_Inflows[t,m,s] + Demand[t,m]
                    end
                    if t<T
                        for kp=1:(Iter-1), sp=1:S
                            buc[2*NS+(kp-1)*S+sp] = Inf
                            blc[2*NS+(kp-1)*S+sp] = Thetas[t][(kp-1)*S+sp]
                            compteur = 1
                            for m=1:NS 
                                for cp=1:StateOrder[t+1,m]
                                    if t+1-cp >= 1
                                        blc[2*NS+(kp-1)*S+sp] = blc[2*NS+(kp-1)*S+sp] + Alphas[t][(kp-1)*S+sp][compteur+cp-1] * Sampled_Inflows[t+1-cp,m,s]
                                    else
                                        blc[2*NS+(kp-1)*S+sp] = blc[2*NS+(kp-1)*S+sp] + Alphas[t][(kp-1)*S+sp][compteur+cp-1] * Recent_History_Inflows[-t+1+cp,m]
                                    end
                                end
                                compteur += StateOrder[t+1,m]
                            end
                        end
                    end
                    
                    ForwardPass = Model(HiGHS.Optimizer)
                    @variable(ForwardPass, blx[i] .<= x[i=1:nbvar] .<= bux[i])
                    if t == T || Iter == 1
                        A = sparse(subi,subj,valij,2*NS,nbvar)
                    else
                        A = sparse(subi,subj,valij,2*NS+(Iter-1)*S,nbvar)
                    end
                    @constraint(ForwardPass, blc .<= A*x .<= buc)
                    @objective(ForwardPass, Min, c'*x)

                    set_silent(ForwardPass)
                    optimize!(ForwardPass)
                    if termination_status(ForwardPass) == JuMP.INFEASIBLE
                        print("Unfeasible problem in forward pass")
                        readline()
                    end
                    sol = JuMP.value.(x)
                    if t < T
                        Trial_States[t,:,s] = sol[1+Nb_Thermal:Nb_Thermal+NS]'
                        ct = ct + sol' * c - sol[nbvar]
                    else
                        ct = ct + sol' * c
                    end
                end #for t
                if (Iter-1)*S+s <= N
                    Costs[(Iter-1)*S+s] = ct
                else
                    Costs = [Costs[2:N]; ct]
                end
            end #for s

            if Iter*S >= N
                Mean_Cost = mean(Costs)
                Sigma_Cost = sqrt(var(Costs))
                zsup = Mean_Cost + talpha * Sigma_Cost / sqrt(N)
                Zsups = [Zsups; zsup]
            end

            #Backward phase of SDDP
            for ta=1:T-1
                t=T-ta+1
                for s=1:S
                    Thetas[t-1] = [Thetas[t-1]; 0]
                    beta = zeros(NS)
                    Alphas[t-1] = [Alphas[t-1]; [zeros(taille[t])]]
                    for j=1:M
                        if t == 1
                            Inflowst = Recent_History_Inflows[1,:]
                        else
                            Sampled_Noise = Inflow_Noises[t-1,:,j]
                            Inflowst = zeros(NS)
                            for m=1:NS
                                Inflowst[m] = Tendancy_Inflows[periode[t,m],m]
                                for p=1:TabP[m][periode[t,m]]
                                    if t-p <= 1
                                        Inflowst[m] = Inflowst[m] + phi[m][periode[t,m]][p] * Recent_History_Inflows[-t+p+2,m]
                                    else
                                        Inflowst[m] = Inflowst[m] + phi[m][periode[t,m]][p] * Sampled_Inflows[t-p,m,s]
                                    end
                                end
                                Inflowst[m] = Inflowst[m] + sigmaInflows[periode[t,m],m] * Sampled_Noise[m]
                                if Inflowst[m] < 0
                                    print("Negative inflows")
                                    readline()
                                end
                            end
                        end

                        blx = []
                        bux = []
                        for m=1:NS
                            blx = [blx; vminT[t,m]]
                            bux = [bux; vmaxT[t,m]]
                        end
                        blx = [blx; lambda[t,:].*xmax[t,:]; vminH[t,:]; zeros(NS)]
                        bux = [bux; xmax[t,:]; vmaxH[t,:]; Inf*ones(NS)]
                        if t < T
                            blx = [blx; -Inf]
                            bux = [bux; Inf]
                        end

                        if t == T
                            nbvar = Nb_Thermal + 3 * NS
                            c = zeros(nbvar)
                            blc = zeros(2*NS)
                            buc = Inf * ones(2*NS)
                        else
                            nbvar = 1 + Nb_Thermal + 3 * NS
                            c = zeros(nbvar)
                            c[nbvar] = 1
                            blc = zeros(2*NS+S*Iter)
                            buc = Inf * ones(2*NS+S*Iter)
                        end

                        if t == T
                            subi = FixedSubi
                            subj = FixedSubj
                            valij = FixedValij
                        else
                            subi = [FixedSubi; DynamicSubi[t]]
                            subj = [FixedSubj; DynamicSubj[t]]
                            valij = [FixedValij; DynamicValij[t]]
                        end

                        for m=1:NS
                            c[sum(Nb[1:(m-1)])+1:sum(Nb[1:m])] = Thermal_Costs[t,m]         
                            if t > 1
                                blc[m] = gamma[t,m] * Inflowst[m] + Trial_States[t-1,m,s] - evap[t,m]
                                buc[m] = gamma[t,m] * Inflowst[m] + Trial_States[t-1,m,s] - evap[t,m]
                            else
                                blc[m] = gamma[t,m] * Inflowst[m] + x0[m] - evap[t,m]
                                buc[m] = gamma[t,m] * Inflowst[m] + x0[m] - evap[t,m]
                            end
                            blc[m+NS] = -(1-gamma[t,m]) * Inflowst[m] + Demand[t,m]
                        end

                        if t < T
                            for kp=1:Iter, sp=1:S
                                blc[2*NS+(kp-1)*S+sp] = Thetas[t][(kp-1)*S+sp]
                                compteur = 1
                                for m=1:NS
                                    for cp=1:StateOrder[t+1,m]
                                        if t+1-cp <= 1
                                            blc[2*NS+(kp-1)*S+sp] = blc[2*NS+(kp-1)*S+sp] + Alphas[t][(kp-1)*S+sp][compteur+cp-1] * Recent_History_Inflows[-t+1+cp,m]
                                        elseif cp == 1
                                            blc[2*NS+(kp-1)*S+sp] = blc[2*NS+(kp-1)*S+sp] + Alphas[t][(kp-1)*S+sp][compteur+cp-1] * Inflowst[m]
                                        else
                                            blc[2*NS+(kp-1)*S+sp] = blc[2*NS+(kp-1)*S+sp] + Alphas[t][(kp-1)*S+sp][compteur+cp-1] * Sampled_Inflows[t+1-cp,m,s]
                                        end
                                    end
                                    compteur += StateOrder[t+1,m]
                                end
                            end
                        end

                        BackwardPass = Model(HiGHS.Optimizer)
                        @variable(BackwardPass, blx[i] .<= x[i=1:nbvar] .<= bux[i])
                        if t == T
                            A = sparse(subi,subj,valij,2*NS,nbvar)
                        else
                            A = sparse(subi,subj,valij,2*NS+Iter*S,nbvar)
                        end
                        @constraint(BackwardPass, BackwardUpper, -Inf .<= A*x .<= buc)
                        @constraint(BackwardPass, BackwardLower, blc .<= A*x .<= Inf)
                        @objective(BackwardPass, Min, c'*x)

                        set_silent(BackwardPass)
                        optimize!(BackwardPass)
                        if termination_status(BackwardPass) == JuMP.INFEASIBLE
                            print("Unfeasible primal problem in backward pass")
                            readline()
                        elseif termination_status(BackwardPass) == JuMP.DUAL_INFEASIBLE
                            print("Primal infinite optimal value in backward pass")
                            readline()
                        end

                        ThetaAux = objective_value(BackwardPass)
                        dualL = [dual(i) for i in BackwardLower]
                        dualU = [dual(i) for i in BackwardUpper]
                        lambda1 = dualL[1:NS] - dualU[1:NS]
                        lambda2 = dualL[NS+1:2*NS]
                        ThetaAux -= lambda1' * Trial_States[t-1,:,s]
                        beta += Probabilities[t-1,j] * lambda1
                        compteur = 1
                        compteurbis = 1
                        AlphasAux = zeros(1,taille[t])
                        for m=1:NS
                            Aux = (lambda1[m] * gamma[t,m] - lambda2[m] * (1-gamma[t,m])) * (phi[m][periode[t,m]])
                            AlphasAux[compteur:compteur+TabP[m][periode[t,m]]-1] = Aux
                            if t < T
                                for kp=1:Iter, sp=1:S
                                    AlphasAux[compteur:(compteur+TabP[m][periode[t,m]]-1)] += Alphas[t][(kp-1)*S+sp][compteurbis] * phi[m][periode[t,m]]' * dualL[2*NS+(kp-1)*S+sp]
                                    for w=2:StateOrder[t+1,m]
                                        AlphasAux[compteur+w-2] += Alphas[t][(kp-1)*S+sp][compteurbis+w-1] * dualL[2*NS+(kp-1)*S+sp]
                                    end
                                end
                                compteurbis += StateOrder[t+1,m]
                            end
                            for Ind=1:StateOrder[t,m]
                                if t-Ind >= 1
                                    ThetaAux -= AlphasAux[compteur+Ind-1] * Sampled_Inflows[t-Ind,m,s]
                                else
                                    ThetaAux -= AlphasAux[compteur+Ind-1] * Recent_History_Inflows[-t+Ind+2,m] 
                                end
                            end
                            compteur += StateOrder[t,m]
                        end
                        Thetas[t-1][(Iter-1)*S+s] += ThetaAux * Probabilities[t-1,j]
                        Alphas[t-1][(Iter-1)*S+s] .+= Probabilities[t-1,j] * AlphasAux'
                    end #for j
                    DynamicSubi[t-1] = [DynamicSubi[t-1]; (2*NS+(Iter-1)*S+s)*ones(NS+1)]
                    DynamicSubj[t-1] = [DynamicSubj[t-1]; Nb_Thermal+3*NS+1; collect(Nb_Thermal+1:1:Nb_Thermal+NS)]
                    DynamicValij[t-1] = [DynamicValij[t-1]; 1; -beta]
                end #for s
            end #for ta
            #zinf

            t = 1
            Inflowst = Recent_History_Inflows[1,:]

            blx = []
            bux = []
            for m=1:NS
                blx=[blx; vminT[m]]
                bux=[bux; vmaxT[m]]
            end
            blx = [blx; (lambda[1,:].*xmax[1,:]); vminH[1,:]; zeros(NS)]
            bux = [bux; xmax[1,:]; vmaxH[1,:]; Inf*ones(NS)]
            blx = [blx; -Inf]
            bux = [bux; Inf]

            nbvar = 1 + Nb_Thermal + 3 * NS
            c = zeros(nbvar)
            c[nbvar] = 1
            blc = zeros(2*NS+S*Iter)
            buc = Inf * ones(2*NS+S*Iter)
            
            subi = [FixedSubi; DynamicSubi[1]]
            subj = [FixedSubj; DynamicSubj[1]]
            valij = [FixedValij; DynamicValij[1]]

            for m=1:NS
                c[sum(Nb[1:(m-1)])+1:sum(Nb[1:m])] = Thermal_Costs[m]
                
                c[Nb_Thermal+2*NS+m] = 0
                c[Nb_Thermal+NS+m] = 0
                
                blc[m] = gamma[1,m] * Inflowst[m] + x0[m] - evap[1,m]
                buc[m] = gamma[1,m] * Inflowst[m] + x0[m] - evap[1,m]
                
                blc[m+NS] =- (1-gamma[1,m]) * Inflowst[m] + Demand[1,m]
            end

            for kp=1:Iter, sp=1:S
                blc[2*NS+(kp-1)*S+sp] = Thetas[t][(kp-1)*S+sp]
                compteur = 1
                for m=1:NS
                    for cp = 1:StateOrder[t+1,m]
                        blc[2*NS+(kp-1)*S+sp] += Alphas[t][(kp-1)*S+sp][compteur+cp-1] * Recent_History_Inflows[-t+1+cp,m]
                    end
                    compteur += StateOrder[t+1,m]
                end
            end

            MSLP = Model(HiGHS.Optimizer)
            @variable(MSLP, blx[i] .<= x[i=1:nbvar] .<= bux[i])
            A = sparse(subi,subj,valij,2*NS+Iter*S,nbvar)
            @constraint(MSLP, blc .<= A*x .<= buc)
            @objective(MSLP, Min, c'*x)
            set_silent(MSLP)
            optimize!(MSLP)
            sol = JuMP.value.(x)
            if termination_status(MSLP) == JuMP.INFEASIBLE
                print("Unfeasible primal problem in backward pass")
                readline()
            elseif termination_status(MSLP) == JuMP.DUAL_INFEASIBLE
                print("Primal infinite optimal value in backward pass")
                readline()
            end

            #Update zinf
        
            zinf = objective_value(MSLP)
            Zinfs = [Zinfs; zinf]
            #zinf
            #zsup
            if Iter == Itermax
                End_Algo = false
            elseif Iter * S >= N && abs((zsup-zinf)/zsup) <= tol
                End_Algo = false
            end
            print(zsup, ",")
            print(zinf, "\n")
            Iter += 1
        end
        Iter = Iter - 1
    end
    return Time, Zsups, Zinfs, Iter, DynamicSubi, DynamicSubj, DynamicValij, Alphas, Thetas
end

