using DelimitedFiles

include("compute_decomposition_periodic.jl")

###########################################################################################
#Inputs
###########################################################################################
#depth: depth of the decision rules
###########################################################################################
#T: number of stages
###########################################################################################
#M: number of realizations for each stage
###########################################################################################
#option: if equal to one then s.d. of noise for stage t=1,..,12, and subsystem m
#        is estimated s.d. of noises using real data otherwise
#        sigmasN is used. 
###########################################################################################


###########################################################################################
#The outputs are the inputs of the functions computing the constant depth
#decision rules
###########################################################################################

function init_data_hydro_cddr(depth,T,M,path,AddP,Inflow_Noises)

    ###########################################################################################
    #NS is the number of subsystems
    ###########################################################################################

    NS = 4

    ###########################################################################################
    #Nb(i) is the number of thermal plants in subsystem i
    ###########################################################################################

    Nb = 2 * ones(Int, NS)

    ###########################################################################################
    #gamma(t,i) is the proportion of inflows coming to reservoir i at t.
    ###########################################################################################

    gamma = ones(T) * [0.8567 0.9013 0.9756 1]

    #Demand = 1.25 * ones(T) * [31055, 8297, 7103, 3367]
    Demand = ones(T) * [31055 8297 7103 3367]

    ###########################################################################################
    #Inflow_Noises(t,m,j) is realization j of noise for t+1 and subsytem m, t=1,..,T-1.
    ###########################################################################################

    ###########################################################################################
    #vmaxH: vmaxH(t,i) is the maximal hydraulic available power at t for subsystem i.
    ###########################################################################################

    vmaxH = ones(T) * [40943.0 10292.5 9406.0 6915.7]

    ###########################################################################################
    #vminH: vminH(t,i) is the minimal hydro production at t i
    ###########################################################################################

    vminH = zeros(T,NS)

    ###########################################################################################
    #evap:evap(t,i) is the evaporation for time step t and subsytem i.
    ###########################################################################################

    evap = zeros(T,NS)

    Thermal_Costs = [[] for i=1:T, j=1:NS]
    vmaxT = [[] for i=1:T, j=1:NS]
    vminT = [[] for i=1:T, j=1:NS]

    for t=1:T
        Thermal_Costs[t,1] = [17.48;12.61;937;6.27;100.4;77.46;432.7;150;10.5;42.6;74.4;108;180;395.87;180;254.14;288.51;288.51;589.31;97.15;124.77;110.48;1047.38;197.85;4170.44]
        vmaxT[t,1] = [417.37;1212.6;26.53;342.28;26.32;0;108.72;79.39;384;93.12;186.24;91.28;0;0;0;0;115.64;343.93;152.12;0;409.02;185.32;8;152.03;100000000]
        vminT[t,1] = zeros(25)
        Thermal_Costs[t,2] = [546.4;219;110.48;191.08;193.25;200.17;160.03;155;116.1;596;116;116;249;92.49;4170.44]
        vmaxT[t,2] = [49.26;440.32;0;54.37;7.87;48.50;105.14;219.19;164.66;20.06;76.65;180.53;11.44;204.7;100000000]
        vminT[t,2] = zeros(15)
        Thermal_Costs[t,3] = [161.51;71.29;80.02;87.12;82.72;65;4170.44]
        vmaxT[t,3] = [307.18;51.17;327.26;102.29;211.48;83.94;100000000]
        vminT[t,3] = zeros(7)
        Thermal_Costs[t,4] = [4170.44]
        vmaxT[t,4] = [100000000]
        vminT[t,4] = [0]
    end

    #Update Thermal Costs to have one equivalent TP per subsystem
    for t=1:T
        for n=1:NS
            Thermal_Costs[t,n] = [maximum(Thermal_Costs[t,n][1:Nb[n]-1]); 4170.44]
            vmaxT[t,n] = [0.6*sum(Demand[:,n])/T; Inf]
            vminT[t,n] = [0; 0]
        end
    end

    ###########################################################################################
    #Mus(t,i) is empirical mean of inflow for subsystem i month t
    ###########################################################################################
    #sigmaInflows(t,i) is empirical s.d. of inflow for subsytem i month t
    ###########################################################################################
    #phi: Normalized coefficients phi of PAR models phi{1,i}[t] are
    #coefficients for subsystem i and month t.
    ###########################################################################################

    Mus = zeros(12,NS)
    sigmaInflows = zeros(12,NS)
    pathaux = string(path,"Historique_Apports_Sud_Est.txt")
    Historique_SE = readdlm(pathaux)
    for i=1:12
        Mus[i,1] = mean(Historique_SE[:,i]);
        sigmaInflows[i,1] = sqrt(var(Historique_SE[:,i]))
    end
    pathaux = string(path,"Historique_Apports_Sud.txt")
    Historique_S = readdlm(pathaux)
    for i=1:12
        Mus[i,2] = mean(Historique_S[:,i])
        sigmaInflows[i,2] = sqrt(var(Historique_S[:,i]))
    end
    pathaux = string(path,"Historique_Apports_Nord_Est.txt")
    Historique_NE = readdlm(pathaux)
    for i=1:12
        Mus[i,3] = mean(Historique_NE[:,i])
        sigmaInflows[i,3] = sqrt(var(Historique_NE[:,i]))
    end
    pathaux = string(path,"Historique_Apports_Nord.txt")
    Historique_Nord = readdlm(pathaux)
    for i=1:12
        Mus[i,4] = mean(Historique_Nord[:,i])
        sigmaInflows[i,4] = sqrt(var(Historique_Nord[:,i]))
    end

    ###########################################################################################
    #Recent_History_Inflows(-t+2,:) inflows for time t (t=1,0,-1,-2,...)
    ###########################################################################################

    #Recent_History_Inflows(-t+2,:) inflows for time t (t=1,0,-1,-2,...)
    Recent_History_Inflows = zeros(13,4)
    Recent_History_Inflows[1,1] = Historique_SE[end,12]
    for t=1:12
        Recent_History_Inflows[14-t,1] = Historique_SE[end,t]
    end

    Recent_History_Inflows[1,2] = Historique_S[end,12]
    for t=1:12
        Recent_History_Inflows[14-t,2] = Historique_S[end,t]
    end

    Recent_History_Inflows[1,3] = Historique_NE[end,12]
    for t=1:12
        Recent_History_Inflows[14-t,3] = Historique_NE[end,t]
    end

    Recent_History_Inflows[1,4] = Historique_Nord[end,12]
    for t=1:12
        Recent_History_Inflows[14-t,4] = Historique_Nord[end,t]
    end

    ###########################################################################################
    #TabP: TabP is a cell array of size NS. TabP{1,i} is the list of PAR model
    #i periods.
    ###########################################################################################
    #TabS:TabS(i) is the period of PAR model i.
    ###########################################################################################


    phinp = [[] for _=1:NS]
    TabS = [12,12,12,12]
    TabP = [Int[] for _=1:NS]
    TabP[1] = [1,1,1,2,3,1,3,1,1,3,1,4]
    TabP[2] = [1,1,1,1,1,1,4,1,1,1,1,1]
    TabP[3] = [5,2,1,1,1,1,2,1,3,3,2,5]
    TabP[4] = [1,4,1,1,2,1,3,2,5,3,5,1]

    for m=1:NS
        phinp[m] = [[] for _=1:12]
    end

    phinp[1][1] = [0.593]
    phinp[1][2] = [0.569]
    phinp[1][3] = [0.671]
    phinp[1][4] = [0.601,0.274]
    phinp[1][5] = [0.617,-0.0484,0.332]
    phinp[1][6] = [0.824]
    phinp[1][7] = [0.699,0.0108,0.284]
    phinp[1][8] = [0.838]
    phinp[1][9] = [0.843]
    phinp[1][10] = [0.321,0.130,0.336]
    phinp[1][11] = [0.730]
    phinp[1][12] = [0.708,-0.203,-0.0600,0.387]

    phinp[2][1] = [0.399]
    phinp[2][2] = [0.568]
    phinp[2][3] = [0.665]
    phinp[2][4] = [0.508]
    phinp[2][5] = [0.468]
    phinp[2][6] = [0.642]
    phinp[2][7] = [0.427,0.350,-0.239,0.361]
    phinp[2][8] = [0.457]
    phinp[2][9] = [0.569]
    phinp[2][10] = [0.474]
    phinp[2][11] = [0.54]
    phinp[2][12] = [0.587]

    phinp[3][1] = [0.709,-0.183,0.134,0.281,-0.233]
    phinp[3][2] = [0.775,-0.344]
    phinp[3][3] = [0.777]
    phinp[3][4] = [0.707]
    phinp[3][5] = [0.827]
    phinp[3][6] = [0.948]
    phinp[3][7] = [1.2,-0.248]
    phinp[3][8] = [0.979]
    phinp[3][9] = [1.12,0.0603,-0.247]
    phinp[3][10] = [0.760,0.679,-0.604]
    phinp[3][11] = [0.971,-0.336]
    phinp[3][12] = [0.708,-0.0534,-0.0342,-0.614,0.606]

    phinp[4][1] = [0.739]
    phinp[4][2] = [0.859,-0.449,-0.0423,0.286]
    phinp[4][3] = [0.778]
    phinp[4][4] = [0.774]
    phinp[4][5] = [0.995,-0.227]
    phinp[4][6] = [0.895]
    phinp[4][7] = [1.12,-0.528,0.339]
    phinp[4][8] = [1.18,-0.228]
    phinp[4][9] = [1.27,-0.754,0.132,-0.0252,0.321]
    phinp[4][10] = [0.476,0.672,-0.335]
    phinp[4][11] = [0.639,-0.295,0.822,-0.246,-0.251]
    phinp[4][12] = [0.708]

    for i=1:4 
        for j=1:12
            TabP[i][j] = TabP[i][j] + AddP
            phinp[i][j]=[phinp[i][j]; 0.01 * ones(AddP)]
        end
    end

    phi = [[] for _=1:NS]

    for m=1:NS
        phi[m] = [[] for _=1:TabS[m]]
        for s=1:TabS[m]
            phi[m][s] = zeros(TabP[m][s])
            for j=1:TabP[m][s]
                if s-j >= 1
                    phi[m][s][j] = phinp[m][s][j] * sigmaInflows[s,m] / sigmaInflows[s-j,m]
                else
                    phi[m][s][j] = phinp[m][s][j] * sigmaInflows[s,m] / sigmaInflows[s-j+TabS[m],m]
                end
            end
        end
    end


    #Computation of sigmasNoises in case option=1
    H_SE = []
    H_S = []
    H_NE = []
    H_Nord = []
    for i in axes(Historique_SE,1)
        H_SE=[H_SE; [Historique_SE[i,:]]]
    end
    for i in axes(Historique_S,1)
        H_S = [H_S; [Historique_S[i,:]]]
    end
    for i in axes(Historique_NE,1)
        H_NE=[H_NE; [Historique_NE[i,:]]]
    end
    for i in axes(Historique_Nord,1)
        H_Nord=[H_Nord; [Historique_Nord[i,:]]]
    end

    Echantillons = [H_SE,H_S,H_NE,H_Nord]

    # sigmasNoises=zeros(12,4);
    # for t=1:12
    #     %Computation of sigmasNoises(t,m)
    #     %Computation of sample of noise for month t subsystem m
    #     for m=1:NS
    # Sample=[];
    #         for k=1:((size(Echantillons,1)/12)-1)
    #             noisevalue=(Echantillons(t+12*k,m)-Mus(t,m))/sigmaInflows(t,m);
    #             for j=1:TabP{1,m}(t)
    #                 if (t-j>=1)
    #                     aux=(Echantillons(t+12*k-j,m)-Mus(t-j,m))/sigmaInflows(t-j,m);
    #                 else
    #                     aux=(Echantillons(t+12*k-j,m)-Mus(t-j+TabS(m),m))/sigmaInflows(t-j+TabS(m),m);
    #                 end
    #                 noisevalue=noisevalue-phinp{1,m}{1,t}(j)*aux;
    #             end
    #             Sample=[Sample,noisevalue];
    #         end
    #         if (option==1)
    #             sigmaNoises(t,m)=sqrt(var(Sample));
    #         else
    #             sigmaNoises(t,m)=sigmasN;
    #         end
    #     end
    # end

    probabilities = zeros(T-1,M)
    for t=1:T-1
        for j=1:M
            probabilities[t,j] = 1/M
        end
    end

    ###########################################################################################
    #lambda: xmin(t,m) is lambda(t,m)xmax(t,m)
    ###########################################################################################


    #= lambda = zeros(T,NS)
    auxiliaire = [33;35;38;37;38;36;32;27;20;14;10;10]/100
    i = 1
    lambda[(i-1)*12+1:i*12,1] = auxiliaire
    auxiliaire = [19;20;19;16;13;13;13;13;13;13;13;13]/100
    lambda[(i-1)*12+1:i*12,2] = auxiliaire
    auxiliaire = [33;36;39;40;37;35;30;24;18;13;10;10]/100
    lambda[(i-1)*12+1:i*12,3] = auxiliaire

    lambda = 0.2 * ones(T,NS) =#

    ###########################################################################################
    #probabilities(t,j): probability of realization Inflow_Noises(t,:,j) for
    #                    t+1,t=1,...,T-1.
    ###########################################################################################

    j = 1
    pmaxS = zeros(Int,NS)
    for m=1:NS
        reste = (j+1) % TabS[m]
        if reste == 0
            pmaxS[m] = TabP[m][TabS[m]] - 1
        else
            pmaxS[m] = TabP[m][reste] - 1
        end
        for s=2:TabS[m]
            reste = (j+s) % TabS[m]
            if reste == 0
                val = TabP[m][TabS[m]] - s
            else
                val = TabP[m][reste] - s
            end
            if val > pmaxS[m] 
                pmaxS[m] = val 
            end
        end
    end

    ###########################################################################################
    #xmax: xmax(t,i) is the maximal level of subsystem i reservoir at t.
    #x0:x0(i) is the initial level subsystem i reservoir.
    #Demand: a (T,NS) matrix, demande(t,i) is the deterministic demand for
    ###########################################################################################

    #Compute reduced demand
    x0 = zeros(NS)
    for m=1:NS
        x0[m] = 0.3 * sum(Demand[:,m])
    end
    xmax = 2 * ones(T,1) * x0'

    #x0=[197263.3;18374.0;51806.1;12415.2]

    ds = [1; M*ones(Int,T-1)]
    ps = 9*NS*ones(Int,T)
    qs = 4*NS*ones(Int,T)

    subi_a = [[] for _=1:T]
    subj_a = [[] for _=1:T]
    valij_a = [[] for _=1:T]
    for t=1:T
        subi_a[t] = [[] for _=1:t]
        subj_a[t] = [[] for _=1:t]
        valij_a[t] = [[] for _=1:t]
        
        row = 0
        col = 0
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; ones(NS)]
        
        col += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; ones(NS)]
        
        row += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        col += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        col += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        row += NS
        col = 0
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; ones(NS)]
        
        row += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        row += NS
        col += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; ones(NS)]
        
        row += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        row += NS
        col += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; ones(NS)]
        
        row += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        row += NS
        col += NS
        subi_a[t][t] = [subi_a[t][t]; collect(row+1:1:row+NS)]
        subj_a[t][t] = [subj_a[t][t]; collect(col+1:1:col+NS)]
        valij_a[t][t] = [valij_a[t][t]; -ones(NS)]
        
        if t >= 2
            subi_a[t][t-1] = collect(1:1:NS)
            subj_a[t][t-1] = collect(1:1:NS)
            valij_a[t][t-1] = -ones(NS)
        end
    end

    ds = [ones(Int,depth-1,1);ds]
    productsu = zeros(Int,T,depth)
    for s=1:T
        productsu[s,depth] = ds[s+depth-1]
        for i=1:(depth-1)
            productsu[s,depth-i] = productsu[s,depth-i+1] * ds[s+depth-i-1]
        end
    end
    ds = [1; M*ones(Int,T-1)]

    alpha,beta,theta = compute_decomposition_periodic(Mus,sigmaInflows,phinp,NS,TabS,TabP,1,T,pmaxS)

    betas=[[] for _=1:T]
    for t=1:T
        betas[t]=[[] for _=1:t]
        if t == 1
            betas[1][1] = [(x0+gamma[1,:].*Recent_History_Inflows[1,:])' ((1 .- gamma[1,:]).*Recent_History_Inflows[1,:].-Demand[1,:])' xmax[1,:]' zeros(NS)' vmaxH[1,:]' zeros(NS)' vmaxT[1][1]' vmaxT[2][1]' vmaxT[3][1]' vmaxT[4][1]' zeros(2*NS)']          
        else
            aux = zeros(NS)
            for m=1:NS
                aux[m] = theta[t-1,m]
                for ell=0:pmaxS[m]
                    aux[m] = aux[m]+alpha[m][t-1,ell+1] * Recent_History_Inflows[ell+1,m]
                end
            end
            betas[t][1] = [(aux.*gamma[t,:])' ((1 .- gamma[t,:]).*aux-Demand[t,:])' xmax[t,:]' zeros(NS)' vmaxH[t,:]' zeros(NS)' vmaxT[t,1][1]' vmaxT[t,2][1]' vmaxT[t,3][1]' vmaxT[t,4][1]' zeros(2*NS)']
            for s=2:t
                betas[t][s] = zeros(productsu[s,1],9*NS)
                nb_blocks = Int(productsu[s,1]/M)
                for i=1:nb_blocks
                    for j=1:M
                        aux = zeros(NS)
                        for m=1:NS
                            aux[m] = beta[m,t-1][s-1] * Inflow_Noises[s-1,m,j]
                        end
                        betas[t][s][(i-1)*M+j,:] = [gamma[t,:].*aux; (1 .- gamma[t,:]).*aux; zeros(7*NS)]
                    end
                end
            end
        end
    end

    cost=[[] for _=1:T]
    for t=1:T
        cost[t] = [zeros(2*NS); Thermal_Costs[t,1][1]; Thermal_Costs[t,2][1]; Thermal_Costs[t,3][1]; Thermal_Costs[t,4][1]; Thermal_Costs[t,1][2]; Thermal_Costs[t,2][2]; Thermal_Costs[t,3][2]; Thermal_Costs[t,4][2]]
    end

    return subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities
end