using JuMP
using HiGHS

include("index_to_scenario.jl")
include("scenario_to_index.jl")

###############################################################################
#Inputs
###############################################################################
#Matrix A^{t tau} is given in sparse format by 
#subi_a{1,t}{1,tau}, subj_a{1,t}{1,tau}, valij_a{1,t}{1,tau}. 
#beta{1,t}{1,s}(ell,i) is beta^{t i}_{s scenario(ell)} of the paper where 
#scenario(ell) is the scenario numbered ell.
#cost{1,t} is the cost vector for stage t
#probabilities(t-1,k) is the probability of kth realization for stage t=2,..,T.
#ds(t) is the number of realizations for stage t
#T is the number of stages
#A^{t tau} has size ps(t)*qs(t) 
###############################################################################
#Outputs
###############################################################################
#sol: optimal solution 
#opt_value: optimal value
#out=0: optimal solution found; out=1: infeasible; out=2: primal infinite
#opt value
#time1: time to store the data for the optimization problem
#time2: time to solve the optimization problem
#nvars: number of variables of the problem
#counter: number of constraints of the problem
###############################################################################

function solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
    time1 = @elapsed begin
        
        ds=[ones(depth-1,1);ds]

        productsu=zeros(Int,T,depth)
        productsz=zeros(Int,T,depth-1)
        for s=1:T
            productsu[s,depth]=ds[s+depth-1]
            productsz[s,depth-1]=ds[s+depth-2]
            for i=1:(depth-2)
                productsu[s,depth-i]=productsu[s,depth-i+1]*ds[s+depth-i-1]
                productsz[s,depth-1-i]=productsz[s,depth-i]*ds[s+depth-2-i]
            end
            productsu[s,1]=productsu[s,2]*ds[s]
        end

        CDDR = Model(HiGHS.Optimizer)

        #Variables
        @variable(CDDR, v[t = 1:T, s = 1:t, ell = 1:productsu[s,1], i = 1:ps[t]])
        @variable(CDDR, u[t = 1:T, s = 1:t, ell = 1:productsu[s,1], i = 1:qs[t]])
        @variable(CDDR, z[t = 2:T, s = 2:t, ell = 1:productsz[s,1], i = 1:ps[t]])

        #Constraints
        @constraint(CDDR, v[1,1,1,:] .<= 0)
        
        nconst = 0
        for t=1:T, s=1:t, i=1:ps[t], ell=1:productsu[s,1]
            var = []
            val = []
            for tau=s:t, k in findall(subi_a[t][tau] .== i)
                append!(var, u[tau,s,ell,subj_a[t][tau][k]])
                append!(val, valij_a[t][tau][k])
            end
            @constraint(CDDR, v[t,s,ell,i] - val'*var == -betas[t][s][ell,i])
            nconst += 1
        end
        
        for t=2:T, i=1:ps[t], ell=1:productsu[t,1]
            etas = index_to_scenario(ell, productsu[t,:], depth)
            ellz = scenario_to_index(etas[1:depth-1],productsz[t,:],depth-1)
            
            @constraint(CDDR, z[t,t,ellz,i] >= v[t,t,ell,i])
            nconst += 1
        end
        
        for t=3:T, s=3:t, i=1:ps[t], ell=1:productsu[s-1,1]
            etas = index_to_scenario(ell,productsu[s-1,:],depth)
            etasleft = etas[1:depth-1]
            etasright = etas[2:depth]
                
            ellzl = scenario_to_index(etasleft,productsz[s-1,:],depth-1)
            ellzr = scenario_to_index(etasright,productsz[s,:],depth-1)
                
            @constraint(CDDR, z[t,s-1,ellzl,i] - z[t,s,ellzr,i] >= v[t,s-1,ell,i])
            
            nconst += 1
        end
        
        for t=2:T, i=1:ps[t]
            @constraint(CDDR, z[t,2,1,i] + v[t,1,1,i] <= 0)
            nconst += 1
        end
        
        #Objective
        c = zero(u)
        for t=1:T, s=1:t, ell=1:productsu[s,1]
            etas = index_to_scenario(ell,productsu[s,:],depth)
            aux = 1
            for k=max(1,depth-s+2):depth
                aux *= probabilities[s-depth+k-1,etas[k]]
            end
            for i=1:qs[t]
                c[t,s,ell,i] = aux * cost[t][i]
            end
        end

        @objective(CDDR, Min, sum(c.*u))
    end
    
    time2 = @elapsed begin
        #Solution
        set_silent(CDDR)
        optimize!(CDDR)

        sol = [JuMP.value.(u), JuMP.value.(v), JuMP.value.(z)]

        if termination_status(CDDR) == JuMP.INFEASIBLE
            print("Infeasible problem")
            out=1
            opt_value=Inf
        elseif termination_status(CDDR) == JuMP.DUAL_INFEASIBLE
            print("Primal infinite optimal value")
            out=2
            opt_value=-Inf
        else
            out=0
            opt_value = objective_value(CDDR)
        end
        nvar = num_variables(CDDR)
    end

    return sol, opt_value, out, time1, time2, nvar, nconst
end