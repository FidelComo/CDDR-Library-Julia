using JuMP
using HiGHS

###############################################################################
#Inputs
###############################################################################
#Matrix A^{t tau} is given in sparse format by 
#subi_a[t][tau], subj_a[t][tau], valij_a[t][tau]. 
#beta[t][s][ell,i] is beta^{t i}_{s scenario(ell)} of the paper where 
#scenario(ell) is the scenario numbered ell.
#cost[t] is the cost vector for stage t
#probabilities[t-1,k] is the probability of kth realization for stage t=2,..,T.
#ds[t] is the number of realizations for stage t
#T is the number of stages
#A^{t tau} has size ps[t]*qs[t] 
###############################################################################
#Outputs
###############################################################################
#sol: optimal solution 
#opt_value: optimal value
#out=0: optimal solution found; out=1: infeasible; out=2: primal infinite
#opt value
#time1: time to define the optimization problem
#time2: time to solve the optimization problem
#nvars: number of variables of the problem
#nconst: number of constraints of the problem
###############################################################################

function solve_constant_depth_decision_rules_depth_one(subi_a,subj_a,valij_a,betas,cost,probabilities,ds,T,ps,qs)
    time1 = @elapsed begin

        CDDR1 = Model(HiGHS.Optimizer)

        #Variables
        @variable(CDDR1, v[t = 1:T, i = 1:ps[t], s = 1:t, k = 1:((ds[s]-1)*(t!=1)+1)])
        @variable(CDDR1, u[t = 1:T, i = 1:qs[t], s = 1:t, k = 1:((ds[s]-1)*(t!=1)+1)])
        @variable(CDDR1, z[t = 2:T, i = 1:ps[t], s = 2:t])
        
        
        #Constraints
        @constraint(CDDR1, v[1,:,1,1] .<= 0)

        nconst = 0
        for t=1:T, s=1:t, i=1:ps[t], eta=1:ds[s]
            var = []
            val = []
            for tau=s:t, k in findall(subi_a[t][tau] .== i)
                append!(var, u[tau,subj_a[t][tau][k],s,eta])
                append!(val, valij_a[t][tau][k])
            end
            @constraint(CDDR1, v[t,i,s,eta] - val'*var == -betas[t][s][eta,i])
            nconst += 1
        end

        for t=2:T, i=1:ps[t]
            s = 2:t
            @constraint(CDDR1, v[t,i,1,1] + sum(z[t,i,s]) <= 0)
            nconst += 1
        end

        for t=2:T, i=1:ps[t], s=2:t, k=1:ds[s]
            @constraint(CDDR1, z[t,i,s] >= v[t,i,s,k])
            nconst += 1
        end

        #Objective
        c=zero(u)
        for t=1:T, s=1:t, i=1:qs[t], k=1:ds[s]
            if s >= 2
                c[t,i,s,k] = cost[t][i]*probabilities[s-1,k]
            else
                c[t,i,s,k] = cost[t][i]
            end
        end

        @objective(CDDR1, Min, sum(c.*u))
    end

    time2 = @elapsed begin
        #Solution
        set_silent(CDDR1)
        optimize!(CDDR1)

        sol = [JuMP.value.(u), JuMP.value.(v), JuMP.value.(z)]

        if termination_status(CDDR1) == JuMP.INFEASIBLE
            print("Infeasible problem")
            out=1
            opt_value=Inf
        elseif termination_status(CDDR1) == JuMP.DUAL_INFEASIBLE
            print("Primal infinite optimal value")
            out=2
            opt_value=-Inf
        else
            out=0
            opt_value = objective_value(CDDR1)
        end
        nvar = num_variables(CDDR1)
    end

    return sol, opt_value, out, time1, time2, nvar, nconst
end