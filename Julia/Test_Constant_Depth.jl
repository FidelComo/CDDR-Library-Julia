include("init_data_hydro_sddp.jl")
include("init_data_hydro_cddr.jl")
include("sddp_hydro_interstage_dependent.jl")
include("solve_constant_depth_decision_rules.jl")
include("solve_constant_depth_decision_rules_depth_one.jl")

function Test_Constant_Depth(T,M,path,AddP)
    Times = zeros(5,2)
    Costs = zeros(5,2)
    CV = zeros(4,2)
    sigmaN = 0.2
    phinp,Nb,NS,gamma,Inflow_Noises,sigmaNoises,Probabilities,Thermal_Costs,Mus,sigmaInflows,phi,Recent_History_Inflows,Tendancy_Inflows,TabP,TabS,vminT,vmaxT,vminH,vmaxH,xmax,x0,lambda,evap,Demand = init_data_hydro_sddp(T,M,2,path,sigmaN,AddP)
    lambda = 0 * ones(T,NS)
    talpha = 1.96
    N = 100
    tol = 0.05
    Itermax = 2000
    S = 1

    #SDDP
    #Time,Zsups,Zinfs,Iter,DynamicSubi,DynamicSubj,DynamicValij,Alphas,Thetas = sddp_hydro_interstage_dependent(Nb,talpha,T,S,NS,N,gamma,M,Inflow_Noises,Probabilities,Demand,Thermal_Costs,vminT,vmaxT,vminH,vmaxH,xmax,x0,lambda,evap,phi,TabS,TabP,tol,Tendancy_Inflows,Recent_History_Inflows,sigmaInflows,2000)
    #Times[1,1] = Time
    #Costs[1,1] = Zinfs[length(Zinfs)]
    #Costs[1,2] = Zsups[length(Zsups)]

    #CDDR depth=1
    path = ""
    depth = 1                                           
    subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = init_data_hydro_cddr(depth,T,M,path,AddP,Inflow_Noises)
    #For initialization with b_t of type depth>=1
    #[subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities]=init_data_cddr(2,T,M,1,path,sigmasN,AddP,Inflow_Noises,2)

    sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules_depth_one(subi_a,subj_a,valij_a,betas,cost,probabilities,ds,T,ps,qs)
    Costs[2,1] = opt_value
    CV[1,1] = nvars
    CV[1,2] = counter
    Times[2,1] = time1
    Times[2,2] = time2

    #CDDR depth=2
    depth = 2
    subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = init_data_hydro_cddr(depth,T,M,path,AddP,Inflow_Noises)

    sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
    Costs[depth+1,1] = opt_value
    CV[depth,1] = nvars
    CV[depth,2] = counter
    Times[depth+1,1] = time1
    Times[depth+1,2] = time2

    #CDDR depth=3
    depth = 3
    subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = init_data_hydro_cddr(depth,T,M,path,AddP,Inflow_Noises)

    sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
    Costs[depth+1,1] = opt_value
    CV[depth,1] = nvars
    CV[depth,2] = counter
    Times[depth+1,1] = time1
    Times[depth+1,2] = time2

    #CDDR depth=4
    depth = 4
    subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = init_data_hydro_cddr(depth,T,M,path,AddP,Inflow_Noises)

    sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
    Costs[depth+1,1] = opt_value
    CV[depth,1] = nvars
    CV[depth,2] = counter
    Times[depth+1,1] = time1
    Times[depth+1,2] = time2

    display(Costs)
    display(CV)
    display(Times)

    return Costs,Times,CV
end