using MATLAB
using JSON

mat"addpath('.MATLAB')"
mat"addpath('./MatlabCode')"
mat"addpath('./MatlabCode/MATLAB')"
T = 20
M = 10.0
AddP = 5
path = ""

Times = zeros(5,2)
Costs = zeros(5,2)
CV = zeros(4,2)
sigmaN = 0.2
phinp,Nb,NS,gamma,Inflow_Noises,sigmaNoises,Probabilities,Thermal_Costs,Mus,sigmaInflows,phi,Recent_History_Inflows,Tendancy_Inflows,TabP,TabS,vminT,vmaxT,vminH,vmaxH,xmax,x0,lambda,evap,Demand = mxcall(:init_data_hydro_sddp,24,T,M,2,path,sigmaN,AddP)
NS = Int(NS)
lambda = 0 * ones(T,NS)
talpha = 1.96
N = 100
tol = 0.05
Itermax = 2000
S=1

#CDDR depth=1
path = ""
depth = 1                                           
subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = mxcall(:init_data_hydro_cddr,9,depth,T,M,path,AddP,Inflow_Noises)
#For initialization with b_t of type depth>=1
#[subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities]=init_data_cddr(2,T,M,1,path,sigmasN,AddP,Inflow_Noises,2)

subi_a = JSON.parse(subi_a)
subj_a = JSON.parse(subj_a)
valij_a = JSON.parse(valij_a)

ps = Int.(ps)
qs = Int.(qs)
ds = Int.(ds)

sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules_depth_one(subi_a,subj_a,valij_a,betas,cost,probabilities,ds,T,ps,qs)
Costs[2,1] = opt_value
CV[1,1] = nvars
CV[1,2] = counter
Times[2,1] = time1
Times[2,2] = time2
print(counter,"\n")

#CDDR depth=2
depth = 2
subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = mxcall(:init_data_hydro_cddr,9,depth,T,M,path,AddP,Inflow_Noises)

subi_a = JSON.parse(subi_a)
subj_a = JSON.parse(subj_a)
valij_a = JSON.parse(valij_a)

ps = Int.(ps)
qs = Int.(qs)
ds = Int.(ds)

sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
Costs[depth+1,1] = opt_value
CV[depth,1] = nvars
CV[depth,2] = counter
Times[depth+1,1] = time1
Times[depth+1,2] = time2

#CDDR depth=3
depth = 3
subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = mxcall(:init_data_hydro_cddr,9,depth,T,M,path,AddP,Inflow_Noises)

subi_a = JSON.parse(subi_a)
subj_a = JSON.parse(subj_a)
valij_a = JSON.parse(valij_a)

ps = Int.(ps)
qs = Int.(qs)
ds = Int.(ds)

sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
Costs[depth+1,1] = opt_value
CV[depth,1] = nvars
CV[depth,2] = counter
Times[depth+1,1] = time1
Times[depth+1,2] = time2

#CDDR depth=4
depth = 4
subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities = mxcall(:init_data_hydro_cddr,9,depth,T,M,path,AddP,Inflow_Noises)

subi_a = JSON.parse(subi_a)
subj_a = JSON.parse(subj_a)
valij_a = JSON.parse(valij_a)

ps = Int.(ps)
qs = Int.(qs)
ds = Int.(ds)

sol,opt_value,out,time1,time2,nvars,counter = solve_constant_depth_decision_rules(subi_a,subj_a,valij_a,betas,cost,ds,T,ps,qs,probabilities,depth)
Costs[depth+1,1] = opt_value
CV[depth,1] = nvars
CV[depth,2] = counter
Times[depth+1,1] = time1
Times[depth+1,2] = time2