 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrix A^{t tau} is given in sparse format by 
%subi_a{1,t}{1,tau}, subj_a{1,t}{1,tau}, valij_a{1,t}{1,tau}. 
%beta{1,t}{1,s}(ell,i) is beta^{t i}_{s scenario(ell)} of the paper where 
%scenario(ell) is the scenario numbered ell.
%cost{1,t} is the cost vector for stage t
%probabilities(t-1,k) is the probability of kth realization for stage t=2,..,T.
%ds(t) is the number of realizations for stage t
%T is the number of stages
%A^{t tau} has size ps(t)*qs(t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sol: optimal solution 
%opt_value: optimal value
%out=0: optimal solution found; out=1: infeasible; out=2: primal infinite
%opt value
%time1: time to store the data for the optimization problem
%time2: time to solve the optimization problem
%nvars: number of variables of the problem
%counter: number of constraints of the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sol,opt_value,out,time1,time2,nvars,counter]=solve_constant_depth_decision_rules_depth_one(subi_a,subj_a,valij_a,betas,cost,probabilities,ds,T,ps,qs)

Iv=cell(1,T);
Iu=cell(1,T);
Iz=cell(1,T);

tic

%Index of v^{t i}_{s, k} is Iv{1,t}{1,i}{1,s}(k)

Iv{1,1}=cell(1,ps(1));
for i=1:ps(1)
    Iv{1,1}{1,i}=cell(1,1);
    Iv{1,1}{1,i}{1,1}(1)=i;
end

nvars=ps(1)+1;
for t=2:T
    Iv{1,t}=cell(1,ps(t));
    for i=1:ps(t)
        Iv{1,t}{1,i}=cell(1,t);
        for s=1:t
            Iv{1,t}{1,i}{1,s}=zeros(ds(s),1);
            for k=1:ds(s)
                Iv{1,t}{1,i}{1,s}(k)=nvars;
                nvars=nvars+1;
            end
        end
    end
end

%Index of u^{t i}_{s, k} is Iu{1,t}{1,i}{1,s}(k)

Iu{1,1}=cell(1,qs(1));
for i=1:qs(1)
    Iu{1,1}{1,i}=cell(1,1);
    Iu{1,1}{1,i}{1,1}(1)=nvars+i-1;
end

nvars=nvars+qs(1); 

for t=2:T
    Iu{1,t}=cell(1,qs(t));
    for i=1:qs(t)
        Iu{1,t}{1,i}=cell(1,t);
        for s=1:t
            Iu{1,t}{1,i}{1,s}=zeros(ds(s),1);
            for k=1:ds(s)
                Iu{1,t}{1,i}{1,s}(k)=nvars;
                nvars=nvars+1;
            end
        end
    end
end

%Index of z^{t i}_{s} is Iz{1,t}{1,i}(s)(k)

for t=2:T
    Iz{1,t}=cell(1,ps(t));
    for i=1:ps(t)
        Iz{1,t}{1,i}=zeros(1,s);
        for s=2:t
            Iz{1,t}{1,i}(s)=nvars;
            nvars=nvars+1;
        end
    end
end

nvars=nvars-1;

subi=[];
subj=[];
valij=[];
blx=-inf*ones(nvars,1);
bux=inf*ones(nvars,1);
blc=[];
buc=[];

bux(1:ps(1))=zeros(ps(1),1);

counter=1;
for t=1:T
    for s=1:t
        for i=1:ps(t)
            for eta=1:ds(s)
                blc=[blc;-betas{1,t}{1,s}(eta,i)];
                buc=[buc;-betas{1,t}{1,s}(eta,i)];
                subi=[subi,counter];
                subj=[subj,Iv{1,t}{1,i}{1,s}(eta)];
                valij=[valij,1];
                for tau=s:t
                    [index]=find(subi_a{1,t}{1,tau}==i);   
                    for k=1:length(index)
                        subi=[subi,counter];
                        subj=[subj,Iu{1,tau}{1,subj_a{1,t}{1,tau}(index(k))}{1,s}(eta)];
                        valij=[valij,-valij_a{1,t}{1,tau}(index(k))];
                    end
                end
                counter=counter+1;
            end
        end
    end
end

for t=2:T
    for i=1:ps(t)
        subi=[subi,counter];
        subj=[subj,Iv{1,t}{1,i}{1,1}(1)];
        valij=[valij,1];
        for s=2:t
            subi=[subi,counter];
            subj=[subj,Iz{1,t}{1,i}(s)];
            valij=[valij,1];
        end
        buc=[buc;0];
        blc=[blc;-inf];
        counter=counter+1;
    end
end

for t=2:T
    for i=1:ps(t)
        for s=2:t
            for k=1:ds(s)
                subi=[subi,counter];
                subj=[subj,Iz{1,t}{1,i}(s)];
                valij=[valij,1];
                subi=[subi,counter];
                subj=[subj,Iv{1,t}{1,i}{1,s}(k)];
                valij=[valij,-1];
                blc=[blc;0];
                buc=[buc;inf];
                counter=counter+1;     
            end
        end
    end
end
time1=toc;

obj=zeros(nvars,1);
for t=1:T
    for s=1:t
        for i=1:qs(t)
            for k=1:ds(s)
                if (s>=2)
                    obj(Iu{1,t}{1,i}{1,s}(k))=cost{1,t}(i)*probabilities(s-1,k);
                else
                    obj(Iu{1,t}{1,i}{1,s}(k))=cost{1,t}(i);
                end
            end
        end
    end
end
prob.c=obj;

prob.blx=blx;
prob.bux=bux;
prob.blc=blc;
prob.buc=buc;

tic;
prob.a = sparse(subi,subj,valij,counter-1,nvars);
[~,res]=mosekopt('minimize echo(0)',prob);
sol=res.sol.bas.xx;      
solsta=strcat('MSK_SOL_STA_', res.sol.bas.solsta);
if (strcmp(solsta,'MSK_SOL_STA_PRIMAL_INFEASIBLE_CER')==1)
    disp('Unfeasible problem');
    out=1;
    opt_value=inf;
elseif (strcmp(solsta,'MSK_SOL_STA_DUAL_INFEASIBLE_CER')==1)
    disp('Primal infinite optimal value');
    out=2;
    opt_value=-inf;
else
    out=0;
    opt_value=prob.c'*sol;
end
time2=toc;
counter=counter-1;








