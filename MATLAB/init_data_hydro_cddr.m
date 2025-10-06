
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%depth: depth of the decision rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T: number of stages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M: number of realizations for each stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%option: if equal to one then s.d. of noise for stage t=1,..,12, and subsystem m
%        is estimated s.d. of noises using real data otherwise
%        sigmasN is used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The outputs are the inputs of the functions computing the constant depth
%decision rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [subi_a,subj_a,valij_a,betas,cost,ds,ps,qs,probabilities]=init_data_hydro_cddr(depth,T,M,path,AddP,Inflow_Noises)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nb(i) is the number of thermal plants in subsystem i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nb=[2;2;2;2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NS is the number of subsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NS=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gamma(t,i) is the proportion of inflows coming to reservoir i at t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma=zeros(T,NS);
for t=1:T
    gamma(t,:)=[0.8567,0.9013,0.9756,1];
end

%Demand=1.25*ones(T,1)*[31055,8297,7103,3367];
Demand=ones(T,1)*[31055,8297,7103,3367];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inflow_Noises(t,m,j) is realization j of noise for t+1 and subsytem m, t=1,..,T-1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vmaxH: vmaxH(t,i) is the maximal hydraulic available power at t for subsystem i.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aux=[40943.0,10292.5,9406.0,6915.7];
vmaxH=ones(T,1)*aux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vminH: vminH(t,i) is the minimal hydro production at t i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vminH=zeros(T,NS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%evap:evap(t,i) is the evaporation for time step t and subsytem i.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evap=zeros(T,NS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vmaxT: vmaxT{t,i} is the list of maximal production powers for the thermal
%and nuclear units of subsytem i and for time step t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Thermal_Costs{t,m} is the list of unit thermal costs of thermal plants of
%subsytem m for stage t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Thermal_Costs=cell(T,NS);

for t=1:T
    Thermal_Costs{t,1}=[17.48;12.61;937;6.27;100.4;77.46;432.7;150;10.5;42.6;74.4;108;180;395.87;180;254.14;288.51;288.51;589.31;97.15;124.77;110.48;1047.38;197.85;4170.44];
    vmaxT{t,1}=[417.37;1212.6;26.53;342.28;26.32;0;108.72;79.39;384;93.12;186.24;91.28;0;0;0;0;115.64;343.93;152.12;0;409.02;185.32;8;152.03;100000000];
    vminT{t,1}=zeros(25,1);
    Thermal_Costs{t,2}=[546.4;219;110.48;191.08;193.25;200.17;160.03;155;116.1;596;116;116;249;92.49;4170.44];
    vmaxT{t,2}=[49.26;440.32;0;54.37;7.87;48.50;105.14;219.19;164.66;20.06;76.65;180.53;11.44;204.7;100000000];    
    vminT{t,2}=zeros(15,1);
    Thermal_Costs{t,3}=[161.51;71.29;80.02;87.12;82.72;65;4170.44];
    vmaxT{t,3}=[307.18;51.17;327.26;102.29;211.48;83.94;100000000];
    vminT{t,3}=zeros(7,1);
    Thermal_Costs{t,4}=[4170.44];
    vmaxT{t,4}=[100000000];
    vminT{t,4}=0;
end

%Update Thermal Costs to have one equivalent TP per subsystem
for t=1:T
    Thermal_Costs{t,1}=[max(Thermal_Costs{t,1}(1:Nb(1)-1));4170.44];
    vmaxT{t,1}=[0.6*sum(Demand(:,1))/T;inf];
    vminT{t,1}=[0;0];
    Thermal_Costs{t,2}=[max(Thermal_Costs{t,2}(1:Nb(2)-1));4170.44];
    vmaxT{t,2}=[0.6*sum(Demand(:,2))/T;inf];
    vminT{t,2}=[0;0];
    Thermal_Costs{t,3}=[max(Thermal_Costs{t,3}(1:Nb(3)-1));4170.44];
    vmaxT{t,3}=[0.6*sum(Demand(:,3))/T;inf];
    vminT{t,3}=[0;0];
    Thermal_Costs{t,4}=[max(Thermal_Costs{t,1}(1:Nb(1)-1));4170.44];
    vmaxT{t,4}=[0.6*sum(Demand(:,4))/T;inf];
    vminT{t,4}=[0;0];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mus(t,i) is empirical mean of inflow for subsystem i month t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sigmaInflows(t,i) is empirical s.d. of inflow for subsytem i month t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%phi: Normalized coefficients phi of PAR models phi{1,i}{1,t} are
%coefficients for subsystem i and month t.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mus=zeros(12,NS);
sigmaInflows=zeros(12,NS);
pathaux=strcat(path,'Historique_Apports_Sud_Est.txt');
Historique_SE=dlmread(pathaux);
for i=1:12
    Mus(i,1)=mean(Historique_SE(:,i));
    sigmaInflows(i,1)=sqrt(var(Historique_SE(:,i)));
end
pathaux=strcat(path,'Historique_Apports_Sud.txt');
Historique_S=dlmread(pathaux);
for i=1:12
    Mus(i,2)=mean(Historique_S(:,i));
    sigmaInflows(i,2)=sqrt(var(Historique_S(:,i)));
end
pathaux=strcat(path,'Historique_Apports_Nord_Est.txt');
Historique_NE=dlmread(pathaux);
for i=1:12
    Mus(i,3)=mean(Historique_NE(:,i));
    sigmaInflows(i,3)=sqrt(var(Historique_NE(:,i)));
end
pathaux=strcat(path,'Historique_Apports_Nord.txt');
Historique_Nord=dlmread(pathaux);
for i=1:12
    Mus(i,4)=mean(Historique_Nord(:,i));
    sigmaInflows(i,4)=sqrt(var(Historique_Nord(:,i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recent_History_Inflows(-t+2,:) inflows for time t (t=1,0,-1,-2,...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Recent_History_Inflows(-t+2,:) inflows for time t (t=1,0,-1,-2,...)
Recent_History_Inflows(1,1)=Historique_SE(size(Historique_SE,1),12);
for t=1:12
    Recent_History_Inflows(14-t,1)=Historique_SE(size(Historique_SE,1),t);
end

Recent_History_Inflows(1,2)=Historique_S(size(Historique_S,1),12);
for t=1:12
    Recent_History_Inflows(14-t,2)=Historique_S(size(Historique_S,1),t);
end

Recent_History_Inflows(1,3)=Historique_NE(size(Historique_NE,1),12);
for t=1:12
    Recent_History_Inflows(14-t,3)=Historique_NE(size(Historique_NE,1),t);
end

Recent_History_Inflows(1,4)=Historique_Nord(size(Historique_Nord,1),12);
for t=1:12
    Recent_History_Inflows(14-t,4)=Historique_Nord(size(Historique_Nord,1),t);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TabP: TabP is a cell array of size NS. TabP{1,i} is the list of PAR model
%i periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TabS:TabS(i) is the period of PAR model i.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


phinp=cell(1,NS);
TabS=[12,12,12,12];
TabP=cell(1,NS);
TabP{1,1}=[1,1,1,2,3,1,3,1,1,3,1,4];
TabP{1,2}=[1,1,1,1,1,1,4,1,1,1,1,1];
TabP{1,3}=[5,2,1,1,1,1,2,1,3,3,2,5];
TabP{1,4}=[1,4,1,1,2,1,3,2,5,3,5,1];

for m=1:NS
    phinp{1,m}=cell(1,12);
end

phinp{1,1}{1,1}=[0.593];
phinp{1,1}{1,2}=[0.569];
phinp{1,1}{1,3}=[0.671];
phinp{1,1}{1,4}=[0.601,0.274];
phinp{1,1}{1,5}=[0.617,-0.0484,0.332];
phinp{1,1}{1,6}=[0.824];
phinp{1,1}{1,7}=[0.699,0.0108,0.284];
phinp{1,1}{1,8}=[0.838];
phinp{1,1}{1,9}=[0.843];
phinp{1,1}{1,10}=[0.321,0.130,0.336];
phinp{1,1}{1,11}=[0.730];
phinp{1,1}{1,12}=[0.708,-0.203,-0.0600,0.387];

phinp{1,2}{1,1}=[0.399];
phinp{1,2}{1,2}=[0.568];
phinp{1,2}{1,3}=[0.665];
phinp{1,2}{1,4}=[0.508];
phinp{1,2}{1,5}=[0.468];
phinp{1,2}{1,6}=[0.642];
phinp{1,2}{1,7}=[0.427,0.350,-0.239,0.361];
phinp{1,2}{1,8}=[0.457];
phinp{1,2}{1,9}=[0.569];
phinp{1,2}{1,10}=[0.474];
phinp{1,2}{1,11}=[0.54];
phinp{1,2}{1,12}=[0.587];

phinp{1,3}{1,1}=[0.709,-0.183,0.134,0.281,-0.233];
phinp{1,3}{1,2}=[0.775,-0.344];
phinp{1,3}{1,3}=[0.777];
phinp{1,3}{1,4}=[0.707];
phinp{1,3}{1,5}=[0.827];
phinp{1,3}{1,6}=[0.948];
phinp{1,3}{1,7}=[1.2,-0.248];
phinp{1,3}{1,8}=[0.979];
phinp{1,3}{1,9}=[1.12,0.0603,-0.247];
phinp{1,3}{1,10}=[0.760,0.679,-0.604];
phinp{1,3}{1,11}=[0.971,-0.336];
phinp{1,3}{1,12}=[0.708,-0.0534,-0.0342,-0.614,0.606];

phinp{1,4}{1,1}=[0.739];
phinp{1,4}{1,2}=[0.859,-0.449,-0.0423,0.286];
phinp{1,4}{1,3}=[0.778];
phinp{1,4}{1,4}=[0.774];
phinp{1,4}{1,5}=[0.995,-0.227];
phinp{1,4}{1,6}=[0.895];
phinp{1,4}{1,7}=[1.12,-0.528,0.339];
phinp{1,4}{1,8}=[1.18,-0.228];
phinp{1,4}{1,9}=[1.27,-0.754,0.132,-0.0252,0.321];
phinp{1,4}{1,10}=[0.476,0.672,-0.335];
phinp{1,4}{1,11}=[0.639,-0.295,0.822,-0.246,-0.251];
phinp{1,4}{1,12}=[0.708];

for i=1:4 
    for j=1:12
        TabP{1,i}(j)=TabP{1,i}(j)+AddP;
        phinp{1,i}{1,j}=[phinp{1,i}{1,j},0.01*ones(1,AddP)];
    end
end

phi=cell(1,NS);

for m=1:NS
    phi{1,m}=cell(1,TabS(m));
    for s=1:TabS(m)
            
        for j=1:TabP{1,m}(s)
            if (s-j>=1)
               phi{1,m}{1,s}(j)=((phinp{1,m}{1,s}(j))*sigmaInflows(s,m))/sigmaInflows(s-j,m);
            else
               phi{1,m}{1,s}(j)=((phinp{1,m}{1,s}(j))*sigmaInflows(s,m))/sigmaInflows(s-j+TabS(m),m); 
            end
        end
    end
end


%Computation of sigmasNoises in case option=1
H_SE=[];
H_S=[];
H_NE=[];
H_Nord=[];
for i=1:size(Historique_SE,1)
    H_SE=[H_SE;Historique_SE(i,:)'];
end
for i=1:size(Historique_S,1)
    H_S=[H_S;Historique_S(i,:)'];
end
for i=1:size(Historique_NE,1)
    H_NE=[H_NE;Historique_NE(i,:)'];
end
for i=1:size(Historique_Nord,1)
    H_Nord=[H_Nord;Historique_Nord(i,:)'];
end

Echantillons=[H_SE,H_S,H_NE,H_Nord];

% sigmasNoises=zeros(12,4);
% for t=1:12
%     %Computation of sigmasNoises(t,m)
%     %Computation of sample of noise for month t subsystem m
%     for m=1:NS
%         Sample=[];
%         for k=1:((size(Echantillons,1)/12)-1)
%             noisevalue=(Echantillons(t+12*k,m)-Mus(t,m))/sigmaInflows(t,m);
%             for j=1:TabP{1,m}(t)
%                 if (t-j>=1)
%                     aux=(Echantillons(t+12*k-j,m)-Mus(t-j,m))/sigmaInflows(t-j,m);
%                 else
%                     aux=(Echantillons(t+12*k-j,m)-Mus(t-j+TabS(m),m))/sigmaInflows(t-j+TabS(m),m);
%                 end
%                 noisevalue=noisevalue-phinp{1,m}{1,t}(j)*aux;
%             end
%             Sample=[Sample,noisevalue];
%         end
%         if (option==1)
%             sigmaNoises(t,m)=sqrt(var(Sample));
%         else
%             sigmaNoises(t,m)=sigmasN;
%         end
%     end
% end

for t=1:T-1
        for j=1:M
            probabilities(t,j)=1/M;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambda: xmin(t,m) is lambda(t,m)xmax(t,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lambda=zeros(T,NS);
auxiliaire=[33;35;38;37;38;36;32;27;20;14;10;10]/100;
i=1;
lambda((i-1)*12+1:i*12,1)=auxiliaire;
auxiliaire=[19;20;19;16;13;13;13;13;13;13;13;13]/100;
lambda((i-1)*12+1:i*12,2)=auxiliaire;
auxiliaire=[33;36;39;40;37;35;30;24;18;13;10;10]/100;
lambda((i-1)*12+1:i*12,3)=auxiliaire;

lambda=0.2*ones(T,NS);

% aux=[197263.3,18374.0,51806.1,12415.2];
% xmax=ones(T,1)*aux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%probabilities(t,j): probability of realization Inflow_Noises(t,:,j) for
%                    t+1,t=1,...,T-1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=1;
for m=1:NS
    reste=mod(j+1,TabS(m));
    if (reste==0)
        pmaxS(m)=TabP{1,m}(TabS(m))-1;
    else
        pmaxS(m)=TabP{1,m}(reste)-1;
    end
    for s=2:TabS(m)
        reste=mod(j+s,TabS(m));
        if (reste==0)
            val=TabP{1,m}(TabS(m))-s;
        else
            val=TabP{1,m}(reste)-s;
        end
        if (val>pmaxS(m))
            pmaxS(m)=val;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xmax: xmax(t,i) is the maximal level of subsystem i reservoir at t.
%x0:x0(i) is the initial level subsystem i reservoir.
%Demand: a (T,NS) matrix, demande(t,i) is the deterministic demand for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute reduced demand
x0=zeros(NS,1);
for m=1:NS
    x0(m)=0.3*sum(Demand(:,m));
end
xmax=2*ones(T,1)*x0';

%x0=[197263.3;18374.0;51806.1;12415.2];

ds=[1;M*ones(T-1,1)];
ps=9*NS*ones(T,1);
qs=4*NS*ones(T,1);

subi_a=cell(1,T);
subj_a=cell(1,T);
valij_a=cell(1,T);
for t=1:T
    subi_a{1,t}=cell(1,t);
    subj_a{1,t}=cell(1,t);
    valij_a{1,t}=cell(1,t);
    
    subi_a{1,t}{1,t}=[];
    subj_a{1,t}{1,t}=[];
    valij_a{1,t}{1,t}=[];
    row=0;
    col=0;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},ones(1,NS)];
    
    col=col+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},ones(1,NS)];
    
    row=row+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    col=col+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    col=col+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    row=row+NS;
    col=0;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},ones(1,NS)];
    
    row=row+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    row=row+NS;
    col=col+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},ones(1,NS)];
    
    row=row+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    row=row+NS;
    col=col+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},ones(1,NS)];
    
    row=row+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    row=row+NS;
    col=col+NS;
    subi_a{1,t}{1,t}=[subi_a{1,t}{1,t},[row+1:1:row+NS]];
    subj_a{1,t}{1,t}=[subj_a{1,t}{1,t},[col+1:1:col+NS]];
    valij_a{1,t}{1,t}=[valij_a{1,t}{1,t},-ones(1,NS)];
    
    if (t>=2)
        subi_a{1,t}{1,t-1}=[1:1:NS];
        subj_a{1,t}{1,t-1}=[1:1:NS];
        valij_a{1,t}{1,t-1}=-ones(1,NS);
    end
end

ds=[ones(depth-1,1);ds];
productsu=zeros(T,depth);
for s=1:T
    productsu(s,depth)=ds(s+depth-1);
    for i=1:(depth-1)
        productsu(s,depth-i)=productsu(s,depth-i+1)*ds(s+depth-i-1);
    end
end
ds=[1;M*ones(T-1,1)];

[alpha,beta,theta]=compute_decomposition_periodic(Mus,sigmaInflows,phinp,NS,TabS,TabP,1,T,pmaxS);

betas=cell(1,T);
for t=1:T
    betas{1,t}=cell(1,t);
    if (t==1)
        betas{1,1}{1,1}=[x0+gamma(1,:)'.*Recent_History_Inflows(1,:)';(1-gamma(1,:)').*Recent_History_Inflows(1,:)'-Demand(1,:)';xmax(1,:)';zeros(NS,1);vmaxH(1,:)';zeros(NS,1);vmaxT{1,1}(1);vmaxT{1,2}(1);vmaxT{1,3}(1);vmaxT{1,4}(1);zeros(2*NS,1)]';
    else
        aux=zeros(NS,1);
        for m=1:NS
            aux(m)=theta(t-1,m);
            for ell=0:pmaxS(m)
                aux(m)=aux(m)+alpha{1,m}(t-1,ell+1)*Recent_History_Inflows(ell+1,m);
            end
        end
        betas{1,t}{1,1}=[aux.*gamma(t,:)';(1-gamma(t,:)').*aux-Demand(t,:)';xmax(t,:)';zeros(NS,1);vmaxH(t,:)';zeros(NS,1);vmaxT{t,1}(1);vmaxT{t,2}(1);vmaxT{t,3}(1);vmaxT{t,4}(1);zeros(2*NS,1)]';
        for s=2:t
            betas{1,t}{1,s}=zeros(productsu(s,1),9*NS);
            nb_blocks=productsu(s,1)/M;
            for i=1:nb_blocks
                for j=1:M
                    aux=zeros(1,NS);
                    for m=1:NS
                        aux(m)=beta{m,t-1}(s-1)*Inflow_Noises(s-1,m,j);
                    end
                    betas{1,t}{1,s}((i-1)*M+j,:)=[gamma(t,:).*aux,(1-gamma(t,:)).*aux,zeros(1,7*NS)];
                end
            end
        end
    end
end

cost=cell(1,T);
for t=1:T
    cost{1,t}=[zeros(2*NS,1);Thermal_Costs{t,1}(1);Thermal_Costs{t,2}(1);Thermal_Costs{t,3}(1);Thermal_Costs{t,4}(1);Thermal_Costs{t,1}(2);Thermal_Costs{t,2}(2);Thermal_Costs{t,3}(2);Thermal_Costs{t,4}(2)];
end

subi_a=jsonencode(subi_a);
subj_a=jsonencode(subj_a);
valij_a=jsonencode(valij_a);