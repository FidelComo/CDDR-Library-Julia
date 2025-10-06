

%For the model Z_t(m)=\sum_{j=1}^{p_t(m)} phinp{1,m}(t,j) Z_{t-j}(m)+eta_t(m)
%we first compute alpha, beta, and theta such that, for a given t,
%Z_{t+j}(m)=\sum_{l=0}^{p_{t,j}^{max}(m)} alpha{1,m}(j+pmaxS(m)+1,l+1)Z_{t-l}(m)
%+\sum_{l=1}^{j} beta{m,j}(l) eta_{t+l}(m)
%We then compute (output of the function) the coefficients alpha, beta, and
%theta such that
%\xi_{t+j}(m)=\theta(j,m)+\sum_{l=0}^{p_{t,j}^{max}(m)} alpha{1,m}(j,l+1)\xi_{t-l}(m)
%+\sum_{l=1}^{j} beta{m,j}(l) eta_{t+l}(m)
%To be called with t=1 if observations are available until time 1.

function [alpha,beta,theta]=compute_decomposition_periodic(mu,sigma,phinp,M,TabS,TabP,t,T,pmaxS)

alpha=cell(1,M);
beta=cell(M,T-t);
theta=zeros(T-t,M);

for i=1:M
    alpha{1,i}=zeros(T-t+pmaxS(i)+1,pmaxS(i)+1);
    for p=1:T-t
        beta{i,p}=zeros(1,p);
    end
end

for i=1:M
    beta{i,1}(1)=1;
    for p=0:pmaxS(i) 
        for k=0:pmaxS(i)
	        if (k==p)
	            alpha{1,i}(pmaxS(i)-p+1,k+1)=1;
            else
 	            alpha{1,i}(pmaxS(i)-p+1,k+1)=0;
            end
        end
    end
end

for i=1:M
    for p=1:T-t
        periode=mod(t+p,TabS(i));
        if (periode==0)
            periode=TabS(i);
        end
	    for k=0:pmaxS(i)
	        for j=1:TabP{1,i}(periode)
                alpha{1,i}(p+pmaxS(i)+1,k+1)=alpha{1,i}(p+pmaxS(i)+1,k+1)+(phinp{1,i}{1,periode}(j))*(alpha{1,i}(p+pmaxS(i)+1-j,k+1));
            end
        end
        for k=1:p-1
            for j=1:min(TabP{1,i}(periode),p-k)
			    beta{i,p}(k)=beta{i,p}(k)+(phinp{1,i}{1,periode}(j))*(beta{i,p-j}(k));
            end
        end
        beta{i,p}(p)=1;
    end
end

for i=1:M
    for p=1:T-t
        periode=mod(t+p,TabS(i));
        if (periode==0)
            periode=TabS(i);
        end
        theta(p,i)=mu(periode,i);
        for k=0:pmaxS(i)
            if (t-k>=1)
                periode2=mod(t-k,TabS(i));
                if (periode2==0)
                    periode2=TabS(i);
                end
            else
                periode2=t-k+TabS(i);
            end
            alpha{1,i}(p,k+1)=(alpha{1,i}(p+pmaxS(i)+1,k+1))*(sigma(periode,i))/(sigma(periode2,i));        
            theta(p,i)=theta(p,i)-(alpha{1,i}(p,k+1))*(mu(periode2,i));
        end
        for k=1:p
            beta{i,p}(k)=(beta{i,p}(k))*(sigma(periode,i));
        end
    end
end
