#For the model Z_t(m)=\sum_{j=1}^{p_t(m)} phinp{1,m}(t,j) Z_{t-j}(m)+eta_t(m)
#we first compute alpha, beta, and theta such that, for a given t,
#Z_{t+j}(m)=\sum_{l=0}^{p_{t,j}^{max}(m)} alpha{1,m}(j+pmaxS(m)+1,l+1)Z_{t-l}(m)
#+\sum_{l=1}^{j} beta{m,j}(l) eta_{t+l}(m)
#We then compute (output of the function) the coefficients alpha, beta, and
#theta such that
#\xi_{t+j}(m)=\theta(j,m)+\sum_{l=0}^{p_{t,j}^{max}(m)} alpha{1,m}(j,l+1)\xi_{t-l}(m)
#+\sum_{l=1}^{j} beta{m,j}(l) eta_{t+l}(m)
#To be called with t=1 if observations are available until time 1.

function compute_decomposition_periodic(mu,sigma,phinp,M,TabS,TabP,t,T,pmaxS)

    alpha = Vector{Matrix{Float64}}(undef, M)
    beta = [[] for i=1:M, j=1:T-t]
    theta = zeros(T-t,M)

    for i=1:M
        alpha[i] = zeros(T-t+pmaxS[i]+1,pmaxS[i]+1)
        for p=1:T-t
            beta[i,p] = zeros(p)
        end
    end

    for i=1:M
        beta[i,1][1] = 1
        for p=0:pmaxS[i] 
            for k=0:pmaxS[i]
                if k == p
                    alpha[i][pmaxS[i]-p+1,k+1] = 1
                else
                    alpha[i][pmaxS[i]-p+1,k+1] = 0
                end
            end
        end
    end

    for i=1:M
        for p=1:T-t
            periode = (t+p) % TabS[i]
            if periode == 0
                periode = TabS[i]
            end
            for k=0:pmaxS[i]
                for j=1:TabP[i][periode]
                    alpha[i][p+pmaxS[i]+1,k+1] += phinp[i][periode][j] * alpha[i][p+pmaxS[i]+1-j,k+1]
                end
            end
            for k=1:p-1
                for j=1:min(TabP[i][periode],p-k)
                    beta[i,p][k] += phinp[i][periode][j] * beta[i,p-j][k]
                end
            end
            beta[i,p][p] = 1
        end
    end

    for i=1:M
        for p=1:T-t
            periode = (t+p) % TabS[i]
            if periode == 0
                periode = TabS[i]
            end
            theta[p,i] = mu[periode,i]
            for k=0:pmaxS[i]
                if t-k >= 1
                    periode2 = (t-k) % TabS[i]
                    if periode2 == 0
                        periode2 = TabS[i]
                    end
                else
                    periode2 = t - k + TabS[i]
                end
                alpha[i][p,k+1] = alpha[i][p+pmaxS[i]+1,k+1] * sigma[periode,i] / sigma[periode2,i]        
                theta[p,i] -= alpha[i][p,k+1] * mu[periode2,i]
            end
            for k=1:p
                beta[i,p][k] *= sigma[periode,i]
            end
        end
    end

    return alpha,beta,theta
end