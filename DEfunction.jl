function Kuramoto(du,u,p,t)
    λ,ω=p
    for i = 1:N
        d = 0.
        for j = 1:N
            if SC[i,j]  > 0 
                d += SC[i,j]*sin(u[j] - u[i])
            end
        end
        du[i] = ω[i] + λ*d
    end
end

function fit_r(modelFC,realFC)
    N = size(modelFC,1)
    modelFCLT = zeros(Int((N^2-N)/2))
    realFCLT = zeros(Int((N^2-N)/2))
    c = 1
    for i = 2:N
            for j = 1:i-1
                    modelFCLT[c] = modelFC[i,j]
                    realFCLT[c] = realFC[i,j]
                    c+=1
            end
    end
    return cor(modelFCLT,realFCLT)
end

function sum_eig_vecs(SC,indxs;w=0)
    N=size(SC,1)
    vals,vecs = eigen(SC)
    lcom = zeros(N,N)
    if w == 0
        w=ones(N)
    end
    count = 1
    for i in indxs
        lcom .= lcom .+  w[count]*(vecs[:,end-i] * vecs[:,end-i]')
        count += 1
    end

    return lcom

end


function sum_eig_vecs2(SC,indxs)
    N=size(SC,1)
    vals,vecs = eigen(SC)
    lcom = zeros(N,N)
 
    for i in indxs
        lcom .= lcom .+  vals[end-i]*(vecs[:,end-i] * vecs[:,end-i]')
        
    end

    return lcom

end




