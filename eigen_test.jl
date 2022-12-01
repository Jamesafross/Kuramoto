using JLD,DifferentialEquations,Random,StatsBase,Plots,LinearAlgebra,JLD,DelimitedFiles

include("Global_NetworkSetup.jl")
include("DEfunction.jl")

R2_save = load("/Users/james/Kuramoto_Oscillator_FC/R2.jld","R2_save")


c=13000
SC,dist,lags,N = networksetup(c;digits=5,nSC=1,nFC=1,N=140,density=1,normalise=false)

eigenmode = sum_eig_vecs(SC,[0,1,2,3])



save_indxs=Array{Any,1}(undef,21)
best_eig_fit = zeros(N,N,21)
for j = 1:21
    global indxs = [0]
    global fit0 = fit_r(R2_save[:,:,j],sum_eig_vecs(SC,indxs).^2)
    for i = 1:138
        indxs_test = cat(indxs,i,dims=1)
        fit_test = fit_r(R2_save[:,:,j],sum_eig_vecs(SC,indxs_test).^2)
        if fit_test > fit0
            global indxs = indxs_test
            global fit0 = fit_test
        end


    end
    save_indxs[j] = indxs
    best_eig_fit[:,:,j] = sum_eig_vecs(SC,indxs).^2
end


weight_range = collect(0:0.01:2)
fit1_test = zeros(length(weight_range))
weights_save = Array{Any,1}(undef,21)

for jj = 1:21
    global sInx = size(save_indxs[jj],1)
    global weights0 = ones(sInx)
    
    global weight_test = ones(sInx)
    fit1 = fit_r(R2_save[:,:,jj],sum_eig_vecs(SC,save_indxs[jj];w=weights0).^2)
    for i = 1:sInx
        for j = 1:length(weight_range)
            weight_test[i] = weight_range[j]
            fit1_test[j] = fit_r(R2_save[:,:,jj],sum_eig_vecs(SC,save_indxs[jj];w=weight_test).^2)
        end
        global weight_test[i] = weight_range[findfirst(fit1_test .== maximum(fit1_test))]
    end

    weights_save[jj] = weight_test
    best_eig_fit[:,:,jj] = sum_eig_vecs(SC,save_indxs[jj],w=weights_save[jj]).^2

    open("/Users/james/Kuramoto_Oscillator_FC/saved_data/modelR2_min_$(Int((jj-1)*5)).txt","w") do io ; writedlm(io,best_eig_fit[:,:,jj]);end
    
end

    

println("fit = ",fit0)
println("indexs = ",indxs)





