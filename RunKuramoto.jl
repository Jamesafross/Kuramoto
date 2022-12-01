using JLD,DifferentialEquations,Random,StatsBase,Plots,LinearAlgebra,JLD
savedir="/Users/james/Kuramoto_Oscillator_FC/saved_data"
include("/Users/james/Kuramoto_Oscillator_FC/Global_NetworkSetup.jl")
include("DEfunction.jl")

init = load("$savedir/init.jld","init")

c=13000
const SC,dist,lags,N = networksetup(c;digits=5,nSC=2,nFC=1,N=140,density=1,normalise=false)

FC, missingnodes= get_mean_all_functional_data(;ROI=140,type="control")
FC = FC .- (diagm(ones(N)))

ntrial = 200
R_trials = zeros(N,N,ntrial)
ω = zeros(N)




s = sum(SC,dims=1)
function g(s,i)
    si = s[i]
    sb = maximum(s[s .>0 ])
    sa = minimum(s[s .> 0 ])
    return ((si - sa)/(sb - sa)) ^2
end

for i = 1:N
    ω[i] = (0.1) - (0.1 - 0.01)*( g(s,i))
end

ω[39] = ω[39] +0.00

for i = 1:ntrial
    tspan = (0.,50*100.)
    λ = 0.00112


    p = λ,ω
    prob = ODEProblem(Kuramoto, 0.1*init[:,i], tspan, p)

    println("solving...")

    global sol = solve(prob,BS3(),maxiters = 1e20,saveat=collect(100*10:0.01:tspan[2]))

    sol = cos.(sol[:,:])


    R = zeros(N,N)
    for i = 1:N
        for j = i+1:N
            R[i,j] = cor(sol[i,:],sol[j,:])
            R[j,i] = R[i,j]
        end
    end

    R_trials[:,:,i] = R

end



print(fit_r(mean(R_trials,dims=3)[:,:].^2,FC.^2))

p1 = heatmap(sol,color=:jet)
p2 = heatmap(mean(R_trials,dims=3)[:,:].^2,color=:jet)

plot(p1,p2)