using LinearAlgebra
using DifferentialEquations
using JLD2
include("functions.jl")
include("systemparams.jl")

using .functions
using .systemparams


const Ω1, Ω2, γ_Decay, γ_dephase, V1_nn, V2_nn, δ, Δ1_0, Δ2_0, T_optimal = unpack_params()


# RAP Fidelity with different decay and sweep_rate
decay = 10 .^ range(-5,0,length=6);
sweep_rate = collect(range(0.01,1,length=20));
l = length(decay)
m = length(sweep_rate)

ψ_ideal = α*kron(r,r,r) + β*kron(g,g,g)
fidelity = Array{Float64}(undef,l,m)  # Ensure fidelity is pre-allocated

# Flatten the nested loops into a single loop that can be parallelized
@time begin
    Threads.@threads for idx in 1:l*m
        i = (idx-1) ÷ m + 1  # Calculate the row index
        j = (idx-1) % m + 1  # Calculate the column index

        a = decay[i]
        b = sweep_rate[j]

        p = Parameters(Ω1,Ω2,a,γ_dephase,V1_nn,V2_nn,b)
        tspan = (0.0, T_optimal)

        solution = solve_master_eqn(p, tspan)
        f = tr(ψ_ideal' * solution * ψ_ideal)
        fidelity[i, j] = real(f)
    end
end

@save "trial.jld2" decay sweep_rate fidelity
