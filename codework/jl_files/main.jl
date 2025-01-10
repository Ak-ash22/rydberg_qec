using LinearAlgebra
using DifferentialEquations
using JLD2, FileIO

include("functions.jl")
include("systemparams.jl")
using .functions
using .systemparams

#Saving the output
script_dir = "/home/agfleischhauer/roq68sum/rydberg_qec/codework"
data_folder = joinpath(script_dir, "results_data/$(n_atoms)atoms")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end


const Ω1, Ω2, γ_Decay, γ_dephase, V1_nn, V2_nn, δ, Δ1_0, Δ2_0, T_optimal = unpack_params()


function case1()
    """
    Case 1: Optimal fidelity with respect to decay and sweep rate
    """
    println("Running for Optimal fidelity with respect to decay and sweep rate")


    @time begin
        decay = 10 .^ range(-5,-2,length=4);
        sweep_rate = collect(range(0.01,0.2,length=40));
        l = length(decay)
        m = length(sweep_rate)
    
        ψ_ideal = α*kron(r,r,r) + β*kron(g,g,g)
        fidelity = Array{Float64}(undef,l,m)  # Ensure fidelity is pre-allocated
    
        Threads.@threads for idx in 1:l*m
            i = (idx-1) ÷ m + 1  # Calculate the row index
            j = (idx-1) % m + 1  # Calculate the column index

            a = decay[i]
            b = sweep_rate[j]

            p = Parameters(Ω1,Ω2,a,γ_dephase,V1_nn,V2_nn,b)
            tspan = (0.0, T_optimal)

            @time solution = solve_master_eqn(p, tspan)
            f = tr(ψ_ideal' * solution * ψ_ideal)
            fidelity[i, j] = real(f)
        end
    end
 
    @save "$(data_folder)/optimal_decay.jld2" decay sweep_rate fidelity
end


function main()
    println("Choose a case to run: ")
    println("1: Optimal fidelity with respect to decay and sweep rate")
    println("Enter the case number: ")

    choice = readline()

    try
        choice = parse(Int64, choice)
        if choice == 1
            case1()
        else
            println("Bruh! Enter a valid choice")
        end
    catch e
        println("Bruh! Enter a valid choice")
    end
end

main()