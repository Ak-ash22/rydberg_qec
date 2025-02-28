include("functions.jl")

#Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/shor_code_data/"
# script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
script_dir = "/scratch/roq68sum/5atoms_code"
data_folder = joinpath(script_dir, "5_atom_work")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """
    println("Running the simulation with N_trajectories = $N_trajectories")

    start_time = time()
    ψ0 = initialize_system()
    println("The system has been initialized.")
    
    ψ = Vector{Vector}(undef, N_trajectories)

    full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 

    # --- Preallocate Arrays ---
    num_timesteps = length(tspan)
    population_data = Dict(key => zeros(num_timesteps) for key in (:a, :b,  :c, :a1, :a2))
    
    println("Starting the simulation...")

    for i in 1:1
        @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f,maxiters=1e9,seed=(N_trajectories*1000 + i))

        population_data[:a] .+= real(expect(n_a, ψt))
        population_data[:b] .+= real(expect(n_b, ψt))
        population_data[:c] .+= real(expect(n_c, ψt))
        population_data[:a1] .+= real(expect(n_1, ψt))
        population_data[:a2] .+= real(expect(n_2, ψt))

    end
    print("Trajectory $N_trajectories.\n")

    end_time = time() - start_time

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" population_data end_time
    println("Data saved.")
end



# --- Parse command-line arguments ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    N_trajectories = parse(Int, ARGS[1])
    main(N_trajectories)
end
