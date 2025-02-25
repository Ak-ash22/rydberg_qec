include("functions.jl")

# Saving the output
script_dir = "/scratch/roq68sum/shor_code_data"

data_folder = joinpath(script_dir, "driving_abc9")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """
    println("Running the simulation for $N_trajectories")

    start_time = time()

    ψ0 = initialize_system()
    println("The system has been initialized.")

    full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 
  
    # --- Preallocate Arrays ---
    num_timesteps = length(tspan)
    population_data = Dict(
        :a => zeros(num_timesteps),
        :c => zeros(num_timesteps),
        :ac => zeros(num_timesteps),
        :p1 => zeros(num_timesteps),
        :p2 => zeros(num_timesteps),
        :p3 => zeros(num_timesteps),
        :p4 => zeros(num_timesteps),
        :p5 => zeros(num_timesteps),
        :p6 => zeros(num_timesteps)
    )
    
    println("Starting the simulation...")

    for i in 1:10
        @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f;maxiters=1e9,seed=(N_trajectories*100 + i))

        population_data[:a] .+= real(expect(n_a, ψt))
        population_data[:c] .+= real(expect(n_c, ψt))
        population_data[:p1] .+= real(expect(n_1, ψt))
        population_data[:p2] .+= real(expect(n_2, ψt))
        population_data[:p3] .+= real(expect(n_3, ψt))
        population_data[:p4] .+= real(expect(n_4, ψt))
        population_data[:p5] .+= real(expect(n_5, ψt))
        population_data[:p6] .+= real(expect(n_6, ψt))
    end

    println("Simulation complete. Calculating the average wavefunction...")

    println("Trajectory $N_trajectories.\n")
    
    # Normalize population data
    for key in keys(population_data)
        population_data[key] .*= 1/10
    end

    end_time = time() - start_time

    output_file = joinpath(data_folder, "N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2")
    @save output_file population_data end_time

    println("Data saved successfully!")
end


# --- Parse command-line arguments ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    N_trajectories = parse(Int, ARGS[1])
    main(N_trajectories)
end
