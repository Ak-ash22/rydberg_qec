include("functions.jl")

#Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/"
# # script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
script_dir = "/scratch/roq68sum/5atoms_code"
data_folder = joinpath(script_dir, "5_atom_correction")

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
    
    max_length = max(length(tspan2), length(tspan3))  # Choose longest possible time span
    fidelity_data = zeros(max_length)  # Initialize with zeros
    
    println("Starting the simulation...")

    for i in 1:1
        @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f,maxiters=1e9,seed=(N_trajectories*1000 + i))

        ancilla1_population = real(expect(n_1, ψt))
        ancilla2_population = real(expect(n_2, ψt))

        population_data[:a] .+= real(expect(n_a, ψt))
        population_data[:b] .+= real(expect(n_b, ψt))
        population_data[:c] .+= real(expect(n_c, ψt))
        population_data[:a1] .+= ancilla1_population
        population_data[:a2] .+= ancilla2_population

        println("Error Detection Commencing...")

        rand_float = rand()
        if rand_float ≤ ancilla1_population[end]
            println("Ancilla 1 error detected. Correcting Atom A")
            ancilla1_population[end] = 1.0

        elseif rand_float ≤ ancilla2_population[end]
            println("Ancilla 2 error detected. Correcting Atom C")
            ancilla2_population[end] = 1.0

        elseif rand_float ≤ ancilla1_population[end] && rand_float ≤ ancilla2_population[end]
            println("Both Ancilla errors detected. Correcting Atom B")
            ancilla1_population[end] = 1.0
            ancilla2_population[end] = 1.0
        end

        println("Error Correction Commencing...")

        #Correcting Atom A
        if ancilla1_population[end] == 1.0

            f1 = f_correct_factory(1)
            @time tout, ψt = timeevolution.mcwf_dynamic(tspan2,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*1000 + i))

            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        
        #Correcting Atom C
        elseif ancilla2_population[end] == 1.0
            
            f1 = f_correct_factory(3)
            @time tout, ψt = timeevolution.mcwf_dynamic(tspan3,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*1000 + i))

            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        
        #Correcting Atom B
        elseif ancilla1_population[end] == 1.0 && ancilla2_population[end] == 1.0

            f1 = f_correct_factory(2)
            @time tout, ψt = timeevolution.mcwf_dynamic(tspan2,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*1000 + i))
           
            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        end
    end

    print("Trajectory $N_trajectories Complete.\n")

    end_time = time() - start_time

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" population_data fidelity_data end_time
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
