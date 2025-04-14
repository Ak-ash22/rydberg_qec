include("functions.jl")

#Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/"
# # script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
script_dir = "/scratch/roq68sum/5atoms_code"
data_folder = joinpath(script_dir, "5_atom_correction/decay_1e_5")

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

    full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 

    # --- Preallocate Arrays ---
    num_timesteps = length(tspan)
    population_data = Dict(key => zeros(num_timesteps) for key in (:a, :b,  :c, :a1, :a2, :abc))
    # print(num_timesteps)
    # print(length(population_data[:a]))
    
    max_length = max(length(tspan2), length(tspan3))  # Choose longest possible time span
    fidelity_data = zeros(max_length)  # Initialize with zeros

    corrected_population_data = Dict(key => zeros(max_length) for key in (:a, :b,  :c, :a1, :a2))
    has_error = false
    detected_error = false
    has_correction_error = false
    
    println("Starting the simulation...")

    for i in 1:1
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f,maxiters=1e9,seed=(N_trajectories*10000),display_jumps=true)

        #Track Jumps Info
        has_error = length(jumps) > 0

        ancilla1_population = real(expect(n_1, ψt))
        ancilla2_population = real(expect(n_2, ψt))

        population_data[:a] .+= real(expect(n_a, ψt))
        population_data[:b] .+= real(expect(n_b, ψt))
        population_data[:c] .+= real(expect(n_c, ψt))
        population_data[:abc] .+= real(expect(n_abc, ψt))
        population_data[:a1] .+= ancilla1_population
        population_data[:a2] .+= ancilla2_population

        println("Error Detection Commencing...")

        rand_float = round(rand();digits=1)

        if rand_float < round(ancilla1_population[end];digits=1) && rand_float < round(ancilla2_population[end];digits=1)
            println("Both Ancilla errors detected. Correcting Atom B")
            ancilla1 = 1.0
            ancilla2 = 1.0
            detected_error = true
            
        elseif rand_float < round(ancilla2_population[end];digits=1)
            println("Ancilla 2 error detected. Correcting Atom C")
            ancilla1 = 0.0
            ancilla2 = 1.0
            detected_error = true

        elseif rand_float < round(ancilla1_population[end];digits=1)
            println("Ancilla 1 error detected. Correcting Atom A")
            ancilla1 = 1.0
            ancilla2 = 0.0
            detected_error = true
        
        else 
            println("No errors detected.")
            ancilla1 = 0.0
            ancilla2 = 0.0
        end

        println("Error Correction Commencing...")

        ψ1 = ψt[end]

        #Correcting Atom B
        if ancilla1 == 1.0 && ancilla2 == 1.0

            f1 = f_correct_factory(2)
            @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψ1,f1,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
            has_correction_error = length(jumps) > 0
        
            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
            corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
            corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
            corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
            corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
            corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))

        #Correcting Atom A
        elseif ancilla1 == 1.0

            f1 = f_correct_factory(1)
            @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψ1,f1,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
            has_correction_error = length(jumps) > 0

            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
            corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
            corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
            corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
            corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
            corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
        
        #Correcting Atom C
        elseif ancilla2 == 1.0
            
            f1 = f_correct_factory(3)
            @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψ1,f1,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
            has_correction_error = length(jumps) > 0

            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
            corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
            corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
            corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
            corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
            corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
        
        
        #No Errors Detected
        else
            f1 = f_correct_factory(0)
            @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψ1,f1,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
            has_correction_error = length(jumps) > 0

            fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
            corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
            corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
            corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
            corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
            corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
        end
    end

    print("Trajectory $N_trajectories Complete.\n")

    end_time = time() - start_time

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" has_error detected_error has_correction_error population_data corrected_population_data fidelity_data end_time
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
