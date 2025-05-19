include("functions.jl")

#Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/"
script_dir = "/Users/akashmalemath/Documents/master_work/rydberg_qec/codework/phaseflip_code"
# script_dir = "/scratch/roq68sum/5atoms_code/phaseflip_code"
data_folder = joinpath(script_dir, "results_data")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """
    println("Running the simulation with N_trajectories = $N_trajectories")

    i = N_trajectories

    start_time = time()
    ψ0 = initialize_system()
    println("The system has been initialized.")

    full_basis = CompositeBasis([NLevelBasis(3) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 
    
    ψ_ideal = α .* reduce(kron,[b,b,b,a,a]) + β .* reduce(kron,[a,a,a,a,a])
    ψ_target = ψ_ideal / norm(ψ_ideal);

    # # --- Preallocate Arrays ---
    population_data = Dict(key => zeros(length(tspan1)) for key in (:a, :b, :c))
    # σx_exp = Dict(key => [] for key in (:a,:b,:c))


    # corrected_σx_exp = Dict(key => [] for key in (:a,:b,:c))
    # has_error = false
    # detected_error = false
    # has_correction_error = false
    
    println("Starting the simulation...")
 
    ################################################################################################################ Protocol Starts

    println("Encoding Commencing...")
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan1,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)  # pass through to `solve`
    println("Encoding complete.")

    evolved_fidelity = abs(ψt[end].data' * ψ_target)^2
    println("Fidelity of the target state with the final state is $(evolved_fidelity)\n")

    # ######################################################## Dynamical Phase Correction
    # Δ_dyn,_,_ = compute_dynamical_phase(tout)
    # ψt_end_corrected = apply_dynamical_phase(Δ_dyn, ψt[end])
    # ψt_end_dyn = ψt_end_corrected.data

    # encoding_fidelity = abs(ψt_end_dyn' * ψ_target)^2
    # println("Adiabatic sweep successfull with encoding fidelity $(encoding_fidelity)\n")

    #Track Jumps Info
    # has_error = length(jumps) > 0
    population_data[:a][1:length(tout)] .+= real(expect(n1_atom1, ψt))
    population_data[:b][1:length(tout)] .+= real(expect(n1_atom2, ψt))
    population_data[:c][1:length(tout)] .+= real(expect(n1_atom3, ψt))

    println("Phaseflip code encoding done successfully.\n")

    # println("Applying hadamards...")
    # @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψt_end_corrected,f2,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)

    # append!(σx_exp[:a], real(expect(σx_small_a, ψt)))
    # append!(σx_exp[:b], real(expect(σx_small_b, ψt)))
    # append!(σx_exp[:c], real(expect(σx_small_c, ψt)))
    # append!(population_data[:a1], zeros(length(tspan2)))
    # append!(population_data[:a2], zeros(length(tspan2)))

    # println("Hadamards applied.\n")

    # println("Driving Ancillas ...")
    # ψ3_0 = kron(ψt_end_dyn, g,g)
    # full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    # ψ3_ket = Ket(full_basis, ComplexF32.(ψ3_0)) 

    # @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψ3_ket,f3,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)

    # append!(σx_exp[:a], real(expect(σx_a, ψt)))
    # append!(σx_exp[:b], real(expect(σx_b, ψt)))
    # append!(σx_exp[:c], real(expect(σx_c, ψt)))
    # ancilla1_population = real(expect(n_1, ψt))
    # ancilla2_population = real(expect(n_2, ψt))

    # append!(population_data[:a1], ancilla1_population)
    # append!(population_data[:a2], ancilla2_population)
    # println("Ancillas drive complete.")

    # # has_error = length(jumps) > 0

    # ####################################################################################################### Error Detection
    # println("Error Detection Commencing...")

    # rand_float = round(rand();digits=1)

    # if rand_float < round(ancilla1_population[end];digits=1) && rand_float < round(ancilla2_population[end];digits=1)
    #     println("Both Ancilla errors detected. Correcting Atom B")
    #     ancilla1 = 1.0
    #     ancilla2 = 1.0
    #     detected_error = true
        
    # elseif rand_float < round(ancilla2_population[end];digits=1)
    #     println("Ancilla 2 error detected. Correcting Atom C")
    #     ancilla1 = 0.0
    #     ancilla2 = 1.0
    #     detected_error = true

    # elseif rand_float < round(ancilla1_population[end];digits=1)
    #     println("Ancilla 1 error detected. Correcting Atom A")
    #     ancilla1 = 1.0
    #     ancilla2 = 0.0
    #     detected_error = true

    # else 
    #     println("No errors detected.")
    #     ancilla1 = 0.0
    #     ancilla2 = 0.0
    # end
    # ################################################################################################### Error Correction
    # println("Error Correction Commencing...")

    # ψ1 = ψt[end]

    # #Correcting Atom B
    # if ancilla1 == 1.0 && ancilla2 == 1.0

    #     fb = f_correct_factory(2)
    #     @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fb,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    #     # has_correction_error = length(jumps) > 0

    #     # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
    #     append!(corrected_σx_exp[:a], real(expect(σx_a, ψt)))
    #     append!(corrected_σx_exp[:b], real(expect(σx_b, ψt)))
    #     append!(corrected_σx_exp[:c], real(expect(σx_c, ψt)))

    #     # corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
    #     # corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
    #     # corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
    #     # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
    #     # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))

    # #Correcting Atom A
    # elseif ancilla1 == 1.0

    #     fa = f_correct_factory(1)
    #     @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fa,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    #     # has_correction_error = length(jumps) > 0

    #     append!(corrected_σx_exp[:a], real(expect(σx_a, ψt)))
    #     append!(corrected_σx_exp[:b], real(expect(σx_b, ψt)))
    #     append!(corrected_σx_exp[:c], real(expect(σx_c, ψt)))

    #     # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
    #     # corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
    #     # corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
    #     # corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
    #     # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
    #     # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))

    # #Correcting Atom C
    # elseif ancilla2 == 1.0
        
    #     fc = f_correct_factory(3)
    #     @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fc,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    #     # has_correction_error = length(jumps) > 0

    #     append!(corrected_σx_exp[:a], real(expect(σx_a, ψt)))
    #     append!(corrected_σx_exp[:b], real(expect(σx_b, ψt)))
    #     append!(corrected_σx_exp[:c], real(expect(σx_c, ψt)))
    #     # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
    #     # corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
    #     # corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
    #     # corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
    #     # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
    #     # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))


    # #No Errors Detected
    # else
    #     f0 = f_correct_factory(0)
    #     @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,f0,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    #     # has_correction_error = length(jumps) > 0

    #     # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
    #     append!(corrected_σx_exp[:a], real(expect(σx_a, ψt)))
    #     append!(corrected_σx_exp[:b], real(expect(σx_b, ψt)))
    #     append!(corrected_σx_exp[:c], real(expect(σx_c, ψt)))
    #     # corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
    #     # corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
    #     # corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
    #     # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
    #     # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
    # end

    # println("Applying hadamards back...")
    # @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan6,ψt[end],f_end,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)

    # append!(corrected_σx_exp[:a], real(expect(σx_a, ψt)))
    # append!(corrected_σx_exp[:b], real(expect(σx_b, ψt)))
    # append!(corrected_σx_exp[:c], real(expect(σx_c, ψt)))
    # println("Hadamards applied.")


    print("Trajectory $N_trajectories Complete.\n")

    end_time = time() - start_time

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_dephase=$(γ_dephase)_Ntraj=$(N_trajectories).jld2" population_data end_time
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
