include("functions.jl")

# Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/"
# script_dir = "/scratch/roq68sum/5atoms_code/phaseflip_code/"
# data_folder = joinpath(script_dir, "dephase_$(γ_dephase)/")

# if !isdir(data_folder)
#     println("Directory does not exist. Creating directory...: $data_folder")
#     mkpath(data_folder)
# end

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """

    ϕ = 0.915

    println("Running the simulation with trajectory number = $(N_trajectories) and phase = $(ϕ)")
    
    i = N_trajectories
    start_time = time()
    ψ0 = initialize_system(ϕ)
    println("The system has been initialized.")

    full_basis = CompositeBasis([NLevelBasis(3) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 


    # # --- Preallocate Arrays ---
    # population_data = Dict(key => zeros(length(tspan1)+2*length(tspan2)+length(tspan3)) for key in (:a,:b,:c,:a1,:a2))
    # final_step_population_data = Dict(key => zeros(length(tspan5)) for key in (:a, :b, :c,))
    # S1_data = zeros(2*length(tspan2)+length(tspan3)+length(tspan4)+length(tspan5))
    # S2_data = zeros(2*length(tspan2)+length(tspan3)+length(tspan4)+length(tspan5))


    # corrected_population_data = Dict(key => zeros(length(tspan4)) for key in (:a,:b,:c,:a1,:a2))
    # has_encoding_error = false
    # has_error = false
    # detected_error = false
    # has_correction_error = false
    
    println("Starting the simulation...")
 
    ################################################################################################################ Protocol Starts

    println("Encoding Commencing...")
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan1,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)  # pass through to `solve`

    #Track Jumps Info
    # has_encoding_error = length(jumps) > 0
    # population_data[:a][1:length(tout)] .+= real(expect(n1_atom1, ψt))
    # population_data[:b][1:length(tout)] .+= real(expect(n1_atom2, ψt))
    # population_data[:c][1:length(tout)] .+= real(expect(n1_atom3, ψt))

    # println("Applying hadamards...")
    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψt_end,f2,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
    #Track Jumps Info
    # has_encoding_error = length(jumps) > 0
    # population_data[:a][length(tspan1)+1:length(tspan1)+length(tout)] .+= real(expect(n1_atom1,ψt))
    # population_data[:b][length(tspan1)+1:length(tspan1)+length(tout)] .+= real(expect(n1_atom2,ψt))
    # population_data[:c][length(tspan1)+1:length(tspan1)+length(tout)] .+= real(expect(n1_atom3,ψt))
    # S1_data[1:length(tout)] .+= real(expect(S1,ψt))
    # S2_data[1:length(tout)] .+= real(expect(S2,ψt))


    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)
    plus = sqrt(1/2) .* a + sqrt(1/2) .* b
    minus = sqrt(1/2) .* a - sqrt(1/2) .* b
    ψ_target2 = α .* reduce(kron,[plus,plus,plus,a,a]) + (exp(1im * ϕ) * β) .* reduce(kron,[minus,minus,minus,a,a])
    fidelity_after_encoding = abs((dagger(Ket(full_basis,ψ_target2)) * ψt_end)^2)


    println("Phaseflip code encoding done successfully with fidelity $(fidelity_after_encoding).\n")

    ################ Storage of qubits
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan6,ψt_end,f_storage,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    fidelity_after_storage = abs((dagger(Ket(full_basis,ψ_target2)) * ψt[end]/norm(ψt_end))^2)

    println("Fidelity after the storage time is $(fidelity_after_storage)")

    ########################## Syndrome Measurement starts
    println("Applying hadamards...")
    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψt_end,f2,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
    # population_data[:a][length(tspan1)+length(tspan2)+1:length(tspan1)+length(tspan2)+length(tout)] .+= real(expect(n1_atom1,ψt))
    # population_data[:b][length(tspan1)+length(tspan2)+1:length(tspan1)+length(tspan2)+length(tout)] .+= real(expect(n1_atom2,ψt))
    # population_data[:c][length(tspan1)+length(tspan2)+1:length(tspan1)+length(tspan2)+length(tout)] .+= real(expect(n1_atom3,ψt))
    # S1_data[1+length(tspan2):length(tspan2)+length(tout)] .+= real(expect(S1,ψt))
    # S2_data[1+length(tspan2):length(tspan2)+length(tout)] .+= real(expect(S2,ψt))

    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)

    ψ_target3 = α .* reduce(kron,[b,b,b,a,a]) - (exp(1im * ϕ) * β) .* reduce(kron, [a,a,a,a,a])
    fidelity = abs((dagger(Ket(full_basis,ψ_target3)) * ψt_end)^2)

    println("Fidelity after hadamard is $(fidelity)")

    # println("Driving Ancillas ...")
    # ψ3_0 = α .* reduce(kron,[b,a,b,a,a]) - β .* reduce(kron, [a,b,a,a,a])
    # ψ3_0 = Ket(full_basis, ComplexF32.(ψ3_0))

    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψt_end,f3,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
    #Track Jumps Info
    has_error = length(jumps) > 0
    ancilla1_population = real(expect(n1_ancilla1, ψt))
    ancilla2_population = real(expect(n1_ancilla2,ψt))

    # population_data[:a][length(tspan1)+2*length(tspan2)+1:length(tspan1)+2*length(tspan2)+length(tout)] = real(expect(n1_atom1,ψt))
    # population_data[:b][length(tspan1)+2*length(tspan2)+1:length(tspan1)+2*length(tspan2)+length(tout)] = real(expect(n1_atom2,ψt))
    # population_data[:c][length(tspan1)+2*length(tspan2)+1:length(tspan1)+2*length(tspan2)+length(tout)] = real(expect(n1_atom3,ψt))
    # population_data[:a1][length(tspan1)+2*length(tspan2)+1:length(tspan1)+2*length(tspan2)+length(tout)] = ancilla1_population
    # population_data[:a2][length(tspan1)+2*length(tspan2)+1:length(tspan1)+2*length(tspan2)+length(tout)] = ancilla2_population
    # S1_data[2*length(tspan2)+1:2*length(tspan2)+length(tout)] .+= real(expect(S1,ψt))
    # S2_data[2*length(tspan2)+1:2*length(tspan2)+length(tout)] .+= real(expect(S2,ψt))

    println("Ancillas drive complete.\n")

    ####################################################################################################### Error Detection
    println("Error Detection Commencing...")

    dims = Int[length(b) for b in ψt[end].basis.bases]
    ψt_end = reshape(ψt[end].data, dims...)

    rand_float1 = round(rand();digits=1)
    rand_float2 = round(rand();digits=1)

    if rand_float1 < round(ancilla1_population[end];digits=1) && rand_float2 < round(ancilla2_population[end];digits=1)#################
        println("Both Ancilla errors detected.")
        ancilla1 = 1.0
        ancilla2 = 1.0
        detected_error = true

        ψ_abc = @view ψt_end[:,:,:,2,2]
        ψ_abc = Ket(full_basis, vec(copy(ψ_abc)))
        ψ_abc = ψ_abc / norm(ψ_abc)
        
    elseif rand_float1 < round(ancilla2_population[end];digits=1)
        println("Ancilla 2 error detected.")
        ancilla1 = 0.0
        ancilla2 = 1.0
        detected_error = true

        ψ_abc = @view ψt_end[:,:,:,1,2]
        ψ_abc = Ket(full_basis, vec(copy(ψ_abc)))
        ψ_abc = ψ_abc / norm(ψ_abc)

    elseif rand_float1 < round(ancilla1_population[end];digits=1)
        println("Ancilla 1 error detected.")
        ancilla1 = 1.0
        ancilla2 = 0.0
        detected_error = true

        ψ_abc = @view ψt_end[:,:,:,2,1]
        ψ_abc = Ket(full_basis, vec(copy(ψ_abc)))
        ψ_abc = ψ_abc / norm(ψ_abc)

    else 
        println("No errors detected.")
        ancilla1 = 0.0
        ancilla2 = 0.0

        ψ_abc = @view ψt_end[:,:,:,1,1]
        ψ_abc = Ket(full_basis, vec(copy(ψ_abc)))
        ψ_abc = ψ_abc / norm(ψ_abc)
    end
    ################################################################################################### Error Correction


    ψ1 = ψ_abc

    #Correcting Atom B
    if ancilla1 == 1.0 && ancilla2 == 1.0

        println("Error Correction commencing for Atom B")

        fb = f_correct_factory(2)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fb,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        # corrected_population_data[:a][1:length(tout)] .+= real(expect(n1_atom1, ψt))
        # corrected_population_data[:b][1:length(tout)] .+= real(expect(n1_atom2, ψt))
        # corrected_population_data[:c][1:length(tout)] .+= real(expect(n1_atom3, ψt))
        # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n1_ancilla1, ψt))
        # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n1_ancilla2, ψt))
        # S1_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S1,ψt))
        # S2_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S2,ψt))

    #Correcting Atom A
    elseif ancilla1 == 1.0

        println("Error Correction commencing for Atom A")

        fa = f_correct_factory(1)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fa,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        # corrected_population_data[:a][1:length(tout)] .+= real(expect(n1_atom1, ψt))
        # corrected_population_data[:b][1:length(tout)] .+= real(expect(n1_atom2, ψt))
        # corrected_population_data[:c][1:length(tout)] .+= real(expect(n1_atom3, ψt))
        # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n1_ancilla1, ψt))
        # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n1_ancilla2, ψt))
        # S1_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S1,ψt))
        # S2_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S2,ψt))

    #Correcting Atom C
    elseif ancilla2 == 1.0

        println("Error Correction commencing for Atom C")
        
        fc = f_correct_factory(3)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fc,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        # corrected_population_data[:a][1:length(tout)] .+= real(expect(n1_atom1, ψt))
        # corrected_population_data[:b][1:length(tout)] .+= real(expect(n1_atom2, ψt))
        # corrected_population_data[:c][1:length(tout)] .+= real(expect(n1_atom3, ψt))
        # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n1_ancilla1, ψt))
        # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n1_ancilla2, ψt))
        # S1_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S1,ψt))
        # S2_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S2,ψt))


    #No Errors Detected
    else
        f0 = f_correct_factory(0)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,f0,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        # fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        
        # corrected_population_data[:a][1:length(tout)] .+= real(expect(n1_atom1, ψt))
        # corrected_population_data[:b][1:length(tout)] .+= real(expect(n1_atom2, ψt))
        # corrected_population_data[:c][1:length(tout)] .+= real(expect(n1_atom3, ψt))
        # corrected_population_data[:a1][1:length(tout)] .+= real(expect(n1_ancilla1, ψt))
        # corrected_population_data[:a2][1:length(tout)] .+= real(expect(n1_ancilla2, ψt))
        # S1_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S1,ψt))
        # S2_data[2*length(tspan2)+length(tspan3)+1:2*length(tspan2)+length(tspan3)+length(tout)] .+= real(expect(S2,ψt))
    end

    println("Applying hadamards back...")
    ψt_end = ψt[end]/norm(ψt[end])
    ψt_end = virtual_z_full * (ψt_end)
    ψt_end /= norm(ψt_end)
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan5,ψt_end,f5,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
    # final_step_population_data[:a][1:length(tout)] .+= real(expect(n1_atom1,ψt))
    # final_step_population_data[:b][1:length(tout)] .+= real(expect(n1_atom2,ψt))
    # final_step_population_data[:c][1:length(tout)] .+= real(expect(n1_atom3,ψt))
    # S1_data[2*length(tspan2)+length(tspan3)+length(tspan4)+1:2*length(tspan2)+length(tspan3)+length(tspan4)+length(tout)] .+= real(expect(S1,ψt))
    # S2_data[2*length(tspan2)+length(tspan3)+length(tspan4)+1:2*length(tspan2)+length(tspan3)+length(tspan4)+length(tout)] .+= real(expect(S2,ψt))

    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)

    ρ_final = ptrace(ψt_end, [4,5])


    if ancilla1 == 0.0 && ancilla2 == 0.0
        ψ_target = α .* reduce(kron,[plus,plus,plus,a,a]) + (exp(1im * ϕ) * β) .* reduce(kron,[minus,minus,minus,a,a])
        ρ_target = ptrace(Ket(full_basis,ψ_target), [4,5])  
    else
        ψ_target = α .* reduce(kron,[plus,plus,plus,a,a]) - (exp(1im * ϕ) * β) .* reduce(kron,[minus,minus,minus,a,a])
        ρ_target = ptrace(Ket(full_basis,ψ_target), [4,5])
    end

    fidelity_after_correction = real(tr(sqrt(sqrt(ρ_final.data)*ρ_target.data*sqrt(ρ_final.data)))^2)
    println("Fidelity after correction is $(fidelity_after_correction)")


    print("Trajectory $N_trajectories Complete.\n")

    end_time = time() - start_time

    println("Simulation complete. Saving data...")
    # --- Save the data ---
    # full_data_folder = joinpath(data_folder, "phase_$(ϕ)")
    # if !isdir(full_data_folder)
    #     println("Directory does not exist. Creating directory...: $full_data_folder")
    #     mkpath(full_data_folder)
    # end
    # file_path = "$(full_data_folder)/N_atoms=$(total_qubits)_γ_dephase=$(γ_dephase)_phase=$(ϕ)_Ntraj=$(N_trajectories).jld2"
    # if isfile(file_path)
    #     println("File already exists. Overwriting the file: $file_path")
    # else
    #     println("Saving data to: $file_path")
    #     @save file_path population_data corrected_population_data final_step_population_data S1_data S2_data fidelity_after_encoding fidelity_after_correction end_time
    # end
    # @save "$(data_folder)/N_atoms=$(total_qubits)_γ_dephase=$(γ_dephase)_phase=$(ϕ)_Ntraj=$(N_trajectories).jld2" population_data corrected_population_data final_step_population_data S1_data S2_data fidelity_after_encoding fidelity_after_correction end_time 
    # @save "$(data_folder)/N_atoms=$(total_qubits)_γ_dephase=$(γ_dephase)_case1_Ntraj=$(N_trajectories).jld2" fidelity_after_encoding fidelity_after_storage fidelity_after_correction

    println("Data saved.")
end



# --- Parse command-line arguments ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    id = parse(Int, ARGS[1])
    N_trajectories = (id % 1000) + 1
    main(N_trajectories)
end
