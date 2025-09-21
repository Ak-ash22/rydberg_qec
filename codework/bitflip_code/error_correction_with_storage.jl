include("functions.jl")

#Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/"
# # script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
script_dir = "/scratch/roq68sum/5atoms_code/bitflip_code/"

function main(N_trajectories::Int,s)
    """
    Main function to run the simulation
    """
    data_folder = joinpath(script_dir, "decay_$(γ_Decay)/s$(s)")

    if !isdir(data_folder)
        println("Directory does not exist. Creating directory...: $data_folder")
        mkpath(data_folder)
    end

    T_storage = s*798
    tspan_s = [0.0: 0.1: T_storage;]
    
    println("Running the simulation with N_trajectories = $N_trajectories")

    i = N_trajectories
    start_time = time()
    ψ0 = initialize_system()
    println("The system has been initialized.")

    full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 

    ############################################# --- Preallocate Arrays ---
    num_timesteps = length(tspan1)+length(tspan2)+length(tspan_s)
    population_data = Dict(key => zeros(num_timesteps) for key in (:a, :b,  :c, :a1, :a2, :abc))
    # print(num_timesteps)
    # print(length(population_data[:a]))
    
    max_length = max(length(tspan3), length(tspan4))  # Choose longest possible time span
    fidelity_data = zeros(max_length)  # Initialize with zeros

    corrected_population_data = Dict(key => zeros(max_length) for key in (:a, :b,  :c, :a1, :a2))
    has_error = false
    detected_error = false
    has_correction_error = false
    
    ################################################################################################################ Protocol Starts
  
    println("Starting the simulation. Encoding commencing...")
    
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan1,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*10000),display_jumps=true)

    #Track Jumps Info
    has_error = length(jumps) > 0

    ancilla1_population = real(expect(n_1, ψt))
    ancilla2_population = real(expect(n_2, ψt))

    population_data[:a][1:length(tspan1)] .+= real(expect(n_a, ψt))
    population_data[:b][1:length(tspan1)] .+= real(expect(n_b, ψt))
    population_data[:c][1:length(tspan1)] .+= real(expect(n_c, ψt))
    population_data[:abc][1:length(tspan1)] .+= real(expect(n_abc, ψt))
    population_data[:a1][1:length(tspan1)] .+= ancilla1_population
    population_data[:a2][1:length(tspan1)] .+= ancilla2_population

    ψt_end = ψt[end]/norm(ψt[end])
    ψ_target1 = α .* reduce(kron,[r,r,r,g,g]) + β .* reduce(kron,[g,g,g,g,g])
    fidelity_after_encoding = abs((dagger(Ket(full_basis,ψ_target1)) * ψt_end)^2)
    println("Bitflip code encoding done successfully with fidelity $(fidelity_after_encoding)")


    ################################################################################################################ Storage Starts
    println("Storing the data qubits...")
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan_s,ψt_end,f_storage,maxiters=1e9,seed=(N_trajectories*10000),display_jumps=true)

    #Track Jumps Info
    has_storage_error = length(jumps) > 0

    ancilla1_population = real(expect(n_1, ψt))
    ancilla2_population = real(expect(n_2, ψt))

    population_data[:a][length(tspan1)+1:length(tspan1)+length(tspan_s)] .+= real(expect(n_a, ψt))
    population_data[:b][length(tspan1)+1:length(tspan1)+length(tspan_s)] .+= real(expect(n_b, ψt))
    population_data[:c][length(tspan1)+1:length(tspan1)+length(tspan_s)] .+= real(expect(n_c, ψt))
    population_data[:abc][length(tspan1)+1:length(tspan1)+length(tspan_s)] .+= real(expect(n_abc, ψt))
    population_data[:a1][length(tspan1)+1:length(tspan1)+length(tspan_s)] .+= ancilla1_population
    population_data[:a2][length(tspan1)+1:length(tspan1)+length(tspan_s)] .+= ancilla2_population


    ψt_end = ψt[end]/norm(ψt[end])
    ψ_target2 = α .* reduce(kron,[r,r,r,g,g]) + β .* reduce(kron,[g,g,g,g,g])
    fidelity_after_storage = abs((dagger(Ket(full_basis,ψ_target2)) * ψt_end)^2)
    println("Fidelity after storage is $(fidelity_after_storage)")

    ################################################################################################################ Driving ancillas 1,2,3
    println("Storage done. Driving Ancillas...")

    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψt_end,f2,maxiters=1e9,seed=(N_trajectories*10000),display_jumps=true)

    #Track Jumps Info
    has_error = length(jumps) > 0

    ancilla1_population = real(expect(n_1, ψt))
    ancilla2_population = real(expect(n_2, ψt))

    population_data[:a][length(tspan1)+length(tspan_s)+1:length(tspan1)+length(tspan_s)+length(tspan2)] .+= real(expect(n_a, ψt))
    population_data[:b][length(tspan1)+length(tspan_s)+1:length(tspan1)+length(tspan_s)+length(tspan2)] .+= real(expect(n_b, ψt))
    population_data[:c][length(tspan1)+length(tspan_s)+1:length(tspan1)+length(tspan_s)+length(tspan2)] .+= real(expect(n_c, ψt))
    population_data[:abc][length(tspan1)+length(tspan_s)+1:length(tspan1)+length(tspan_s)+length(tspan2)] .+= real(expect(n_abc, ψt))
    population_data[:a1][length(tspan1)+length(tspan_s)+1:length(tspan1)+length(tspan_s)+length(tspan2)] .+= ancilla1_population
    population_data[:a2][length(tspan1)+length(tspan_s)+1:length(tspan1)+length(tspan_s)+length(tspan2)] .+= ancilla2_population

   ################################################################################################################ Error Detection Starts

    println("Error Detection Commencing...")

    bA, bB, bC, b1, b2 = ψt[end].basis.bases
    ψt_end = normalize(ψt[end])

    rng = MersenneTwister(N_trajectories*10000 + i)
    rand_float1 = rand(rng,1)[1]
    rand_float2 = rand(rng,1)[1]

    if rand_float1 < ancilla1_population[end] && rand_float2 < ancilla2_population[end]
        println("Both Ancilla errors detected. Correcting Atom B")
        ancilla1 = 1.0
        ancilla2 = 1.0
        detected_error = true

        ancilla1_2 = projector(basisstate(b1,2))
        ancilla1_2 = Operator(ancilla1_2.basis_l, ancilla1_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_2.data))

        ancilla2_2 = projector(basisstate(b2,2))
        ancilla2_2 = Operator(ancilla2_2.basis_l, ancilla2_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_2.data))

        P11 = full_operator([ancilla1_2, ancilla2_2], 5, [4,5])
        ψ_full = (P11 * ψt_end) / norm(P11 * ψt_end);

        ψ_proj_target = α .* reduce(kron,[r,g,r,r,r]) + β .* reduce(kron,[g,r,g,r,r])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)\n")

        
    elseif rand_float1 < ancilla2_population[end]
        println("Ancilla 2 error detected. Correcting Atom C")
        ancilla1 = 0.0
        ancilla2 = 1.0
        detected_error = true

        ancilla1_1 = projector(basisstate(b1,1))
        ancilla1_1 = Operator(ancilla1_1.basis_l, ancilla1_1.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_1.data))

        ancilla2_2 = projector(basisstate(b2,2))
        ancilla2_2 = Operator(ancilla2_2.basis_l, ancilla2_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_2.data))

        P01 = full_operator([ancilla1_1, ancilla2_2], 5, [4, 5])
        ψ_full = (P01 * ψt_end) / norm(P01 * ψt_end);

        ψ_proj_target = α .* reduce(kron,[r,r,g,g,r]) + β .* reduce(kron,[g,g,r,g,r])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)\n")


    elseif rand_float1 < ancilla1_population[end]
        println("Ancilla 1 error detected. Correcting Atom A")
        ancilla1 = 1.0
        ancilla2 = 0.0
        detected_error = true

        ancilla1_2 = projector(basisstate(b1,2))
        ancilla1_2 = Operator(ancilla1_2.basis_l, ancilla1_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_2.data))

        ancilla2_1 = projector(basisstate(b2,1))
        ancilla2_1 = Operator(ancilla2_1.basis_l, ancilla2_1.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_1.data))

        P10 = full_operator([ancilla1_2, ancilla2_1], 5, [4, 5])
        ψ_full = (P10 * ψt_end) / norm(P10 * ψt_end);

        ψ_proj_target = α .* reduce(kron,[g,r,r,r,g]) + β .* reduce(kron,[r,g,g,r,g])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)\n")


    else 
        println("No errors detected.")
        ancilla1 = 0.0
        ancilla2 = 0.0

        ancilla1_1 = projector(basisstate(b1,1))
        ancilla1_1 = Operator(ancilla1_1.basis_l, ancilla1_1.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_1.data))

        ancilla2_1 = projector(basisstate(b2,1))
        ancilla2_1 = Operator(ancilla2_1.basis_l, ancilla2_1.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_1.data))

        P00 = full_operator([ancilla1_1, ancilla2_1], 5, [4, 5])
        ψ_full = (P00 * ψt_end) / norm(P00 * ψt_end);

        ψ_proj_target = α .* reduce(kron,[r,r,r,g,g]) + β .* reduce(kron,[g,g,g,g,g])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)\n")


    end
    ################################################################################################################ Error Correction

    ψ1 = ψ_full

    #Correcting Atom B
    if ancilla1 == 1.0 && ancilla2 == 1.0

        println("Error Correction Commencing for Atom B")

        f_correct2 = f_correct_factory(2)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,f_correct2,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0
    
        fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
        corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
        corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
        corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
        corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))

    #Correcting Atom A
    elseif ancilla1 == 1.0

        println("Error Correction Commencing for Atom A")

        f_correct1 = f_correct_factory(1)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψ1,f_correct1,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
        corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
        corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
        corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
        corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
    
    #Correcting Atom C
    elseif ancilla2 == 1.0

        println("Error Correction Commencing for Atom C")
        
        f_correct3 = f_correct_factory(3)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψ1,f_correct3,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
        corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
        corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
        corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
        corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
    
    
    #No Errors Detected
    else
        f_no_correct = f_correct_factory(0)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψ1,f_no_correct,maxiters=1e9,seed=(N_trajectories*1000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

        fidelity_data[1:length(tout)] .+= real(expect(n_abc, ψt))
        corrected_population_data[:a][1:length(tout)] .+= real(expect(n_a, ψt))
        corrected_population_data[:b][1:length(tout)] .+= real(expect(n_b, ψt))
        corrected_population_data[:c][1:length(tout)] .+= real(expect(n_c, ψt))
        corrected_population_data[:a1][1:length(tout)] .+= real(expect(n_1, ψt))
        corrected_population_data[:a2][1:length(tout)] .+= real(expect(n_2, ψt))
    end

    ψt_end = ψt[end]/norm(ψt[end])
    ρ_final = ptrace(ψt_end, [1,2])

    ψ_target = α .* reduce(kron,[r,r,r,g,g]) - β .* reduce(kron,[g,g,g,g,g])
    ρ_target = ptrace(Ket(full_basis,ψ_target), [1,2])
    fidelity_after_correction = real(tr(sqrt(sqrt(ρ_final.data)*ρ_target.data*sqrt(ρ_final.data)))^2)
    println("Fidelity after correction is $(fidelity_after_correction)")


    print("Trajectory $N_trajectories Complete.\n")

    end_time = time() - start_time

    println("Simulation complete in $(end_time). Saving data...")
    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" has_error detected_error has_correction_error population_data corrected_population_data fidelity_data end_time
    println("Data saved.")
end



# --- Parse command-line arguments ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    id = parse(Int, ARGS[1])
    N_trajectories = (id % 10000) + 1
    s = fld(id,10000)*10 + 10
    main(N_trajectories,s)
end

