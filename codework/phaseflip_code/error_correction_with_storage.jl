include("functions.jl")

# Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/master_work/"
script_dir = "/scratch/roq68sum/5atoms_code/phaseflip_code/"
data_folder = joinpath(script_dir, "dephase_$(γ_dephase)/")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

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

    detected_error = false
    println("Starting the simulation...")
 
    ################################################################################################################ Protocol Starts

    println("Encoding Commencing...")
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan1,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)  # pass through to `solve`

    # println("Applying hadamards...")
    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψt_end,f2,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)
    plus = sqrt(1/2) .* a + sqrt(1/2) .* b
    minus = sqrt(1/2) .* a - sqrt(1/2) .* b
    ψ_target2 = α .* reduce(kron,[plus,plus,plus,a,a]) + (exp(1im * ϕ) * β) .* reduce(kron,[minus,minus,minus,a,a])
    fidelity_after_encoding = abs((dagger(Ket(full_basis,ψ_target2)) * ψt_end)^2)


    println("Phaseflip code encoding done successfully with fidelity $(fidelity_after_encoding).\n")

    ################################################################################################################ Storage of qubits

    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan6,ψt_end,f_storage,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    fidelity_after_storage = abs((dagger(Ket(full_basis,ψ_target2)) * ψt[end]/norm(ψt_end))^2)

    println("Fidelity after the storage time is $(fidelity_after_storage)")

    ################################################################################################################ Syndrome Measurement starts

    println("Applying hadamards...")
    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan2,ψt_end,f2,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    

    ψt_end = virtual_z_full * (ψt[end]/norm(ψt[end]))
    ψt_end /= norm(ψt_end)

    ψ_target3 = α .* reduce(kron,[b,b,b,a,a]) - (exp(1im * ϕ) * β) .* reduce(kron, [a,a,a,a,a])
    fidelity = abs((dagger(Ket(full_basis,ψ_target3)) * ψt_end)^2)

    println("Fidelity after hadamard is $(fidelity)")


    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan3,ψt_end,f3,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
    #Track Jumps Info
    has_error = length(jumps) > 0
    ancilla1_population = real(expect(n1_ancilla1, ψt))
    ancilla2_population = real(expect(n1_ancilla2,ψt))

    println("Ancillas drive complete.\n")

    ################################################################################################################ Error Detection
    println("Error Detection Commencing...")

    bA, bB, bC, b1, b2 = ψt[end].basis.bases
    ψt_end = normalize(ψt[end])

    rand_float1 = round(rand();digits=1)
    rand_float2 = round(rand();digits=1)

    if rand_float1 < round(ancilla1_population[end];digits=1) && rand_float2 < round(ancilla2_population[end];digits=1)
        println("Both Ancilla errors detected.")
        ancilla1 = 1.0
        ancilla2 = 1.0
        detected_error = true

        ancilla1_2 = projector(basisstate(b1,2))
        ancilla1_2 = Operator(ancilla1_2.basis_l, ancilla1_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_2.data))

        ancilla2_2 = projector(basisstate(b2,2))
        ancilla2_2 = Operator(ancilla2_2.basis_l, ancilla2_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_2.data))

        P11 = full_operator([ancilla1_2, ancilla2_2], 5, [4, 5])
        ψ_full = (P11 * ψt_end) / norm(P11 * ψt_end);

        ψ_proj_target = α .* reduce(kron,[b,a,b,b,b]) + (exp(1im * ϕ) * β) .* reduce(kron,[a,b,a,b,b])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)")

    elseif rand_float1 < round(ancilla2_population[end];digits=1)
        println("Ancilla 2 error detected.")
        ancilla1 = 0.0
        ancilla2 = 1.0
        detected_error = true

        ancilla1_1 = projector(basisstate(b1,1))
        ancilla1_1 = Operator(ancilla1_1.basis_l, ancilla1_1.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_1.data))

        ancilla2_2 = projector(basisstate(b2,2))
        ancilla2_2 = Operator(ancilla2_2.basis_l, ancilla2_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_2.data))

        P01 = full_operator([ancilla1_1, ancilla2_2], 5, [4, 5])
        ψ_full = (P01 * ψt_end) / norm(P01 * ψt_end);

        ψ_proj_target = α .* reduce(kron,[b,b,a,a,b]) + (exp(1im * ϕ) * β) .* reduce(kron,[a,a,b,a,b])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)")

    elseif rand_float1 < round(ancilla1_population[end];digits=1)
        println("Ancilla 1 error detected.")
        ancilla1 = 1.0
        ancilla2 = 0.0
        detected_error = true

        ancilla1_2 = projector(basisstate(b1,2))
        ancilla1_2 = Operator(ancilla1_2.basis_l, ancilla1_2.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla1_2.data))

        ancilla2_1 = projector(basisstate(b2,1))
        ancilla2_1 = Operator(ancilla2_1.basis_l, ancilla2_1.basis_r, SparseMatrixCSC{ComplexF32, Int64}(ancilla2_1.data))

        P10 = full_operator([ancilla1_2, ancilla2_1], 5, [4, 5])
        ψ_full = (P10 * ψt_end) / norm(P10 * ψt_end);

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

        ψ_proj_target = α .* reduce(kron,[b,b,b,a,a]) + (exp(1im * ϕ) * β) .* reduce(kron,[a,a,a,a,a])
        ψ_proj_target = Ket(full_basis,ψ_proj_target)
        ψ_proj_target /= norm(ψ_proj_target)
        fidelity_after_projection = abs((dagger(ψ_proj_target) * ψ_full)^2)
        print("Error Detection and Projection done with fidelity $(fidelity_after_projection)")
    end
    ################################################################################################################ Error Correction

    ψ1 = ψ_full

    #Correcting Atom B
    if ancilla1 == 1.0 && ancilla2 == 1.0

        println("Error Correction commencing for Atom B")

        fb = f_correct_factory(2)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fb,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0


    #Correcting Atom A
    elseif ancilla1 == 1.0

        println("Error Correction commencing for Atom A")

        fa = f_correct_factory(1)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fa,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0


    #Correcting Atom C
    elseif ancilla2 == 1.0

        println("Error Correction commencing for Atom C")
        
        fc = f_correct_factory(3)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,fc,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0



    #No Errors Detected
    else
        f0 = f_correct_factory(0)
        @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan4,ψ1,f0,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
        has_correction_error = length(jumps) > 0

    end

    println("Applying hadamards back...")
    ψt_end = ψt[end]/norm(ψt[end])
    ψt_end = virtual_z_full * (ψt_end)
    ψt_end /= norm(ψt_end)

    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan5,ψt_end,f5,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)
    
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

    println("Simulation complete in $(end_time). Saving data...")

    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_dephase=$(γ_dephase)_case1_Ntraj=$(N_trajectories).jld2" fidelity_after_encoding fidelity_after_storage fidelity_after_correction detected_error

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
