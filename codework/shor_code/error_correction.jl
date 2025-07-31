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
    ϕ = 0

    println("Running the simulation with trajectory number = $(N_trajectories)")
    
    i = N_trajectories
    start_time = time()
    ψ0 = initialize_system()
    println("The system has been initialized.")

    full_basis = CompositeBasis([NLevelBasis(3) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 
    
    println("Starting the simulation...")
 
    ################################################################################################################ Protocol Starts

    println("Encoding Commencing...")
    @time tout, ψt, jumps = timeevolution.mcwf_dynamic(tspan1,ψ0_ket,f1,maxiters=1e9,seed=(N_trajectories*10000 + i),display_jumps=true)  # pass through to `solve`
    
    println("Encoding done! Calculating Fidelity...")

    plus = sqrt(1/2) .* a + sqrt(1/2) .* b
    minus = sqrt(1/2) .* a - sqrt(1/2) .* b
    ψ_target2 = α .* get_full_wavefunction([plus,plus,plus],[1,2,3]) + (exp(1im * ϕ) * β) .* get_full_wavefunction([minus,minus,minus],[1,2,3])
    fidelity_after_encoding = abs((dagger(Ket(full_basis,ψ_target2)) * ψt_end)^2)
    end_time = time()

    println("Phaseflip code encoding done successfully with fidelity $(fidelity_after_encoding).\n")
    println("Time taken for Encoding Step is $(end_time-start_time)")

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
