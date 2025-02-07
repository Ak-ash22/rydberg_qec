include("dependencies.jl")
include("system_params.jl")
include("functions.jl")

#Saving the output
script_dir = "/home/agfleischhauer/roq68sum/rydberg_qec/codework"
# script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
data_folder = joinpath(script_dir, "shor_code_data/benchmark_data")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """
    println("Running the simulation with N_trajectories = $N_trajectories")

    @time begin
        ψ0 = initialize_system()
        println("The system has been initialized.")
        
        ψ = Vector{Vector}(undef, N_trajectories)

        full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
        ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 

        println("Starting the simulation...")

        @sync Threads.@threads for i in 1:N_trajectories
            local ψt  # Local variable per thread
            @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f,dt=1e-2, abstol=1e-6, reltol=1e-4,maxiters=1e9,save_everystep=false)
            ψ[i] = ψt
            print("Trajectory $i/$N_trajectories.\n")
            GC.gc()  # Force garbage collection to prevent memory overflow
        end


        m = length(tspan)
        l = N_trajectories
        ρt = Matrix{Operator}(undef, l, m)


        println("Computing the density matrices...\n")
        @sync Threads.@threads for idx in 1:l*m
            i = (idx-1) ÷ m + 1  # Calculate the row index
            j = (idx-1) % m + 1  # Calculate the column index

            wavefunction = ψ[i][j] / norm(ψ[i][j])
            ρ = dm(wavefunction)

            ρt[i,j] = ρ
        end
        println("Density matrices computed.\n")


        ρ_avg = Vector{Operator}(undef, length(tspan))
       
        @sync Threads.@threads for j in 1:length(tspan)
            ρ_sum = zero(ρt[1, j])  # Initialize sum with a zero matrix of the same type
            for i in 1:N_trajectories
                ρ_sum += ρt[i, j]  # Sum all density matrices at time step j
            end
            ρ_avg[j] = ρ_sum / N_trajectories  # Compute the average
        end
        println("Averaged density matrix computation complete.\n")

    end


    println("Simulation complete. Saving data...")
    @save "$(data_folder)/mcwf_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" ρ_avg 
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
