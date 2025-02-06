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

    ρt = []

    @time begin
        ψ0 = initialize_system()
        println("The system has been initialized.")
  
        full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
        ψ0_ket = Ket(full_basis, ComplexF64.(ψ0)) 

        println("Starting the simulation...")
        Threads.@threads for i in 1:N_trajectories

            @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f;alg=QNDF(autodiff=false),abstol=1e-8,reltol=1e-6,maxiters=1e7)

            for j in 1:length(tout)
                ψt[j] = ψt[j] / norm(ψt[j])
                if i==1
                    push!(ρt, dm(ψt[j]))
                else
                    ρt[j] .+= dm(ψt[j])
                end
            end
        print("Trajectory $i/$N_trajectories.")
        end

    for i in 1:length(tspan)
        ρt[i] ./= N_trajectories
    end

    end

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/mcwf_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" ρt
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
