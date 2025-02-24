include("functions.jl")

# Saving the output
# script_dir = "/scratch/roq68sum/shor_code_data"

# data_folder = joinpath(script_dir, "driving_abc")

# if !isdir(data_folder)
#     println("Directory does not exist. Creating directory...: $data_folder")
#     mkpath(data_folder)
# end

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """
    println("Running the simulation for $N_trajectories")

    start_time = time()

    basis = NLevelBasis(2)
    k = transition(basis,2,2)
    n_a = full_operator(k,total_qubits, [1])
    n_c = full_operator(k, total_qubits, [3])
    n_ac = full_operator(k, total_qubits, [1,3])

    ψ0 = initialize_system()
    println("The system has been initialized.")

    full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 

    println("Starting the simulation...")
    
    population_a = zeros(length(tspan))
    population_c = zeros(length(tspan)) 
    population_ac = zeros(length(tspan))

    ψ = Vector{Vector{Ket}}(undef,10)
    
    for i in 1:1
        @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f;maxiters=1e9,seed=(N_trajectories*100 + i))
        ψ[i] = ψt
        population_a .+= real(expect(n_a, ψt))
        population_c .+= real(expect(n_c, ψt))
        population_ac .+= real(expect(n_ac, ψt))
    end

    println("Simulation complete. Calculating the average wavefunction...")

    # ψ_avg = Vector{Ket}(undef,length(tspan))
    # for j in 1:length(tspan)
    #     ψ_sum = zero(ψ[1][j])  # Initialize sum with a zero matrix of the same type
    #     for i in 1:50
    #         ψ_sum .+= ψ[i][j]
    #     end
    #     ψ_avg[j] = ψ_sum ./ norm(ψ_sum)
    #     ψ_avg[j] = ψ_avg[j] ./ sqrt(50)
    # end


    # m = length(tspan)
    # l = 10
    # ρ_avg = Vector{Matrix}(undef,m)
    # for j in 1:m
    #     ρ_sum = zero(ψ[1][j].data * ψ[1][j].data')  # Initialize sum with a zero matrix of the same type
    #     for i in 1:l
    #         ρ_sum .+= (ψ[i][j].data * ψ[i][j].data')
    #     end
    #     ρ_avg[j] = ρ_sum / N_trajectories
    # end

    println("Trajectory $N_trajectories.\n")
    
    population_a .*= 1/10
    population_c .*= 1/10
    population_ac .*= 1/10

    end_time = time() - start_time

    println("Simulation complete in $(end_time). Saving data...")
    # @save "$(data_folder)/N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" population_a population_c population_ac ρ_avg end_time
    println("Data saved.")
end


# --- Parse command-line arguments ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    elseif length(ARGS) > 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    N_trajectories = parse(Int, ARGS[1])
    main(N_trajectories)
end
