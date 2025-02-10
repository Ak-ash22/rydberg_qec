
#Saving the output
script_dir = "/home/agfleischhauer/roq68sum/rydberg_qec/codework"
# script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
data_folder = joinpath(script_dir, "shor_code_data/driving_abc_parallel")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

@everywhere begin
    include("functions.jl")   
end

# addprocs(5)

@everywhere function parallel_f(N::Int)
    """
    Main function to run the simulation
    """
    ψ0 = initialize_system()
    
    full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
    ψ0_ket = Ket(full_basis, ComplexF32.(ψ0)) 

    rng = MersenneTwister(N)
    rng_state = timeevolution.JumpRNGState(rng,0.5)
    @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f;maxiters=1e9,rng_state=rng_state)

    println("Trajectory $N.\n")
        
    return ψt
end


function main(N_trajectories::Int)

    println("Running the simulation with N_trajectories = $N_trajectories")
   

    @time begin

        basis = NLevelBasis(2)
        k = transition(basis,2,2)
        n_a = full_operator(k,total_qubits, [1])
        n_c = full_operator(k, total_qubits, [3])
        n_ac = full_operator(k, total_qubits, [1,3])

        population_a = zeros(length(tspan))
        population_c = zeros(length(tspan))
        population_ac = zeros(length(tspan))

        # Run trajectories in parallel using pmap (multiprocessing)
        ψ = pmap(parallel_f, 1:N_trajectories; batch_size=10)
        # Profile.print(format=:flat)

        m = length(tspan)
        l = N_trajectories
        ρ_avg = Vector{Matrix}(undef,m)
        for j in 1:m
            ρ_sum = zero(ψ[1][j].data * ψ[1][j].data')  # Initialize sum with a zero matrix of the same type
            for i in 1:l
                ρ_sum .+= (ψ[i][j].data * ψ[i][j].data')
            end
            ρ_avg[j] = ρ_sum / N_trajectories
        end

        for i in 1:length(tspan)
            population_a[i] = real(tr(n_a.data*ρ_avg[i]))
            population_c[i] = real(tr(n_c.data*ρ_avg[i]))
            population_ac[i] = real(tr(n_ac.data*ρ_avg[i]))
        end
    end

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/N_atoms=$(total_qubits)_γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories).jld2" ψ population_a population_c population_ac
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
