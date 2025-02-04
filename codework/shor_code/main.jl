include("dependencies.jl")
include("system_params.jl")
include("functions.jl")

#Saving the output
script_dir = "/home/agfleischhauer/roq68sum/rydberg_qec/codework"
# script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework"
data_folder = joinpath(script_dir, "shor_code_data/driving_abc_atoms")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end

const σ_x, n, Π_g, n, I, σ_minus, σ_plus, σ_z = system_constants()
const Ω, γ_Decay, γ_dephase, V_nn, δ, Δ1_0, T_optimal = unpack_params()

function main(N_trajectories::Int)
    """
    Main function to run the simulation
    """
    println("Running the simulation with N_trajectories = $N_trajectories")

    basis = NLevelBasis(2)
    k = transition(basis,2,2)
    n_a = full_operator(k,total_qubits, [1])
    n_c = full_operator(k, total_qubits, [3])
    n_ac = full_operator(k, total_qubits, [1,3])

    @time begin
        ψ0 = initialize_system()
        println("The system has been initialized.")
        # p1 = qubit_parameters(Ω, γ_Decay, γ_dephase, V_nn, δ)
        tspan = [0.0:1:T_optimal;]

        full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
        ψ0_ket = Ket(full_basis, ComplexF64.(ψ0)) 

        population_a = zeros(length(tspan))
        population_c = zeros(length(tspan))
        population_ac = zeros(length(tspan))

        println("Starting the simulation...")
        Threads.@threads for i in 1:N_trajectories
            
            @time tout, ψt = timeevolution.mcwf_dynamic(tspan,ψ0_ket,f;alg=Rodas3(autodiff=false),maxiters=1e7)
            println("Trajectory $i")
            population_a .+= real(expect(n_a, ψt))
            population_c .+= real(expect(n_c, ψt))
            population_ac .+= real(expect(n_ac, ψt))
        end

        population_a ./= N_trajectories
        population_c ./= N_trajectories
        population_ac ./= N_trajectories
    end

    println("Simulation complete. Saving data...")
    @save "$(data_folder)/γ_decay=$(γ_Decay)_Ntraj=$(N_trajectories)_4atoms.jld2" population_a population_c population_ac
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
