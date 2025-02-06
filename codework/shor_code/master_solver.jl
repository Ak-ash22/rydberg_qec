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

function main()
    """
    Main function to run the simulation
    """
    ρt = []

    @time begin
        ψ0 = initialize_system()
        println("The system has been initialized.")
  
        full_basis = CompositeBasis([NLevelBasis(2) for _ in 1:total_qubits]...)
        ψ0_ket = Ket(full_basis, ComplexF64.(ψ0)) 

        println("Starting the simulation...")
        @time tout, ρt = timeevolution.master_dynamic(tspan,ψ0_ket,f;alg=QNDF(autodiff=false),maxiters=1e7)

        println("Simulation complete. Saving data...")
    end
    @save "$(data_folder)/master_γ_decay=$(γ_Decay).jld2" ρt
end

main()