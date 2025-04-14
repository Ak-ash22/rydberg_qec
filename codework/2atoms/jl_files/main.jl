include("functions.jl")

#Saving the output
# script_dir = "/home/agfleischhauer/roq68sum/rydberg_qec/codework"
script_dir = "C:/Users/14aka/OneDrive/Documents/rydberg_qec/codework/$(n_atoms)atoms"
data_folder = joinpath(script_dir, "results_data")

if !isdir(data_folder)
    println("Directory does not exist. Creating directory...: $data_folder")
    mkpath(data_folder)
end


const Ω1, Ω2, γ_Decay, γ_dephase, V1_nn, V2_nn, δ, Δ1_0, Δ2_0, T_optimal = unpack_params()


function case1()
    """
    Case 1: Optimal fidelity with respect to decay and sweep rate
            - Returns the end fidelity for a range of decay and sweep rate
            - And also returns the maximum fidelity obtained for a range of decay and sweep rate

    """
    println("Running for Optimal fidelity with respect to decay and sweep rate")


    @time begin
        decay = 10 .^ range(-5,-1,length=5);
        sweep_rate = collect(range(0.01,0.50,length=1000));
        l = length(decay)
        m = length(sweep_rate)
    
        ψ_ideal = α*kron(r,r,r) + β*kron(g,g,g)
        fidelity_end = Array{Float64}(undef,l,m)  # Ensure fidelity is pre-allocated
        max_fidelity = Array{Float64}(undef,l,m)  # Ensure max fidelity is pre-allocated
        Threads.@threads for idx in 1:l*m
            i = (idx-1) ÷ m + 1  # Calculate the row index
            j = (idx-1) % m + 1  # Calculate the column index

            a = decay[i]
            b = sweep_rate[j]

            p = Parameters(Ω1,Ω2,a,γ_dephase,V1_nn,V2_nn,b)
            T_optimal = 2*(Δ1_0+V1_nn)/b
            tspan = (0.0, T_optimal)

            f = []
            @time solution = solve_master_eqn(p, tspan)

            for ρ in solution.u
                push!(f, tr(ψ_ideal' * ρ * ψ_ideal))
            end
            fidelity_end[i, j] = real(f[end])
            max_fidelity[i, j] = maximum(real(f))
        end
    end
 
    @save "$(data_folder)/decay_endmax_fidelity.jld2" decay sweep_rate fidelity_end max_fidelity
end


function case2()
    """
    Case 2: Optimal fidelity with respect to dephase and sweep rate
            - Returns the end fidelity for a range of dephase and sweep rate
            - And also returns the maximum fidelity obtained for a range of dephase and sweep rate

    """
    println("Running for Optimal fidelity with respect to dephase and sweep rate")


    @time begin
        dephase = 10 .^ range(-5,-1,length=5);
        sweep_rate = collect(range(0.01,0.50,length=1000));
        l = length(dephase)
        m = length(sweep_rate)
    
        ψ_ideal = α*kron(r,r,r) + β*kron(g,g,g)
        fidelity_end = Array{Float64}(undef,l,m)  # Ensure fidelity is pre-allocated
        max_fidelity = Array{Float64}(undef,l,m)  # Ensure max fidelity is pre-allocated
        Threads.@threads for idx in 1:l*m
            i = (idx-1) ÷ m + 1  # Calculate the row index
            j = (idx-1) % m + 1  # Calculate the column index

            a = dephase[i]
            b = sweep_rate[j]

            p = Parameters(Ω1,Ω2,γ_Decay,a,V1_nn,V2_nn,b)
            T_optimal = 2*(Δ1_0+V1_nn)/b
            tspan = (0.0, T_optimal)

            f = []
            @time solution = solve_master_eqn(p, tspan)

            for ρ in solution.u
                push!(f, tr(ψ_ideal' * ρ * ψ_ideal))
            end
            fidelity_end[i, j] = real(f[end])
            max_fidelity[i, j] = maximum(real(f))
        end
    end
 
    @save "$(data_folder)/dephase_endmax_fidelity.jld2" dephase sweep_rate fidelity_end max_fidelity
end

function main()
    println("Choose a case to run: ")
    println("1: Optimal fidelity with respect to decay and sweep rate")
    println("2: Optimal fidelity with respect to dephase and sweep rate")
    println("Enter the case number: ")

    choice = readline()

    # try
    choice = parse(Int64, choice)
    if choice == 1
        case1()
    elseif choice == 2
        case2()
    else
        println("Bruh! Enter a valid choice")
    end
    # catch e
    #     println("An Error has occured in the case selected")
    # end
end

main()