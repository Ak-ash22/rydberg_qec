using JLD2

function average_populations()
    # Path to files
    base_path = "/scratch/roq68sum/shor_code_data/driving_abc12/"
    file_pattern = "N_atoms=3_γ_decay=0.001_Ntraj="

    num_files = 100  # Since Ntraj goes from 1 to 100

    # Initialize accumulators
    population_a_total = nothing
    population_c_total = nothing
    population_ac_total = nothing

    # Loop through each file
    for i in 1:num_files
        file_path = base_path * file_pattern * string(i) * ".jld2"
        
        # Load data
        @load file_path population_a population_c population_ac
        
        # Initialize accumulators on first iteration
        if isnothing(population_a_total)
            population_a_total = zeros(length(population_a))
            population_c_total = zeros(length(population_c))
            population_ac_total = zeros(length(population_ac))
        end

        # Sum the elements
        population_a_total .+= population_a
        population_c_total .+= population_c
        population_ac_total .+= population_ac
    end

    # Divide by 100 to get the average
    population_a_avg = population_a_total ./ num_files
    population_c_avg = population_c_total ./ num_files
    population_ac_avg = population_ac_total ./ num_files

    # Save final lists
    final_file_path = base_path * "N_atoms=3_γ_decay=0.001_Ntraj=1000_avg.jld2"
    @save final_file_path population_a_avg population_c_avg population_ac_avg

    println("Averaged populations saved to ", final_file_path)
end

# Call the function
average_populations()
