using JLD2

function average_populations()
    # --- Path to data ---
    # base_path = "/scratch/roq68sum/shor_code_data/driving_abc9/"
    base_path = "/scratch/roq68sum/5atoms_code/5_atom_work"
    file_pattern = "N_atoms=5_γ_decay=0.001_Ntraj="

    num_files = 100  # Number of files to process

    # --- Load the first file to get available keys and array size dynamically ---
    first_file_path = base_path * file_pattern * "1.jld2"
    @load first_file_path population_data  # Load dictionary from file

    num_timesteps = length(population_data[:a])  # Auto-detect array size

    # --- Initialize accumulators for all population types ---
    avg_population_data = Dict(key => zeros(num_timesteps) for key in keys(population_data))

    println("Processing $num_files files...")

    # --- Parallelized Loop for File Processing ---
    for i in 1:num_files
        file_path = base_path * file_pattern * string(i) * ".jld2"

        try
            @load file_path population_data  # Load the dictionary

            # Accumulate population data for all keys dynamically
            for key in keys(population_data)
                avg_population_data[key] .+= population_data[key]
            end
        catch e
            @warn "Skipping missing or corrupted file: $file_path ($e)"
        end
    end

    # --- Compute Averages ---
    for key in keys(avg_population_data)
        avg_population_data[key] .*= 1 / num_files
    end

    # --- Save Averaged Data ---
    final_file_path = base_path * "N_atoms=5_γ_decay=0.001_Ntraj=$(num_files)_avg.jld2"
    @save final_file_path avg_population_data

    println("Averaged populations saved to $final_file_path")
end

# --- Run the function ---
average_populations()
