using JLD2, FileIO

function average_populations()
    # --- Path to data ---
    # base_path = "/scratch/roq68sum/shor_code_data/driving_abc9/"
    base_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/"
    file_pattern = "N_atoms=5_γ_dephase=0.0001_Ntraj="

    num_files = 1000  # Number of files to process

    # --- Load the first file to get available keys and array size dynamically ---
    first_file_path = base_path * file_pattern * "1.jld2"
    @load first_file_path population_data corrected_population_data #fidelity_data # Load dictionary from file

    num_timesteps = length(population_data[:a])  # Auto-detect array size
    # print(num_timesteps)

    # --- Initialize accumulators for all population types ---
    avg_population_data = Dict(key => zeros(num_timesteps) for key in keys(population_data))
    avg_corrected_population_data = Dict(key => zeros(length(corrected_population_data[:a])) for key in keys(corrected_population_data))

    print(length(population_data[:a]))
    print(length(avg_population_data[:a]))    # avg_fidelity_data = zeros(length(fidelity_data))
    print(type(population_data[:a]))
    print(type(avg_population_data[:a]))


    println("Processing $num_files files...")

    # --- Parallelized Loop for File Processing ---
    for i in 1:num_files
        file_path = base_path * file_pattern * string(i) * ".jld2"

        
        @load file_path population_data corrected_population_data #fidelity_data # Load the dictionary

        # Accumulate population data for all keys dynamically
        for key in keys(population_data)
            avg_population_data[key] .+= population_data[key]
        end

            # #Accumulate error corrected population data
            # for key in keys(corrected_population_data)
            #     avg_corrected_population_data[key] .+= corrected_population_data[key]
            # end

            # Accumulate fidelity data
            # avg_fidelity_data .+= fidelity_data

        # catch e
        #     @warn "Skipping missing or corrupted file in function 1: $file_path ($e)"
        # end
    end

    # --- Compute Averages ---
    for key in keys(avg_population_data)
        avg_population_data[key] .*= 1 / num_files
    end

    # for key in keys(avg_corrected_population_data)
    #     avg_corrected_population_data[key] .*= 1 / num_files
    # end

    # avg_fidelity_data .*= 1 / num_files

    # --- Save Averaged Data ---
    final_file_path = base_path * "N_atoms=5_γ_dephase=1e-4_Ntraj=$(num_files)_avg.jld2"
    @save final_file_path avg_population_data avg_corrected_population_data # avg_fidelity_data

    println("Averaged populations saved to $final_file_path")
end


function save_jump_files()
    base_path = "/scratch/roq68sum/5atoms_code/5_atom_correction/decay_1e_5/"
    file_pattern = "N_atoms=5_γ_decay=1.0e-5_Ntraj="
    num_files = 1000  # Number of files to process

    jump_folder = base_path * "jump_files/"

    if !isdir(jump_folder)
        println("Directory does not exist. Creating directory...: $jump_folder")
        mkpath(jump_folder)
    end

    for i in 1:num_files
        file_path = base_path * file_pattern * string(i) * ".jld2"
        try
            # Load the file
            @load file_path has_error detected_error has_correction_error population_data corrected_population_data fidelity_data end_time
            
            if has_error && detected_error
                # Save the jumps data to a new file
                jump_file_path = jump_folder * "jumps_corrected_" * string(i) * ".jld2"
                @save jump_file_path has_error population_data corrected_population_data fidelity_data end_time
            end

            if has_correction_error
                # Save the jumps data to a new file
                jump_file_path = jump_folder * "jumps_correction_error_" * string(i) * ".jld2"
                @save jump_file_path has_correction_error population_data corrected_population_data fidelity_data end_time
            end

            if !has_error && detected_error
                # Save the jumps data to a new file
                jump_file_path = jump_folder * "jumps_detection_error_" * string(i) * ".jld2"
                @save jump_file_path detected_error population_data corrected_population_data fidelity_data end_time
            end
                
            println("Saved jumps data Trajectories")

        catch e
            @warn "Skipping missing or corrupted file in function 2: $file_path ($e)"
        end
    end

end


# --- Run the function ---
average_populations()
# save_jump_files()