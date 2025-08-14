using JLD2, FileIO

function average_populations(s)
    # --- Path to data ---
    # base_path = "/scratch/roq68sum/shor_code_data/driving_abc9/"
    # base_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/dephase_1.0e-5/"
    dephase = 1.0e-5
    ϕ = 0.915

    base_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/dephase_$(dephase)/s$(s)/"
    file_pattern = "N_atoms=5_γ_dephase=$(dephase)_case1_Ntraj="

    num_files = 1000  # Number of files to process

    # --- Load the first file to get available keys and array size dynamically ---
    first_file_path = base_path * file_pattern * "1.jld2"

    @load first_file_path fidelity_after_encoding fidelity_after_storage fidelity_after_correction detected_error

    # --- Initialize accumulators for all population types ---
    avg_fidelity_after_encoding = fidelity_after_encoding
    avg_fidelity_after_storage = fidelity_after_storage
    avg_fidelity_after_correction = fidelity_after_correction


    println("Processing $num_files files...")

    # --- Parallelized Loop for File Processing ---
    fid_after_enc = []
    fid_after_store = []
    fid_after_correct = []

    for i in 2:num_files
        file_path = base_path * file_pattern * string(i) * ".jld2"

        @load file_path fidelity_after_encoding fidelity_after_storage fidelity_after_correction detected_error
        
        avg_fidelity_after_encoding += fidelity_after_encoding
        avg_fidelity_after_storage += fidelity_after_storage
        avg_fidelity_after_correction += fidelity_after_correction

        push!(fid_after_enc, fidelity_after_encoding)
        push!(fid_after_store, fidelity_after_storage)
        push!(fid_after_correct, fidelity_after_correction)
    end

    # --- Compute Averages ---
    avg_fidelity_after_encoding *= 1 / num_files
    avg_fidelity_after_storage *= 1 / num_files
    avg_fidelity_after_correction *= 1 / num_files

    # --- Save Averaged Data ---
    final_file_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/dephase_$(dephase)/N_atoms=5_γ_dephase=$(dephase)_case1_phase=$(ϕ)_Ntraj=$(num_files)_avg_s$(s).jld2"
    @save final_file_path avg_fidelity_after_encoding avg_fidelity_after_storage avg_fidelity_after_correction fid_after_enc fid_after_store fid_after_correct

    println("Averaged populations saved to $final_file_path")
end


function save_jump_files(s)
    dephase = 1.0e-5
    ϕ = 0.915

    base_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/dephase_$(dephase)/s$(s)/"
    file_pattern = "N_atoms=5_γ_dephase=$(dephase)_case1_Ntraj="

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
            @load file_path fidelity_after_encoding fidelity_after_storage fidelity_after_correction detected_error
            
            if detected_error
                # Save the jumps data to a new file
                jump_file_path = jump_folder * "jumps_detected_" * string(i) * ".jld2"
                @save jump_file_path fidelity_after_encoding fidelity_after_storage fidelity_after_correction detected_error
            end
            println("Saved jumps data Trajectories")

        catch e
            @warn "Skipping missing or corrupted file in function 2: $file_path ($e)"
        end
    end

end







# --- Run the function ---
# --- Parse command-line arguments ---
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    id = parse(Int, ARGS[1])
    average_populations(id)
    save_jump_files(id)
end
#######Run slurm batch over 200 jobs#############
