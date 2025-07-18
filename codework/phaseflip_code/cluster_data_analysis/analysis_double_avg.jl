using JLD2

dephase_list = [0.1, 0.01, 0.001, 0.0001, 1.0e-5]
case_list = [1,2,3,4]

function average_population_again(dephase,case)
    base_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/avg_analysis/dephase_$(dephase)/case$(case)/"
    file_pattern = "N_atoms=5_γ_dephase=$(dephase)_"


    phase_list = [0.915, 3.234, 2.458, 2.769, 5.551, 0.359, 1.227, 4.389, 3.466, 1.899]
    num_files = 10

    first_file_path = base_path * file_pattern * "phase_$(phase_list[1])_Ntraj=1000_avg.jld2"
    # @load first_file_path population_data corrected_population_data #fidelity_data # Load dictionary from file
    @load first_file_path population_data corrected_population_data final_step_population_data S1_data S2_data fidelity_after_encoding fidelity_after_correction

    # --- Initialize accumulators for all population types ---
    avg_population_data = Dict(key => copy(population_data[key]) for key in (:a,:b,:c,:a1,:a2))
    avg_corrected_population_data = Dict(key => copy(corrected_population_data[key]) for key in (:a,:b,:c,:a1,:a2))
    avg_final_step_population_data = Dict(key => copy(final_step_population_data[key]) for key in (:a,:b,:c))
    avg_S1_data = copy(S1_data)
    avg_S2_data = copy(S2_data)
    
    avg_fidelity_after_encoding = fidelity_after_encoding
    avg_fidelity_after_correction = fidelity_after_correction


    println("Processing $num_files files...")

    # --- Parallelized Loop for File Processing ---
    for i in 2:num_files
        file_path = base_path * file_pattern * "phase_$(phase_list[i])_Ntraj=1000.jld2"

        
        @load file_path population_data corrected_population_data final_step_population_data S1_data S2_data fidelity_after_encoding fidelity_after_correction   #fidelity_data # Load the dictionary

        # Accumulate population data and corrrected populaiton data for all keys dynamically
	    for key in keys(avg_population_data)
	        avg_population_data[key] .+= population_data[key]	
            avg_corrected_population_data[key] .+= corrected_population_data[key]	
        end

        #Accumulate final step population data
        for key in keys(avg_final_step_population_data)
    	    avg_final_step_population_data[key] .+= final_step_population_data[key]
    	end

        avg_S1_data .+= S1_data
        avg_S2_data .+= S2_data

        avg_fidelity_after_encoding += fidelity_after_encoding
        avg_fidelity_after_correction += fidelity_after_correction
    end

    # --- Compute Averages ---
    for key in keys(avg_population_data)
        avg_population_data[key] .*= 1 / num_files
        avg_corrected_population_data[key] .*= 1 / num_files
    end

    for key in keys(avg_final_step_population_data)
        avg_final_step_population_data[key] .*= 1 / num_files
    end

    avg_S1_data .*= 1 / num_files
    avg_S2_data .*= 1 / num_files

    avg_fidelity_after_encoding *= 1 / num_files
    avg_fidelity_after_correction *= 1 / num_files

    # --- Save Averaged Data ---
    final_file_path = "/scratch/roq68sum/5atoms_code/phaseflip_code/avg_analysis/dephase_$(dephase)/case$(case)/N_atoms=5_γ_dephase=$(dephase)_Nphase=$(num_files)_avg.jld2"
    @save final_file_path avg_population_data avg_corrected_population_data avg_final_step_population_data avg_S1_data avg_S2_data avg_fidelity_after_encoding avg_fidelity_after_correction # avg_fidelity_data

    println("Average of Averaged populations saved to $final_file_path")

end



if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia main.jl <N_trajectories>")
        exit(1)
    end
    id = parse(Int, ARGS[1])
    dephase = dephase_list[fld(id,4)+1]
    case = case_list[(id % 4) + 1]
    average_populations_again(dephase, case)
end

#########Run the slurm batch over 20 jobs############
