# ===================================================================
# PSA SIMULATOR - STANDALONE RUN SCRIPT
# ===================================================================
#
# This script runs a single PSA simulation based on settings
# defined in the accompanying `config.yaml` file.
#
# ===================================================================

println("\n" * "="^60)
println(" PSA SIMULATOR - CONFIGURABLE RUN ")
println("="^60 * "\n")

# ===================================================================
# SETUP AND IMPORTS
# ===================================================================

using Pkg
# Activate the project environment and add YAML dependency if needed
Pkg.activate(@__DIR__)

using PSASimulator
using DataFrames
using CSV
using Dates
using YAML

# Include the configuration data module
include("config.jl")
using .PSAConfigData

println("✓ All modules loaded successfully")

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

function get_material_index(material_name)
    for material in MATERIALS_LIST
        if material.name == material_name
            return material.index
        end
    end
    error("Material '$(material_name)' not found in the dataset.")
end

function get_material_data(material_name)
    idx = get_material_index(material_name)
    properties = SIMULATION_PARAMETERS[idx, :]
    isotherm = ISOTHERM_PARAMETERS[idx, :]
    return properties, isotherm
end

function get_opt_vars(material_name, scenario_key)
    idx = get_material_index(material_name)
    if scenario_key == "MaxPurity90"
        return OPT_VARS_PURITY[idx, :]
    elseif scenario_key == "MaxPurity95"
        return OPT_VARS_RECOVERY[idx, :]
    elseif scenario_key == "MaxProductivity"
        return OPT_VARS_PRODUCTIVITY[idx, :]
    elseif scenario_key == "MinEnergy"
        return OPT_VARS_ENERGY[idx, :]
    else
        error("Optimization scenario key '$(scenario_key)' not recognized.")
    end
end

function save_simulation_data(traj, material_name, scenario_name, N, timestamp_dir)
    if traj === nothing
        println("⚠️ No trajectory data to save.")
        return
    end
    println("  Saving data to: $(timestamp_dir)")

    headers = ["Time"]
    append!(headers, ["P_Node$(i)" for i in 1:N+2])
    append!(headers, ["y_Node$(i)" for i in 1:N+2])
    append!(headers, ["x1_Node$(i)" for i in 1:N+2])
    append!(headers, ["x2_Node$(i)" for i in 1:N+2])
    append!(headers, ["T_Node$(i)" for i in 1:N+2])

    step_map = [
        ("a_storage", "t1_storage", "1_Co-current_Pressurization"),
        ("b_storage", "t2_storage", "2_Adsorption"),
        ("c_storage", "t3_storage", "3_Heavy_Reflux"),
        ("d_storage", "t4_storage", "4_Counter-current_Depressurization"),
        ("e_storage", "t5_storage", "5_Light_Reflux")
    ]

    if !haskey(traj, :a_storage) || isempty(traj[:a_storage])
        println("⚠️ Trajectory history not found in results. Nothing to save.")
        return
    end

    num_cycles = length(traj[:a_storage])
    println("  Found data for $(num_cycles) cycles.")

    for cycle_idx in 1:num_cycles
        cycle_dir = joinpath(timestamp_dir, "cycle_$(cycle_idx)")
        mkpath(cycle_dir)

        for (data_key, time_key, step_name) in step_map
            data_storage = traj[Symbol(data_key)]
            time_storage = traj[Symbol(time_key)]

            if length(data_storage) >= cycle_idx && length(time_storage) >= cycle_idx
                data = data_storage[cycle_idx]
                time = time_storage[cycle_idx]
                
                df = DataFrame(hcat(time, data), headers)
                filename = "$(step_name).csv"
                CSV.write(joinpath(cycle_dir, filename), df)
            end
        end
    end
end

# ===================================================================
# MAIN EXECUTION
# ===================================================================

function execute_simulation()
    # --- 1. LOAD CONFIGURATION ---
    config = YAML.load_file("config.yaml")
    sim_settings = config["simulation_settings"]
    proc_vars_config = config["process_variables"]
    fault_settings = config["fault_injection"]

    # Extract parameters from config
    N = sim_settings["N"]
    material_name = sim_settings["material_name"]
    scenario_name = sim_settings["scenario_name"]
    max_iterations = sim_settings["max_iterations"]
    optimization_scenario = sim_settings["optimization_scenario"]
    run_type = "ProcessEvaluation" # This could also be in the config

    # --- Create Timestamped Directory ---
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    base_dir = "simulation_output"
    material_dir = joinpath(base_dir, replace(material_name, r"[/: ]" => "_"))
    scenario_dir = joinpath(material_dir, replace(scenario_name, r"[/: ]" => "_"))
    timestamp_dir = joinpath(scenario_dir, timestamp)
    mkpath(timestamp_dir)

    log_file_path = joinpath(timestamp_dir, "simulation.log")
    open(log_file_path, "w") do log_file_stream
        original_stdout = stdout
        redirect_stdout(log_file_stream)
        try
            println("="^60)
            println(" PSA SIMULATOR RUN ")
            println("="^60)
            println("Timestamp: $(timestamp)")
            println("\n--- SIMULATION CONFIGURATION ---")
            for (key, value) in sim_settings
                println("  $(key): $(value)")
            end
            println("\n--- PROCESS VARIABLES ---")
            for (key, value) in proc_vars_config
                println("  $(key): $(value)")
            end
            println("\n--- FAULT INJECTION ---")
            for (key, value) in fault_settings
                println("  $(key): $(value)")
            end
            println("\n" * "="^60 * "\n")


            # --- 2. GET MATERIAL AND PROCESS VARIABLES ---
            material_properties, isotherm_params = get_material_data(material_name)
            material_data = (material_properties, isotherm_params)
            opt_vars = get_opt_vars(material_name, optimization_scenario)

            process_vars = [
                proc_vars_config["bed_length"],            # L [m]
                opt_vars[1],                                # P_0 [Pa]
                opt_vars[1] * opt_vars[4] / 8.314 / 313.15, # n_dot_0 [mol/s]
                opt_vars[2],                                # t_ads [s]
                opt_vars[3],                                # alpha [-]
                opt_vars[5],                                # beta [-]
                proc_vars_config["intermediate_pressure"],# P_I [Pa]
                opt_vars[6]                                 # P_l [Pa]
            ]

            # --- 3. RUN THE SIMULATION ---
            println("🚀 Running simulation... (This may take a moment)")
            result = PSASimulator.psacycle(process_vars, material_data; N=N, run_type=Symbol(run_type), it_disp=true, max_iters=max_iterations)

            # --- 4. SAVE THE RESULTS ---
            if result.traj !== nothing
                redirect_stdout(original_stdout) # Switch back to console for user messages
                println("\n💾 Simulation finished. Saving data...")
                save_simulation_data(result.traj, material_name, scenario_name, N, timestamp_dir)
                println("\n" * "="^60)
                println("✅ SIMULATION COMPLETE & DATA SAVED")
                println("="^60 * "\n")
            else
                redirect_stdout(original_stdout)
                println("❌ SIMULATION FAILED. No data to save.")
            end
        catch e
            redirect_stdout(original_stdout)
            println(stderr, "An error occurred during simulation: ", e)
            showerror(stderr, e, catch_backtrace())
        finally
            redirect_stdout(original_stdout) # Ensure it's always restored
        end
    end
    println("✓ Simulation log saved to: $(log_file_path)")
end

# --- RUN THE SCRIPT ---
execute_simulation()
