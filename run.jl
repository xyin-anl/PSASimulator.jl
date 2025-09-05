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
Pkg.add("YAML")

using PSASimulator
using DataFrames
using CSV
using Dates
using YAML

println("✓ All modules loaded successfully")

# ===================================================================
# SELF-CONTAINED DATA AND HELPER FUNCTIONS
# ===================================================================

# --- DATA (Copied from demo/demo_data.jl) ---

const ISOTHERM_PARAMETERS = [
    7.022572e+00 6.326170e+00 2.647463e-12 2.269914e-11 -3.375370e+04 -3.641874e+04 7.022572e+00 0.000000e+00 8.228112e-09 0.000000e+00 -1.200000e+04 0.000000e+00 0.000000e+00;
    9.999929e+00 8.838786e+00 7.416992e-11 1.296440e-10 -2.493171e+04 -2.493171e+04 1.000000e+01 0.000000e+00 2.368616e-08 0.000000e+00 -5.800000e+03 0.000000e+00 0.000000e+00;
    1.000000e+00 1.892000e+01 1.395845e-09 5.911813e-12 -3.300000e+04 -3.300000e+04 3.200000e+00 0.000000e+00 1.040135e-08 0.000000e+00 -1.200000e+04 0.000000e+00 0.000000e+00;
    6.800000e+00 9.900000e+00 2.440000e-11 1.390000e-10 -4.200000e+04 -2.400000e+04 1.400000e+01 0.000000e+00 4.960000e-10 0.000000e+00 -1.800000e+04 0.000000e+00 0.000000e+00;
    4.800000e+01 0.000000e+00 8.060000e-10 0.000000e+00 -1.400000e+04 -1.400000e+04 4.800000e+01 0.000000e+00 1.110000e-09 0.000000e+00 -1.000000e+04 0.000000e+00 0.000000e+00;
    6.212000e+00 7.150000e+00 4.377195e-11 4.608793e-13 -3.785300e+04 -3.785300e+04 1.190000e+01 0.000000e+00 1.383811e-10 0.000000e+00 -1.943460e+04 0.000000e+00 0.000000e+00;
    2.855382e+01 0.000000e+00 5.938529e-11 0.000000e+00 -2.500000e+04 -2.500000e+04 2.855382e+01 0.000000e+00 1.220402e-09 0.000000e+00 -1.200000e+04 0.000000e+00 0.000000e+00;
    1.535652e+00 2.674846e+00 5.853452e-13 4.523258e-12 -3.250000e+04 -3.250000e+04 1.535652e+00 0.000000e+00 5.010345e-11 0.000000e+00 -1.800000e+04 0.000000e+00 0.000000e+00;
    6.816578e+00 0.000000e+00 8.436744e-11 0.000000e+00 -3.200000e+04 -3.200000e+04 3.890598e+00 0.000000e+00 3.827598e-09 0.000000e+00 -1.200000e+04 0.000000e+00 0.000000e+00;
    2.727278e+00 2.495993e-01 2.891278e-12 2.007587e-21 -5.209388e+04 -9.023652e+04 2.727278e+00 0.000000e+00 6.520000e-09 0.000000e+00 -1.200000e+04 0.000000e+00 0.000000e+00;
    4.540935e+00 0.000000e+00 9.931281e-13 0.000000e+00 -4.200000e+04 -4.200000e+04 2.181907e+00 0.000000e+00 1.538952e-10 0.000000e+00 -1.900000e+04 0.000000e+00 0.000000e+00;
    5.000000e+00 3.000000e+00 9.460000e-11 6.150000e-16 -3.300000e+04 -4.800000e+04 1.270000e+01 0.000000e+00 4.290000e-10 0.000000e+00 -1.230000e+04 0.000000e+00 0.000000e+00;
    7.265438e+00 0.000000e+00 1.215474e-09 0.000000e+00 -1.612044e+04 -1.612044e+04 7.265438e+00 0.000000e+00 2.345892e-10 0.000000e+00 -7.020028e+03 0.000000e+00 0.000000e+00;
    7.000000e+00 5.100000e+00 3.207140e-10 2.148947e-11 -2.750000e+04 -2.750000e+04 6.000000e+00 0.000000e+00 5.189774e-09 0.000000e+00 -1.000000e+04 0.000000e+00 0.000000e+00;
    6.570000e+00 3.130000e+00 1.440000e-07 9.410000e-07 -2.924300e+04 -3.080000e+04 9.606000e+00 0.000000e+00 7.350000e-06 0.000000e+00 -1.290000e+04 0.000000e+00 1.000000e+00;
    3.090000e+00 2.540000e+00 8.650000e-07 2.630000e-08 -3.664121e+04 -3.569066e+04 5.840000e+00 0.000000e+00 2.500000e-06 0.000000e+00 -1.580000e+04 0.000000e+00 1.000000e+00
];

const SIMULATION_PARAMETERS = [
    1.169000e+03 -3.450000e+04 -1.200000e+04;
    7.890000e+02 -2.490000e+04 -5.800000e+03;
    7.820000e+02 -3.300000e+04 -1.200000e+04;
    9.100000e+02 -3.130000e+04 -1.800000e+04;
    4.270000e+02 -1.400000e+04 -1.000000e+04;
    1.200000e+03 -3.785300e+04 -1.943460e+04;
    6.260000e+02 -2.500000e+04 -1.200000e+04;
    1.236000e+03 -3.250000e+04 -1.800000e+04;
    8.500000e+02 -3.200000e+04 -1.200000e+04;
    1.620000e+03 -5.530000e+04 -1.200000e+04;
    1.770000e+03 -4.200000e+04 -1.900000e+04;
    1.659000e+03 -3.300000e+04 -1.200000e+04;
    9.500000e+02 -1.610000e+04 -7.000000e+03;
    1.179000e+03 -2.750000e+04 -1.000000e+04;
    1.230000e+03 -3.030000e+04 -1.290000e+04;
    1.130000e+03 -3.600000e+04 -1.580000e+04
];

const MATERIALS_LIST = [
    (index=1, name="Co-MOF-74"), (index=2, name="Cu-BTTRi"), (index=3, name="Cu-TDPAT"),
    (index=4, name="Mg-MOF-74"), (index=5, name="MOF-177"), (index=6, name="Ni-MOF-74"),
    (index=7, name="NTU-105"), (index=8, name="Sc2BDC3"), (index=9, name="SIFSIX-2-Cu-i"),
    (index=10, name="SIFSIX-3-Ni"), (index=11, name="Ti-MIL-91"), (index=12, name="USTA-16"),
    (index=13, name="UiO-66(OH)2"), (index=14, name="ZIF-8"), (index=15, name="Zn-MOF-74"),
    (index=16, name="Zeolite_13X")
];

const OPT_VARS_PURITY = [
    1.000000e+05 8.082500e+02 2.400000e-01 2.600000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 6.843800e+02 1.400000e-01 1.000000e-01 1.000000e+00 1.000000e+04;
    1.580000e+05 5.806300e+02 1.800000e-01 2.000000e-01 1.000000e+00 1.000000e+04;
    2.050000e+05 5.110000e+02 1.400000e-01 2.700000e-01 1.000000e+00 1.000000e+04;
    1.600000e+05 9.138000e+01 1.900000e-01 1.000000e-01 1.000000e+00 1.000000e+04;
    1.080000e+05 8.660400e+02 1.400000e-01 2.800000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 4.740900e+02 1.900000e-01 1.000000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 2.471600e+02 1.600000e-01 1.000000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 5.192200e+02 1.800000e-01 3.600000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 4.343700e+02 1.100000e-01 1.000000e-01 1.000000e+00 1.000000e+04;
    1.120000e+05 5.268500e+02 1.800000e-01 3.500000e-01 1.000000e+00 1.000000e+04;
    1.420000e+05 8.283200e+02 1.100000e-01 3.500000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 6.921900e+02 2.000000e-01 3.600000e-01 1.000000e+00 1.000000e+04;
    1.160000e+05 2.012900e+02 1.100000e-01 1.000000e-01 1.000000e+00 1.000000e+04;
    1.000000e+05 8.599000e+02 2.200000e-01 3.000000e-01 1.000000e+00 1.000000e+04;
    1.640000e+05 4.865300e+02 1.000000e-01 3.000000e-01 1.000000e+00 1.000000e+04
];

# --- Helper functions to access the data ---

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
    else
        error("Scenario key '$(scenario_key)' not recognized.")
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
    fault_settings = config["fault_injection"]

    # Extract parameters from config
    N = sim_settings["N"]
    material_name = sim_settings["material_name"]
    scenario_name = sim_settings["scenario_name"]
    max_iterations = sim_settings["max_iterations"]
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
            println("\n--- FAULT INJECTION ---")
            for (key, value) in fault_settings
                println("  $(key): $(value)")
            end
            println("\n" * "="^60 * "\n")


            # --- 2. GET MATERIAL AND PROCESS VARIABLES ---
            material_properties, isotherm_params = get_material_data(material_name)
            material_data = (material_properties, isotherm_params)
            opt_vars = get_opt_vars(material_name, "MaxPurity90") # This is still hardcoded

            process_vars = [
                1.0,                                        # L [m]
                opt_vars[1],                                # P_0 [Pa]
                opt_vars[1] * opt_vars[4] / 8.314 / 313.15, # n_dot_0 [mol/s]
                opt_vars[2],                                # t_ads [s]
                opt_vars[3],                                # alpha [-]
                opt_vars[5],                                # beta [-]
                1.0e4,                                      # P_I [Pa]
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