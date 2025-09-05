# Fault Simulation Guide for PSASimulator.jl

This document details the procedures used to simulate various fault conditions within the `PSASimulator.jl` environment. The goal is to provide a step-by-step guide for reproducing the fault data.

## 1. Prerequisites

Before starting, ensure the following:
*   You have `PSASimulator.jl` installed and its environment activated.
*   The `run.jl` script (which we created and committed) is present in the project root directory. This script is our primary tool for running simulations.
*   Familiarize yourself with `fault-injection.md` for conceptual understanding of fault types.

## 2. General Procedure for Fault Injection

Our approach involves using the `run.jl` script, which is self-contained. For certain fault types, we need to temporarily modify the core simulator file (`src/PSACycle.jl`).

**Important Considerations:**
*   **Temporary Modifications:** When `src/PSACycle.jl` is modified, it's crucial to **back up** the original file and **restore** it immediately after the simulation to avoid permanent changes to the simulator's source code.
*   **Output Organization:** All simulation results will be saved in the `simulation_output/` directory, organized by material, scenario, and a timestamp for each run.

## 3. Simulating Normal Operation (Baseline)

Before injecting faults, it's useful to generate a baseline "normal" operation data set.

**Procedure:**
1.  Ensure `run.jl` is in its default state (no fault injection code active, `scenario_name` set to `"90%_Recovery_(Max_Purity)"`).
    ```julia
    # Part of run.jl (ensure this section is as follows for normal run)
    material_name = "Zeolite_13X"
    scenario_name = "90%_Recovery_(Max_Purity)" # Default, non-faulty scenario name
    # ...
    # No fault injection code active in this section
    # material_properties, isotherm_params = get_material_data(material_name)
    # material_data = (material_properties, isotherm_params)
    ```
2.  Run the simulation:
    ```bash
    julia run.jl
    ```
**Expected Output:**
Data will be saved in `simulation_output/Zeolite_13X/90%_Recovery_(Max_Purity)/YYYY-MM-DD_HH-MM-SS/`.

## 4. Fault Type 1: Adsorbent Degradation

**Description:** This fault simulates a reduction in the adsorbent material's capacity over time, e.g., due to aging or poisoning. We will simulate a 20% loss in CO₂ adsorption capacity.

**Procedure:**
1.  **Modify `run.jl`:**
    *   Change the `scenario_name` to `"Adsorbent_Degradation"`.
    *   Inject code to reduce the CO₂ adsorption capacity parameters (`isotherm_params[1]` and `isotherm_params[2]`) by 20%.

    Locate the `execute_simulation` function in `run.jl` and modify the section where `material_data` is defined:

    ```julia
    # run.jl (modified section)
    function execute_simulation()
        # ...
        material_name = "Zeolite_13X"
        scenario_name = "Adsorbent_Degradation" # Changed for this fault
        run_type = "ProcessEvaluation"
        
        println("Setting up simulation for:")
        println("  - Material: $(material_name)")
        println("  - Scenario: $(scenario_name)\n")

        # --- 2. GET MATERIAL AND PROCESS VARIABLES ---
        material_properties, initial_isotherm_params = get_material_data(material_name)

        # --- FAULT INJECTION: ADSORBENT DEGRADATION ---
        println("⚠️ INJECTING FAULT: Adsorbent Degradation (20% capacity loss)")
        degraded_isotherm_params = copy(initial_isotherm_params)
        degradation_factor = 0.8 # 20% loss
        degraded_isotherm_params[1] *= degradation_factor # q_s_b for CO2
        degraded_isotherm_params[2] *= degradation_factor # q_s_d for CO2
        # --- END FAULT INJECTION ---

        material_data = (material_properties, degraded_isotherm_params)
        opt_vars = get_opt_vars(material_name, "MaxPurity90")
        # ... rest of the function ...
    end
    ```
2.  **Run the simulation:**
    ```bash
    julia run.jl
    ```
**Expected Output:**
Data will be saved in `simulation_output/Zeolite_13X/Adsorbent_Degradation/YYYY-MM-DD_HH-MM-SS/`.

**Cleanup:**
*   **Revert `run.jl`:** After generating the data, it's crucial to revert `run.jl` to its original state. You can do this by manually undoing the changes or by using a version control system if you've committed the original `run.jl`.

## 5. Fault Type 2: Stuck Valve (Heavy Reflux)

**Description:** This fault simulates a heavy reflux valve being stuck closed, preventing flow during that step. This requires a temporary modification to the simulator's core logic.

**Procedure:**
1.  **Ensure `run.jl` is in its normal state:** Make sure you have reverted any changes from previous fault simulations.
2.  **Backup `src/PSACycle.jl`:** Before making any changes, save a copy of the original file.
    ```bash
    cp src/PSACycle.jl src/PSACycle.jl.bak
    ```
3.  **Modify `src/PSACycle.jl`:**
    *   Locate the section where `ndot_HR` is defined for `Params[35]`. This occurs after the Light Reflux step (Step 5) and before the Heavy Reflux step (Step 3 in the next cycle).
    *   Inject code to set `Params[35]` to `0.0` (zero flow).

    Locate the `psacycle` function in `src/PSACycle.jl` and find the lines where `Params[33]`, `Params[34]`, `Params[35]` are updated (around line 330 in the original file). Modify this section:

    ```julia
    # src/PSACycle.jl (modified section)
    # ...
                # Update parameters for Heavy Reflux step
                y_HR = CO2Front / totalFront
                T_HR = TFront
                ndot_HR = totalFront * beta / t_HR
                Params[33] = y_HR
                Params[34] = T_HR
                Params[35] = ndot_HR

                # --- FAULT INJECTION: HEAVY REFLUX VALVE STUCK ---
                println("⚠️ INJECTING FAULT: Heavy Reflux Valve Stuck Closed")
                Params[35] = 0.0 # Set ndot_HR to zero
                # --- END FAULT INJECTION ---
    # ...
    ```
4.  **Modify `run.jl`:**
    *   Change the `scenario_name` to `"Stuck_Valve_Fault"`.

    ```julia
    # run.jl (modified section)
    function execute_simulation()
        # ...
        material_name = "Zeolite_13X"
        scenario_name = "Stuck_Valve_Fault" # Changed for this fault
        run_type = "ProcessEvaluation"
        # ...
    end
    ```
5.  **Run the simulation:**
    ```bash
    julia run.jl
    ```
**Expected Output:**
Data will be saved in `simulation_output/Zeolite_13X/Stuck_Valve_Fault/YYYY-MM-DD_HH-MM-SS/`. You will also see `⚠️ INJECTING FAULT: Heavy Reflux Valve Stuck Closed` messages during the simulation run.

**Cleanup:**
*   **Restore `src/PSACycle.jl`:** Immediately restore the original file from your backup.
    ```bash
    cp src/PSACycle.jl.bak src/PSACycle.jl
    rm src/PSACycle.jl.bak # Optional: remove the backup file
    ```
*   **Revert `run.jl`:** Change `scenario_name` back to `"90%_Recovery_(Max_Purity)"`.

## 6. Fault Type 3: Column Leak

**Description:** This fault simulates a leak in the column, causing a pressure drop. We will simulate a 5% pressure loss across all nodes after the Adsorption step.

**Procedure:**
1.  **Ensure `run.jl` is in its normal state:** Make sure you have reverted any changes from previous fault simulations.
2.  **Backup `src/PSACycle.jl`:** (If you removed the `.bak` file in the previous step, create it again).
    ```bash
    cp src/PSACycle.jl src/PSACycle.jl.bak
    ```
3.  **Modify `src/PSACycle.jl`:**
    *   Locate the section after the Adsorption step where the `b` state vector is calculated (`b = reduce(hcat, sol2.u)'`).
    *   Inject code to reduce the pressure values in `b` by 5%.

    Locate the `psacycle` function in `src/PSACycle.jl` and find the lines after the Adsorption step (around line 240 in the original file). Modify this section:

    ```julia
    # src/PSACycle.jl (modified section)
    # ...
                t2 = sol2.t
                b = reduce(hcat, sol2.u)'

                # --- FAULT INJECTION: SIMULATE A LEAK ---
                println("⚠️ INJECTING FAULT: Column Leak (5% pressure loss)")
                leakage_factor = 0.95
                b[:, 1:N+2] .*= leakage_factor 
                # --- END FAULT INJECTION ---

                # Correct the output (clean up results from simulation)
                idx = findall(b[:, N+1] .< 1)
    # ...
    ```
4.  **Modify `run.jl`:**
    *   Change the `scenario_name` to `"Column_Leak_Fault"`.

    ```julia
    # run.jl (modified section)
    function execute_simulation()
        # ...
        material_name = "Zeolite_13X"
        scenario_name = "Column_Leak_Fault" # Changed for this fault
        run_type = "ProcessEvaluation"
        # ...
    end
    ```
5.  **Run the simulation:**
    ```bash
    julia run.jl
    ```
**Expected Output:**
Data will be saved in `simulation_output/Zeolite_13X/Column_Leak_Fault/YYYY-MM-DD_HH-MM-SS/`. You will also see `⚠️ INJECTING FAULT: Column Leak (5% pressure loss)` messages during the simulation run.

**Cleanup:**
*   **Restore `src/PSACycle.jl`:** Immediately restore the original file from your backup.
    ```bash
    cp src/PSACycle.jl.bak src/PSACycle.jl
    rm src/PSACycle.jl.bak # Optional: remove the backup file
    ```
*   **Revert `run.jl`:** Change `scenario_name` back to `"90%_Recovery_(Max_Purity)"`.

This guide provides a comprehensive overview of how to reproduce the fault simulations.
