# Simulating Faults and Anomalies in PSASimulator.jl

While `PSASimulator.jl` does not have a built-in feature for explicitly modeling faults, its structure is well-suited for simulating anomalies. You can achieve this by programmatically intervening in the simulation process to mimic the effects of various faults.

## General Approach

The most effective method is to compare a baseline (normal) simulation with a faulty one:

1.  **Run a normal simulation** and save the resulting trajectory (`traj` object).
2.  **Introduce a fault** by modifying the simulation logic or inputs.
3.  **Run the faulty simulation** and save its trajectory.
4.  **Analyze the differences** between the normal and faulty trajectories to develop and test your fault detection algorithms.

---

## How to Simulate Different Faults

Here’s a breakdown of how you could approach this for different types of faults, primarily by modifying the main simulation loop in `src/PSACycle.jl`.

### 1. Actuator Faults (e.g., Stuck or Leaky Valves)

These are faults in the mechanical components that control the process. You can simulate them by altering the parameters between the simulation steps inside the `psacycle` function.

*   **Example**: Heavy Reflux valve stuck closed.
*   **Location to Modify**: In `src/PSACycle.jl`, after the "Adsorption Step" (line 250) and before the "Heavy Reflux Step" (line 283).
*   **Logic**: The flow rate for the heavy reflux step (`ndot_HR`) is calculated based on the output of the previous step. To simulate a stuck-closed valve, you can manually set this flow rate to zero.
*   **Implementation**:

    ```julia
    // Inside PSACycle.jl, around line 330
    ...
    totalFront, CO2Front, TFront = stream_composition_calculator(...)
    push!(c_fin, [c[end, :]; CO2Front; totalFront; CO2End; totalEnd])

    # --- FAULT INJECTION: HEAVY REFLUX VALVE STUCK ---
    y_HR = 0.0 
    T_HR = T_0 # Reset to initial temp
    ndot_HR = 0.0 // Set flow rate to zero
    # --- END FAULT INJECTION ---

    # Update parameters for Heavy Reflux step
    # y_HR = CO2Front / totalFront   // Original line
    # T_HR = TFront                  // Original line
    # ndot_HR = totalFront * beta / t_HR // Original line
    Params[33] = y_HR
    Params[34] = T_HR
    Params[35] = ndot_HR
    ...
    ```

### 2. Process Faults (e.g., Leaks, Channeling, Feed Gas Changes)

These are faults related to the process conditions or the integrity of the column itself.

*   **Example**: A leak in the column during the high-pressure adsorption step.
*   **Location to Modify**: In `src/PSACycle.jl`, within the main loop, after the ODE solver for the "Adsorption Step" (line 240).
*   **Logic**: A leak would cause a pressure drop. You can simulate this by manually reducing the pressure in the state vector `b` (the output of the adsorption step) before it's used as the input for the next step.
*   **Implementation**:

    ```julia
    // Inside PSACycle.jl, around line 241
    ...
    sol2 = solve(prob2, QNDF(autodiff=false); ...)
    t2 = sol2.t
    b = reduce(hcat, sol2.u)'

    # --- FAULT INJECTION: SIMULATE A LEAK ---
    # Reduce the pressure in all nodes by 5%
    leakage_factor = 0.95
    b[:, 1:N+2] .*= leakage_factor 
    # --- END FAULT INJECTION ---

    # Correct the output (clean up results from simulation)
    ...
    ```

### 3. Adsorbent Degradation

This fault represents the aging or poisoning of the adsorbent material, leading to reduced performance. This is one of the easiest faults to simulate.

*   **Example**: The adsorbent has lost 10% of its CO₂ adsorption capacity.
*   **Location to Modify**: In your script that calls `psacycle` (like `demo/demo_psa_simulator.jl`).
*   **Logic**: Adsorption capacity is defined by the isotherm parameters. You can modify these parameters before passing them to the simulator.
*   **Implementation**:

    ```julia
    // In your calling script
    ...
    # Get material data
    material_properties, isotherm_params = get_material_data("Zeolite_13X")

    # --- FAULT INJECTION: ADSORBENT DEGRADATION ---
    # Reduce the saturation capacity (q_s_b and q_s_d) for CO2 by 10%
    degradation_factor = 0.9
    isotherm_params[1] *= degradation_factor // q_s_b for CO2
    isotherm_params[2] *= degradation_factor // q_s_d for CO2
    # --- END FAULT INJECTION ---

    material_data = (material_properties, isotherm_params)

    # Run the simulation with the faulty material
    result = PSASimulator.psacycle(process_vars, material_data)
    ...
    ```
