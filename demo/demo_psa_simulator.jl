# ===================================================================
# PSA SIMULATOR DEMO
# ===================================================================
# 
# This demo tests the PSA Simulator against literature data from:
# Yancy-Caballero et al., "Process-level modelling and optimization 
# to evaluate metal–organic frameworks for post-combustion capture 
# of CO2." Molecular Systems Design & Engineering 5.7 (2020): 1205-1218.
# https://pubs.rsc.org/en/content/articlelanding/2020/me/d0me00060d
#
# Test scenarios:
# - Table S4: Process optimization (max purity, 90% CO2 recovery)
# - Table S5: Process optimization (max purity, 95% CO2 recovery)
# - Table S10: Economic optimization (max productivity)
# - Table S11: Economic optimization (min energy requirement)
# ===================================================================

println("\n" * "="^60)
println(" PSA SIMULATOR DEMO ")
println("="^60 * "\n")

# ===================================================================
# SETUP AND IMPORTS
# ===================================================================

using Pkg
# Use the demo environment, with PSASimulator taken from this checkout
Pkg.activate(@__DIR__)
Pkg.develop(path=dirname(@__DIR__); io=devnull)
Pkg.instantiate()

using PSASimulator
using DataFrames
using PrettyTables
using Statistics

println("✓ All modules loaded successfully")

# ===================================================================
# DATA DEFINITIONS
# ===================================================================

# Include material property data
include("demo_data.jl")

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

"""
    run_psa_simulation(opt_vars, material_data, run_type, N)

Run a single PSA simulation with the given optimization variables and material.

# Arguments
- `opt_vars`: Vector of optimization variables [P_0, t_ads, alpha, v_0, beta, P_l]
- `material_data`: Tuple of (material_properties, isotherm_parameters)
- `run_type`: Either "ProcessEvaluation" or "EconomicEvaluation"
- `N`: Number of discretization points

# Returns
- Tuple of (objectives, constraints)
"""
function run_psa_simulation(opt_vars, material_data, run_type, N=10)
    # Map design variables to process variables (matching MATLAB)
    process_vars = [
        1.0,                                        # L [m] - bed length
        opt_vars[1],                                # P_0 [Pa] - feed pressure  
        opt_vars[1] * opt_vars[4] / 8.314 / 313.15, # n_dot_0 [mol/s] - feed molar flow rate
        opt_vars[2],                                # t_ads [s] - adsorption time
        opt_vars[3],                                # alpha [-] - light reflux ratio
        opt_vars[5],                                # beta [-] - heavy reflux ratio
        1.0e4,                                      # P_I [Pa] - intermediate pressure
        opt_vars[6]                                 # P_l [Pa] - light product pressure
    ]

    try
        result = PSASimulator.psacycle(process_vars, material_data,
            N=N,
            run_type=Symbol(run_type),
            it_disp=false)

        return result.objectives, result.constraints
    catch e
        @warn "Simulation failed: $e"
        return [NaN, NaN], [NaN, NaN, NaN]
    end
end

"""
    parse_results(objectives, constraints, run_type)

Parse simulation results based on the run type.

# Returns
Named tuple with purity, recovery, productivity, and energy values
"""
function parse_results(objectives, constraints, run_type)
    if run_type == "ProcessEvaluation"
        return (
            purity=-objectives[1],
            recovery=-objectives[2],
            productivity=NaN,
            energy=NaN,
            constraints=constraints
        )
    else  # EconomicEvaluation
        return (
            purity=NaN,
            recovery=NaN,
            productivity=-objectives[1],
            energy=objectives[2],
            constraints=constraints
        )
    end
end

"""
    run_scenario_tests(scenario_name, materials_list, opt_vars_matrix, 
                      isotherm_params, sim_params, run_type, N)

Run tests for a complete scenario.
"""
function run_scenario_tests(scenario_name, materials_list, opt_vars_matrix,
    isotherm_params, sim_params, run_type, N=10)

    results = DataFrame(
        Material=String[],
        Purity=Float64[],
        Recovery=Float64[],
        Productivity=Float64[],
        Energy=Float64[],
        Constraint1=Float64[],
        Constraint2=Float64[],
        Constraint3=Float64[]
    )

    n_materials = length(materials_list)
    n_test = min(n_materials, size(opt_vars_matrix, 1))

    println("\nTesting $n_test materials for $scenario_name scenario...")
    println("-"^50)

    for i in 1:n_test
        material_idx = materials_list[i].index
        material_name = materials_list[i].name

        # Extract parameters for this material
        material_iso_params = isotherm_params[material_idx, :]
        material_properties = sim_params[material_idx, :]
        material_data = (material_properties, material_iso_params)

        # Get optimization variables
        opt_vars = opt_vars_matrix[i, :]

        # Run simulation
        objectives, constraints = run_psa_simulation(opt_vars, material_data, run_type, N)

        # Parse results
        res = parse_results(objectives, constraints, run_type)

        # Add to results
        push!(results, (
            material_name,
            res.purity,
            res.recovery,
            res.productivity,
            res.energy,
            res.constraints[1],
            res.constraints[2],
            res.constraints[3]
        ))

        # Display progress
        print("  $(lpad(i, 2)). $(rpad(material_name, 20))")
        if run_type == "ProcessEvaluation"
            println("Purity: $(round(res.purity, digits=4)), " *
                    "Recovery: $(round(res.recovery, digits=4))")
        else
            println("Productivity: $(round(res.productivity, digits=4)), " *
                    "Energy: $(round(res.energy, digits=2))")
        end
    end

    return results
end

"""
    display_scenario_summary(scenario_name, results_df, run_type)

Display formatted summary for a scenario.
"""
function display_scenario_summary(scenario_name, results_df, run_type)
    println("\n📊 $scenario_name Results:")
    println("="^50)

    if run_type == "ProcessEvaluation"
        summary_df = select(results_df, :Material, :Purity, :Recovery)
        pretty_table(summary_df,
            header=["Material", "Purity", "Recovery"],
            formatters=ft_printf("%.4f", [2, 3]),
            alignment=[:l, :r, :r],
            crop=:none)
    else
        summary_df = select(results_df, :Material, :Productivity, :Energy)
        pretty_table(summary_df,
            header=["Material", "Productivity [mol/kg/s]", "Energy [kWh/ton CO₂]"],
            formatters=(ft_printf("%.4f", 2), ft_printf("%.2f", 3)),
            alignment=[:l, :r, :r],
            crop=:none)
    end
end

# ===================================================================
# MAIN EXECUTION
# ===================================================================

function main()
    # Configuration
    N = 10  # Number of discretization points

    # Define test scenarios
    scenarios = [
        (name="90% Recovery (Max Purity)",
            materials=PURITY_RECOVERY_MATERIALS,
            opt_vars=OPT_VARS_PURITY,
            run_type="ProcessEvaluation"), (name="95% Recovery (Max Purity)",
            materials=PURITY_RECOVERY_MATERIALS,
            opt_vars=OPT_VARS_RECOVERY,
            run_type="ProcessEvaluation"), (name="Economic (Max Productivity)",
            materials=ECONOMIC_MATERIALS,
            opt_vars=OPT_VARS_PRODUCTIVITY,
            run_type="EconomicEvaluation"), (name="Economic (Min Energy)",
            materials=ECONOMIC_MATERIALS,
            opt_vars=OPT_VARS_ENERGY,
            run_type="EconomicEvaluation")
    ]

    # Store all results
    all_results = DataFrame()

    println("\n🚀 RUNNING ALL TEST SCENARIOS")
    println("="^60)

    # Run each scenario
    for (i, scenario) in enumerate(scenarios)
        println("\n\n[$i/$(length(scenarios))] $(scenario.name)")
        println("="^60)

        results = run_scenario_tests(
            scenario.name,
            scenario.materials,
            scenario.opt_vars,
            ISOTHERM_PARAMETERS,
            SIMULATION_PARAMETERS,
            scenario.run_type,
            N
        )

        # Add scenario column and append to all results
        results[!, :Scenario] .= scenario.name
        append!(all_results, results)

        # Display scenario summary
        display_scenario_summary(scenario.name, results, scenario.run_type)
    end

    # Display overall summary
    display_overall_summary(all_results)

    return all_results
end

"""
    display_overall_summary(all_results)

Display comprehensive summary statistics for all scenarios.
"""
function display_overall_summary(all_results)
    println("\n\n" * "="^60)
    println("📈 OVERALL SUMMARY STATISTICS")
    println("="^60)

    # Performance metrics summary
    println("\n1. Performance Metrics Summary:")
    println("-"^50)

    perf_stats = DataFrame(
        Metric=String[],
        Value=Float64[],
        Material=String[]
    )

    # Calculate statistics for each metric
    for (metric, col) in [("Purity", :Purity), ("Recovery", :Recovery),
        ("Productivity", :Productivity), ("Energy", :Energy)]
        valid_values = filter(!isnan, all_results[!, col])
        if !isempty(valid_values)
            if metric == "Energy"
                best_val = minimum(valid_values)
                best_idx = findfirst(x -> x == best_val, all_results[!, col])
                metric_label = "Min $metric"
            else
                best_val = maximum(valid_values)
                best_idx = findfirst(x -> x == best_val, all_results[!, col])
                metric_label = "Max $metric"
            end

            push!(perf_stats, (
                metric_label,
                best_val,
                all_results[best_idx, :Material]
            ))

            push!(perf_stats, (
                "Mean $metric",
                mean(valid_values),
                "-"
            ))
        end
    end

    pretty_table(perf_stats,
        header=["Metric", "Value", "Best Material"],
        formatters=ft_printf("%.4f", 2),
        alignment=[:l, :r, :l],
        crop=:none)

    # Constraint violations summary
    println("\n2. Constraint Violations by Scenario:")
    println("-"^50)

    constraint_summary = DataFrame(
        Scenario=String[],
        Tests=Int[],
        C1_Violations=Int[],
        C2_Violations=Int[],
        C3_Violations=Int[]
    )

    for scenario in unique(all_results.Scenario)
        scenario_df = filter(row -> row.Scenario == scenario, all_results)
        push!(constraint_summary, (
            scenario,
            nrow(scenario_df),
            sum(scenario_df.Constraint1 .> 0),
            sum(scenario_df.Constraint2 .> 0),
            sum(scenario_df.Constraint3 .> 0)
        ))
    end

    pretty_table(constraint_summary,
        header=["Scenario", "Tests", "C1", "C2", "C3"],
        title="Constraint Violations",
        alignment=[:l, :r, :r, :r, :r],
        crop=:none)

    # Top performers summary
    println("\n3. Top Performers by Metric:")
    println("-"^50)

    # Best purity
    purity_df = filter(row -> !isnan(row.Purity), all_results)
    if !isempty(purity_df)
        top_purity = sort(purity_df, :Purity, rev=true)[1:min(3, nrow(purity_df)), :]
        println("\n🏆 Top 3 - Purity:")
        pretty_table(select(top_purity, :Material, :Purity, :Scenario),
            header=["Material", "Purity", "Scenario"],
            formatters=ft_printf("%.4f", 2),
            show_row_number=true,
            crop=:none)
    end

    # Best productivity
    prod_df = filter(row -> !isnan(row.Productivity), all_results)
    if !isempty(prod_df)
        top_prod = sort(prod_df, :Productivity, rev=true)[1:min(3, nrow(prod_df)), :]
        println("\n🏆 Top 3 - Productivity:")
        pretty_table(select(top_prod, :Material, :Productivity, :Scenario),
            header=["Material", "Productivity", "Scenario"],
            formatters=ft_printf("%.4f", 2),
            show_row_number=true,
            crop=:none)
    end

    println("\n" * "="^60)
    println("✅ DEMO COMPLETE")
    println("="^60 * "\n")
end

# Run the demo
if abspath(PROGRAM_FILE) == @__FILE__
    results = main()
end