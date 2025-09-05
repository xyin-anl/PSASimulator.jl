module PSACycleDriver
# =============================================================================
#  PSACycleDriver – Complete Julia translation of PSACycle.m 
#  
#  This module implements a faithful port of the MATLAB PSACycle.m with:
#  - 700-iteration convergence algorithm with sophisticated CSS checking
#  - Complete stream composition calculator
#  - Economic evaluation functions 
#  - Velocity correction and cleanup functions
#  - Process evaluation with purity, recovery, and mass balance
#  
#  Matches MATLAB functionality exactly including boundary condition handling,
#  convergence criteria, and post-processing.
# =============================================================================

using DifferentialEquations
using LinearAlgebra
using Statistics
using OrdinaryDiffEq
using SparseArrays

# Get the directory of this file (src/)
const SRC_DIR = dirname(@__FILE__)

# Import required modules
include(joinpath(SRC_DIR, "PSAInput.jl"))
using .PSAInput

include(joinpath(SRC_DIR, "PSAUtils.jl"))
using .PSAUtils

include(joinpath(SRC_DIR, "StepModels.jl"))
using .StepModels

export psacycle

"""
    psacycle(vars, material; x0=nothing, N=10, it_disp=false, run_type=:ProcessEvaluation)

Main PSA cycle simulation function - faithful port of MATLAB PSACycle.m

# Arguments
- `vars::Vector`: Process variables [L, P_0, n_dot, t_ads, alpha, beta, P_I, P_l]
- `material::Tuple`: Material properties (adsorption_data, isotherm_par)
- `x0::Union{Vector,Nothing}=nothing`: Initial conditions (auto-generated if nothing)
- `N::Int=10`: Number of finite volumes
- `it_disp::Bool=false`: Display iteration information
- `run_type::Symbol=:ProcessEvaluation`: :ProcessEvaluation or :EconomicEvaluation

# Returns
- Named tuple with objectives, constraints, and trajectory data
"""
function psacycle(vars::Vector{<:Real}, material::Tuple;
    x0::Union{Vector{<:Real},Nothing}=nothing,
    N::Int=10,
    it_disp::Bool=false,
    run_type::Symbol=:ProcessEvaluation,
    max_iters::Int=700)

    # Initialize objectives and constraints
    objectives = [0.0, 0.0]
    constraints = [0.0, 0.0, 0.0]

    # Process input parameters (matching MATLAB exactly)
    ip = process_input_parameters(vars, material, N)
    Params = ip.Params
    IsothermParams = ip.IsothermParams
    Times = ip.Times
    EconomicParams = ip.EconomicParams

    # Extract key parameters (matching MATLAB variable names)
    N = Int(Params[1])
    ro_s = Params[4]
    T_0 = Params[5]
    epsilon = Params[6]
    r_p = Params[7]
    mu = Params[8]
    R = Params[9]
    v_0 = Params[10]
    q_s0 = Params[11] / ro_s
    P_0 = Params[17]
    L = Params[18]
    MW_CO2 = Params[19]
    MW_N2 = Params[20]
    y_0 = Params[23]
    ndot_0 = vars[3]
    P_l = Params[25]
    P_inlet = Params[26]
    alpha = Params[30]
    beta = Params[31]
    y_HR = Params[33]
    T_HR = Params[34]
    ndot_HR = Params[35]

    # Step times (dimensional and dimensionless)
    t_CoCPres = Times[1]
    t_ads = Times[2]
    t_HR = Times[6]
    t_CnCDepres = Times[3]
    t_LR = Times[4]

    tau_CoCPres = t_CoCPres * v_0 / L
    tau_ads = t_ads * v_0 / L
    tau_HR = t_HR * v_0 / L
    tau_CnCDepres = t_CnCDepres * v_0 / L
    tau_LR = t_LR * v_0 / L

    # Initialize column (matching MATLAB exactly)
    if x0 === nothing
        q = Isotherm([y_0], [P_l], [298.15], IsothermParams)[1, :]
        x0 = zeros(Float64, 5 * N + 10)
        x0[1:N+2] .= P_l / P_0
        x0[N+3] = y_0
        x0[N+4:2*N+4] .= y_0
        x0[2*N+5:3*N+6] .= q[1] / q_s0
        x0[3*N+7:4*N+8] .= q[2] / q_s0
        x0[4*N+9] = 1.0
        x0[4*N+10:5*N+10] .= 298.15 / T_0
    else
        x0 = Float64.(x0)
    end

    # ODE solver options (matching MATLAB tolerances)
    ode_opts = Dict(
        :reltol => 1e-6,  # Match MATLAB tolerances
        :abstol => 1e-6   # Match MATLAB tolerances
    )

    # Storage for trajectory data (matching MATLAB)
    a_storage = []
    b_storage = []
    c_storage = []
    d_storage = []
    e_storage = []
    t1_storage = []
    t2_storage = []
    t3_storage = []
    t4_storage = []
    t5_storage = []

    # Storage for final conditions analysis
    a_fin = []
    b_fin = []
    c_fin = []
    d_fin = []
    e_fin = []

    # Skip simulation if constraints already violated
    if constraints[1] == 0
        # Initialize statesIC for CSS check
        statesIC = nothing

        # Main cyclic loop
        for i in 1:max_iters
            if it_disp
                println("Iteration for CSS condition number: $i")
            end

            # Store initial conditions for analysis (get statesIC at the start of each iteration)
            if i == 1
                statesIC = x0[vcat(2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9)]
            end

            #= ============================================================
            1. Co-current Pressurization Step
            ============================================================ =#
            rhs1 = step_rhs(:CoCPressurization, Params, IsothermParams)
            jac1 = JacPressurization(N)  # Get Jacobian pattern

            f1 = ODEFunction(rhs1; jac_prototype=jac1)
            prob1 = ODEProblem(f1, x0, (0.0, tau_CoCPres))
            sol1 = solve(prob1, QNDF(autodiff=false);  # Use efficient stiff solver without autodiff
                reltol=ode_opts[:reltol],
                abstol=ode_opts[:abstol],
                maxiters=1e5)

            t1 = sol1.t
            a = reduce(hcat, sol1.u)'  # Convert to MATLAB-style matrix

            # Cleanup results (matching MATLAB exactly)
            idx = findall(a[:, 1] .< a[:, 2])
            a[idx, 1] = a[idx, 2]                        # P_1 = P_2
            a[idx, N+3] = a[idx, N+4]                    # y_1 = y_2
            a[idx, 4*N+9] = a[idx, 4*N+10]              # T_1 = T_2
            a[:, 2*N+5] = a[:, 2*N+6]                   # x1_1 = x1_2
            a[:, 3*N+7] = a[:, 3*N+8]                   # x2_1 = x2_2
            a[:, 3*N+6] = a[:, 3*N+5]                   # x1_N+2 = x1_N+1
            a[:, 4*N+8] = a[:, 4*N+7]                   # x2_N+2 = x2_N+1
            a[:, N+3:2*N+4] = max.(min.(a[:, N+3:2*N+4], 1), 0)  # 0 <= y <= 1

            # Stream composition calculation
            totalFront, CO2Front, _ = stream_composition_calculator(t1 * L / v_0, a, Params, "HPEnd")
            totalEnd, CO2End, _ = stream_composition_calculator(t1 * L / v_0, a, Params, "LPEnd")
            push!(a_fin, [a[end, :]; CO2Front; totalFront; CO2End; totalEnd])

            # Prepare initial conditions for Adsorption
            x10 = copy(a[end, :])
            x10[1] = P_inlet
            x10[N+2] = 1.0
            x10[N+3] = y_0
            x10[2*N+4] = x10[2*N+3]
            x10[4*N+9] = 1.0
            x10[5*N+10] = x10[5*N+9]

            #= ============================================================
            2. Adsorption Step
            ============================================================ =#
            rhs2 = step_rhs(:Adsorption, Params, IsothermParams)
            jac2 = JacAdsorption(N)  # Get Jacobian pattern

            f2 = ODEFunction(rhs2; jac_prototype=jac2)
            prob2 = ODEProblem(f2, x10, (0.0, tau_ads))
            sol2 = solve(prob2, QNDF(autodiff=false);  # Use efficient stiff solver without autodiff
                reltol=ode_opts[:reltol],
                abstol=ode_opts[:abstol],
                maxiters=1e5)

            t2 = sol2.t
            b = reduce(hcat, sol2.u)'

            # Correct the output (clean up results from simulation)
            idx = findall(b[:, N+1] .< 1)
            b[idx, N+2] = b[idx, N+1]
            b[:, 2*N+5] = b[:, 2*N+6]
            b[:, 3*N+7] = b[:, 3*N+8]
            b[:, 3*N+6] = b[:, 3*N+5]
            b[:, 4*N+8] = b[:, 4*N+7]
            b[:, N+3:2*N+4] = clamp.(b[:, N+3:2*N+4], 0, 1)

            if Params[end] == 0
                # Velocity cleanup to match MATLAB behavior
                b = velocity_cleanup(collect(b), Params)
            end

            # Stream composition and parameter updates for Light Reflux
            totalFront, CO2Front, _ = stream_composition_calculator(t2 * L / v_0, b, Params, "HPEnd")
            totalEnd, CO2End, TEnd = stream_composition_calculator(t2 * L / v_0, b, Params, "LPEnd")
            push!(b_fin, [b[end, :]; CO2Front; totalFront; CO2End; totalEnd])

            # Update parameters for Light Reflux step
            y_LR = CO2End / totalEnd
            T_LR = TEnd
            ndot_LR = totalEnd / t_ads
            Params[27] = y_LR
            Params[28] = T_LR
            Params[29] = ndot_LR

            # Prepare initial conditions for Heavy Reflux
            x20 = copy(b[end, :])
            x20[1] = x20[2]
            x20[N+2] = 1.0
            x20[N+3] = y_HR
            x20[2*N+4] = x20[2*N+3]
            x20[4*N+9] = T_HR / T_0
            x20[5*N+10] = x20[5*N+9]

            #= ============================================================
            3. Heavy Reflux Step
            ============================================================ =#
            rhs3 = step_rhs(:HeavyReflux, Params, IsothermParams)
            jac3 = JacAdsorption(N)  # Heavy Reflux uses same pattern as Adsorption

            f3 = ODEFunction(rhs3; jac_prototype=jac3)
            prob3 = ODEProblem(f3, x20, (0.0, tau_HR))
            sol3 = solve(prob3, QNDF(autodiff=false);  # Use efficient stiff solver without autodiff
                reltol=ode_opts[:reltol],
                abstol=ode_opts[:abstol],
                maxiters=1e5)

            t3 = sol3.t
            c = reduce(hcat, sol3.u)'

            # Cleanup results
            idx = findall(c[:, N+1] .< 1)
            c[idx, N+2] = c[idx, N+1]
            c[:, 2*N+5] = c[:, 2*N+6]
            c[:, 3*N+7] = c[:, 3*N+8]
            c[:, 3*N+6] = c[:, 3*N+5]
            c[:, 4*N+8] = c[:, 4*N+7]
            c[:, N+3:2*N+4] = max.(min.(c[:, N+3:2*N+4], 1), 0)

            # Velocity correction if enabled
            if Params[end] == 0
                # Velocity correction to match MATLAB behavior
                c = velocity_correction(collect(c), ndot_HR, "HPEnd", Params)
            end

            totalFront, CO2Front, _ = stream_composition_calculator(t3 * L / v_0, c, Params, "HPEnd")
            totalEnd, CO2End, _ = stream_composition_calculator(t3 * L / v_0, c, Params, "LPEnd")
            push!(c_fin, [c[end, :]; CO2Front; totalFront; CO2End; totalEnd])

            # Prepare initial conditions for CnCDepressurization
            x30 = copy(c[end, :])
            x30[1] = x30[2]
            x30[N+2] = x30[N+1]
            x30[N+3] = x30[N+4]
            x30[2*N+4] = x30[2*N+3]
            x30[4*N+9] = x30[4*N+10]
            x30[5*N+10] = x30[5*N+9]

            #= ============================================================
            4. Counter-current Depressurization Step
            ============================================================ =#
            rhs4 = step_rhs(:CnCDepressurization, Params, IsothermParams)
            jac4 = Jac_CnCDepressurization(N)  # Get Jacobian pattern

            f4 = ODEFunction(rhs4; jac_prototype=jac4)
            prob4 = ODEProblem(f4, x30, (0.0, tau_CnCDepres))
            sol4 = solve(prob4, QNDF(autodiff=false);  # Use efficient stiff solver without autodiff
                reltol=ode_opts[:reltol],
                abstol=ode_opts[:abstol],
                maxiters=1e5)

            t4 = sol4.t
            d = reduce(hcat, sol4.u)'

            # Cleanup results
            idx = findall(d[:, 2] .< d[:, 1])
            d[idx, 1] = d[idx, 2]
            d[:, 2*N+5] = d[:, 2*N+6]
            d[:, 3*N+7] = d[:, 3*N+8]
            d[:, 3*N+6] = d[:, 3*N+5]
            d[:, 4*N+8] = d[:, 4*N+7]
            d[:, N+3:2*N+4] = max.(min.(d[:, N+3:2*N+4], 1), 0)

            totalFront, CO2Front, _ = stream_composition_calculator(t4 * L / v_0, d, Params, "HPEnd")
            totalEnd, CO2End, _ = stream_composition_calculator(t4 * L / v_0, d, Params, "LPEnd")
            push!(d_fin, [d[end, :]; CO2Front; totalFront; CO2End; totalEnd])

            # Prepare initial conditions for Light Reflux
            x40 = copy(d[end, :])
            x40[1] = P_l / P_0
            x40[N+3] = x40[N+4]
            x40[2*N+4] = y_LR
            x40[4*N+9] = x40[4*N+10]
            x40[5*N+10] = T_LR / T_0

            #= ============================================================
            5. Light Reflux Step
            ============================================================ =#
            rhs5 = step_rhs(:LightReflux, Params, IsothermParams)
            jac5 = Jac_LightReflux(N)  # Get Jacobian pattern

            f5 = ODEFunction(rhs5; jac_prototype=jac5)
            prob5 = ODEProblem(f5, x40, (0.0, tau_LR))
            sol5 = solve(prob5, QNDF(autodiff=false);  # Use efficient stiff solver without autodiff
                reltol=ode_opts[:reltol],
                abstol=ode_opts[:abstol],
                maxiters=1e5)

            t5 = sol5.t
            e = reduce(hcat, sol5.u)'

            # Correct the output (clean up results from simulation)
            idx = findall(e[:, 2] .< e[:, 1])
            e[idx, 1] = e[idx, 2]
            e[:, 2*N+5] = e[:, 2*N+6]
            e[:, 3*N+7] = e[:, 3*N+8]
            e[:, 3*N+6] = e[:, 3*N+5]
            e[:, 4*N+8] = e[:, 4*N+7]
            e[:, N+3:2*N+4] = clamp.(e[:, N+3:2*N+4], 0, 1)

            # Velocity correction to match MATLAB behavior
            e = velocity_correction(collect(e), ndot_LR * alpha, "LPEnd", Params)

            totalFront, CO2Front, TFront = stream_composition_calculator(t5 * L / v_0, e, Params, "HPEnd")
            totalEnd, CO2End, _ = stream_composition_calculator(t5 * L / v_0, e, Params, "LPEnd")
            push!(e_fin, [e[end, :]; CO2Front; totalFront; CO2End; totalEnd])

            # Update parameters for Heavy Reflux step
            y_HR = CO2Front / totalFront
            T_HR = TFront
            ndot_HR = totalFront * beta / t_HR
            Params[33] = y_HR
            Params[34] = T_HR
            Params[35] = ndot_HR

            # Prepare initial conditions for next cycle
            x0 = copy(e[end, :])
            x0[1] = x0[2]
            x0[N+2] = x0[N+1]
            x0[N+3] = y_0
            x0[2*N+4] = x0[2*N+3]
            x0[4*N+9] = 1.0
            x0[5*N+10] = x0[5*N+9]

            # Store trajectory data
            push!(a_storage, a)
            push!(b_storage, b)
            push!(c_storage, c)
            push!(d_storage, d)
            push!(e_storage, e)
            push!(t1_storage, t1)
            push!(t2_storage, t2)
            push!(t3_storage, t3)
            push!(t4_storage, t4)
            push!(t5_storage, t5)

            #= ============================================================
            Cyclic Steady State (CSS) Check - Matching MATLAB exactly
            ============================================================ =#
            statesFC = e[end, vcat(2:N+1, N+4:2*N+3, 2*N+6:3*N+5, 3*N+8:4*N+7, 4*N+10:5*N+9)]

            if i > 1  # Need at least 2 iterations to check convergence
                CSS_states = norm(statesIC - statesFC)

                # Mass balance calculation
                purity, recovery, massBalance = process_evaluation(
                    a, b, c, d, e, t1, t2, t3, t4, t5, Params)

                # Check CSS condition (matching MATLAB exactly)
                if CSS_states <= 1e-3 && abs(massBalance - 1) <= 0.005
                    if it_disp
                        println("Converged at iteration $i")
                        println("CSS_states = $CSS_states")
                        println("Mass balance = $massBalance")
                    end
                    break
                end

                if it_disp && i % 10 == 0
                    println("Iteration $i: CSS_states = $CSS_states, MB = $massBalance")
                end
            end

            # Update statesIC for next iteration
            statesIC = statesFC
        end
    end

    #= ================================================================
    Process and Economic Evaluation
    ================================================================ =#

    # Use last cycle for evaluation
    a, b, c, d, e = a_storage[end], b_storage[end], c_storage[end],
    d_storage[end], e_storage[end]
    t1, t2, t3, t4, t5 = t1_storage[end], t2_storage[end], t3_storage[end],
    t4_storage[end], t5_storage[end]

    purity, recovery, MB = process_evaluation(a, b, c, d, e, t1, t2, t3, t4, t5, Params)

    # Calculate constraints
    con = recovery / MB - 0.9
    if con < 0
        constraints[2] = abs(con)
    end

    if run_type == :ProcessEvaluation
        objectives[1] = -purity
        objectives[2] = -recovery / MB

        con = purity - y_0
        if con < 0
            constraints[3] = 0.0  # Note: MATLAB sets this to 0 in this case
        end
    else  # EconomicEvaluation
        economic_eval = economic_evaluation(
            a, b, c, d, e, t1, t2, t3, t4, t5, Params, EconomicParams, vars)
        objectives[1] = -economic_eval[:productivity]
        objectives[2] = economic_eval[:energy_requirements]

        con = purity - 0.9
        if con < 0
            constraints[3] = abs(con)
        end
    end

    # Return trajectory data (all cycles)
    traj = Dict(
        :a_storage => a_storage, :b_storage => b_storage, :c_storage => c_storage, :d_storage => d_storage, :e_storage => e_storage,
        :t1_storage => t1_storage, :t2_storage => t2_storage, :t3_storage => t3_storage, :t4_storage => t4_storage, :t5_storage => t5_storage,
        :purity => purity, :recovery => recovery, :mass_balance => MB
    )

    return (objectives=objectives, constraints=constraints, traj=traj)
end

# =============================================================================
# Supporting Functions (faithful MATLAB ports)
# =============================================================================

"""
Stream composition calculator - faithful port of MATLAB StreamCompositionCalculator
"""
function stream_composition_calculator(time, state_vars, Params, ProductEnd="HPEnd")
    N = Int(Params[1])
    L = Params[18]
    R = Params[9]
    T_0 = Params[5]
    P_0 = Params[17]
    MW_CO2 = Params[19]
    MW_N2 = Params[20]
    epsilon = Params[6]
    mu = Params[8]
    r_p = Params[7]

    dz = L / N

    ndot_tot = similar(time)
    ndot_CO2 = similar(time)
    Temp_vec = similar(time)

    for (i, τ) in enumerate(time)
        x = state_vars[i, :]

        if ProductEnd == "HPEnd"
            P = x[1:2] .* P_0
            y = x[N+3]
            T = x[4*N+9] * T_0
        else  # LPEnd
            P = x[N+1:N+2] .* P_0
            y = x[2*N+4]
            T = x[5*N+10] * T_0
        end

        # Gas density and concentrations
        # Use P[1] for HPEnd (node 1), P[2] for LPEnd (node N+2)
        P_node = ProductEnd == "HPEnd" ? P[1] : P[2]
        ro_g = (y * MW_CO2 + (1 - y) * MW_N2) * P_node / R / T
        C_tot = P_node / R / T
        C_CO2 = C_tot * y

        # Pressure gradient and velocity calculation
        dPdz = 2 * (P[2] - P[1]) / dz
        viscous_term = 150 * mu * (1 - epsilon)^2 / 4 / r_p^2 / epsilon^3
        kinetic_term = 1.75 * (1 - epsilon) / 2 / r_p / epsilon^3 * ro_g

        if abs(kinetic_term) > 1e-10
            v = -sign(dPdz) * (-viscous_term + sqrt(viscous_term^2 + 4 * kinetic_term * abs(dPdz))) / 2 / kinetic_term
        else
            v = 0.0
        end

        ndot_tot[i] = abs(v * C_tot)
        ndot_CO2[i] = abs(v * C_CO2)
        Temp_vec[i] = T
    end

    n_tot = trapz(time, ndot_tot)
    n_CO2 = trapz(time, ndot_CO2)

    # Average temperature weighted by mass flow
    energy_flux_tot = ndot_tot .* Temp_vec
    energy_tot = trapz(time, energy_flux_tot)
    Temp = energy_tot / n_tot

    return n_tot, n_CO2, Temp
end

"""
Process evaluation - faithful port of MATLAB ProcessEvaluation
"""
function process_evaluation(a, b, c, d, e, t1, t2, t3, t4, t5, Params)
    # Convert dimensionless times to dimensional
    L = Params[18]
    v_0 = Params[10]
    beta = Params[31]

    t1_dim = t1 .* L ./ v_0
    t2_dim = t2 .* L ./ v_0
    t3_dim = t3 .* L ./ v_0
    t4_dim = t4 .* L ./ v_0
    t5_dim = t5 .* L ./ v_0

    # Calculate stream compositions for all steps
    _, n_CO2_CoCPres_HPEnd, _ = stream_composition_calculator(t1_dim, a, Params, "HPEnd")
    _, n_CO2_ads_HPEnd, _ = stream_composition_calculator(t2_dim, b, Params, "HPEnd")
    _, n_CO2_ads_LPEnd, _ = stream_composition_calculator(t2_dim, b, Params, "LPEnd")
    _, n_CO2_HR_LPEnd, _ = stream_composition_calculator(t3_dim, c, Params, "LPEnd")
    _, n_CO2_HR_HPEnd, _ = stream_composition_calculator(t3_dim, c, Params, "HPEnd")
    n_tot_CnC_HPEnd, n_CO2_CnC_HPEnd, _ = stream_composition_calculator(t4_dim, d, Params, "HPEnd")
    _, n_CO2_LR_LPEnd, _ = stream_composition_calculator(t5_dim, e, Params, "LPEnd")
    n_tot_LR_HPEnd, n_CO2_LR_HPEnd, _ = stream_composition_calculator(t5_dim, e, Params, "LPEnd")

    # Calculate purity, recovery, and mass balance
    purity = (n_CO2_CnC_HPEnd + (1 - beta) * n_CO2_LR_HPEnd) /
             (n_tot_CnC_HPEnd + (1 - beta) * n_tot_LR_HPEnd)

    recovery = (n_CO2_CnC_HPEnd + (1 - beta) * n_CO2_LR_HPEnd) /
               (n_CO2_CoCPres_HPEnd + n_CO2_ads_HPEnd)

    mass_balance = (n_CO2_CnC_HPEnd + n_CO2_ads_LPEnd + n_CO2_HR_LPEnd + n_CO2_LR_LPEnd) /
                   (n_CO2_CoCPres_HPEnd + n_CO2_ads_HPEnd + n_CO2_HR_HPEnd + n_CO2_LR_LPEnd)

    return purity, recovery, mass_balance
end

"""
Velocity cleanup function - faithful port of MATLAB velocitycleanup
"""
function velocity_cleanup(x, Params)
    N = Int(Params[1])
    epsilon = Params[6]
    r_p = Params[7]
    mu = Params[8]
    R = Params[9]
    v_0 = Params[10]
    T_0 = Params[5]
    P_0 = Params[17]
    L = Params[18]
    MW_CO2 = Params[19]
    MW_N2 = Params[20]

    x_new = copy(x)
    numb1 = 150 * mu * (1 - epsilon)^2 / 4 / r_p^2 / epsilon^2
    ro_gent = x[:, 2] .* P_0 / R / T_0
    numb2_ent = ro_gent .* (MW_N2 .+ (MW_CO2 - MW_N2) .* x[:, N+3]) .* (1.75 * (1 - epsilon) / 2 / r_p / epsilon)

    x_new[:, 1] = (numb1 * v_0 .+ numb2_ent .* v_0 .* v_0) .* L ./ P_0 ./ 2 ./ N .+ x[:, 2]

    return x_new
end

"""
Velocity correction function - faithful port of MATLAB VelocityCorrection
"""
function velocity_correction(x, n_hr, CorrectionEnd, Params)
    N = Int(Params[1])
    epsilon = Params[6]
    r_p = Params[7]
    mu = Params[8]
    R = Params[9]
    T_0 = Params[5]
    P_0 = Params[17]
    L = Params[18]
    MW_CO2 = Params[19]
    MW_N2 = Params[20]

    x_new = copy(x)
    dz = L / N

    if CorrectionEnd == "HPEnd"
        T = x[:, 4*N+9] .* T_0
        y = x[:, N+3]
        P = x[:, 2] .* P_0
    elseif CorrectionEnd == "LPEnd"
        T = x[:, 5*N+10] .* T_0
        y = x[:, 2*N+4]
        P = x[:, N+1] .* P_0
    else
        error("CorrectionEnd must be 'HPEnd' or 'LPEnd'")
    end

    MW = MW_N2 .+ (MW_CO2 - MW_N2) .* y

    a_1 = 150 * mu * (1 - epsilon)^2 * dz / 2 / 4 / r_p^2 / epsilon^3 ./ R ./ T
    a_2_1 = 1.75 * (1 - epsilon) / 2 / r_p / epsilon / epsilon / epsilon * dz / 2
    a_2 = a_2_1 ./ R ./ T .* n_hr .* MW

    a_a = a_1 .+ a_2
    b_b = P ./ T / R
    c_c = -n_hr

    vh = (-b_b .+ sqrt.(b_b .^ 2 .- 4 .* a_a .* c_c)) ./ 2 ./ a_a

    a_p = a_1 .* T * R
    b_p = a_2_1 .* MW ./ R ./ T

    if CorrectionEnd == "HPEnd"
        x_new[:, 1] = ((a_p .* vh .+ P) ./ (1 .- b_p .* vh .* vh)) ./ P_0
    elseif CorrectionEnd == "LPEnd"
        x_new[:, N+2] = ((a_p .* vh .+ P) ./ (1 .- b_p .* vh .* vh)) ./ P_0
    end

    return x_new
end

"""
Economic evaluation - faithful port of MATLAB economic calculations
"""
function economic_evaluation(a, b, c, d, e, t1, t2, t3, t4, t5, Params, EconomicParams, vars)
    # Extract parameters
    L = Params[18]
    v_0 = Params[10]
    epsilon = Params[6]
    ro_s = Params[4]
    MW_CO2 = Params[19]
    beta = Params[31]

    # Calculate step times from dimensional time arrays
    # t1, t2, etc. are dimensionless time arrays from ODE solver
    # Convert final times to get step durations in seconds
    t_CoCPres = t1[end] * L / v_0
    t_ads = t2[end] * L / v_0
    t_HR = t3[end] * L / v_0
    t_CnCDepres = t4[end] * L / v_0
    t_LR = t5[end] * L / v_0

    # Calculate cycle time (matching MATLAB)
    cycle_time = t_CoCPres + t_ads + t_HR + t_CnCDepres + t_LR

    # Calculate gas fed during cycle
    n_tot_pres, _, _ = stream_composition_calculator(t1 .* L ./ v_0, a, Params, "HPEnd")
    n_tot_ads, _, _ = stream_composition_calculator(t2 .* L ./ v_0, b, Params, "HPEnd")
    gas_fed = n_tot_pres + n_tot_ads

    # Calculate required radius
    desired_flow = EconomicParams[1]
    radius_inner = sqrt((desired_flow * cycle_time / gas_fed) / π)

    # Calculate energy requirements
    E_pres = compression_energy(t1 .* L ./ v_0, a, 1e5, Params, radius_inner)
    E_feed = compression_energy(t2 .* L ./ v_0, b, 1e5, Params, radius_inner)
    E_HR = compression_energy(t3 .* L ./ v_0, c, 1e5, Params, radius_inner)
    E_bldwn = vacuum_energy(t4 .* L ./ v_0, d, 1e5, "HPEnd", Params, radius_inner)
    E_evac = vacuum_energy(t5 .* L ./ v_0, e, 1e5, "HPEnd", Params, radius_inner)

    energy_per_cycle = E_pres + E_feed + E_HR + E_bldwn + E_evac

    # Calculate CO2 recovered
    _, n_CO2_CnCD, _ = stream_composition_calculator(t4 .* L ./ v_0, d, Params, "HPEnd")
    _, n_CO2_LR, _ = stream_composition_calculator(t5 .* L ./ v_0, e, Params, "HPEnd")

    CO2_recovered_cycle = (n_CO2_CnCD + (1 - beta) * n_CO2_LR) * radius_inner^2 * π * MW_CO2 / 1e3
    CO2_recovered_cycle2 = (n_CO2_CnCD + (1 - beta) * n_CO2_LR) * radius_inner^2 * π

    # Calculate productivity and energy requirements
    mass_adsorbent = L * π * radius_inner^2 * (1 - epsilon) * ro_s
    productivity = CO2_recovered_cycle2 / cycle_time / mass_adsorbent
    energy_requirements = energy_per_cycle / CO2_recovered_cycle

    return Dict(:productivity => productivity, :energy_requirements => energy_requirements)
end

"""
Compression energy calculation - faithful port of MATLAB CompressionEnergy
"""
function compression_energy(time::AbstractVector, state_vars::AbstractMatrix, Patm::Real, Params::AbstractVector, r_in::Real)
    # Extract parameters
    N = Int(Params[1])
    L = Params[18]
    P_0 = Params[17]
    epsilon = Params[6]
    mu = Params[8]
    r_p = Params[7]
    R = Params[9]
    T_0 = Params[5]
    MW_CO2 = Params[19]
    MW_N2 = Params[20]

    # Compressor parameters
    adiabatic_index = 1.4
    compressor_efficiency = 0.72
    dz = L / N

    # Dimensionalize variables at the heavy product end
    P = state_vars[:, 1:2] .* P_0
    y = state_vars[:, N+3]
    T = state_vars[:, 4*N+9] .* T_0

    # Gas density
    ro_g = (y .* MW_CO2 .+ (1 .- y) .* MW_N2) .* P[:, 1] ./ R ./ T

    # Pressure gradient
    dPdz = 2 .* (P[:, 2] .- P[:, 1]) ./ dz

    # Superficial velocity using Ergun equation
    viscous_term = 150 * mu * (1 - epsilon)^2 / 4 / r_p^2 / epsilon^3
    kinetic_term = (1.75 * (1 - epsilon) / 2 / r_p / epsilon^3) .* ro_g

    # Safeguard against zero kinetic term
    safe_kinetic_term = kinetic_term .+ 1e-12
    v = -sign.(dPdz) .* (-viscous_term .+ sqrt.(abs.(viscous_term^2 .+ 4 .* safe_kinetic_term .* abs.(dPdz)))) ./ (2 .* safe_kinetic_term)

    # Calculate compression ratio term
    ratio_term = (P[:, 1] ./ Patm) .^ ((adiabatic_index - 1) / adiabatic_index) .- 1
    ratio_term = max.(ratio_term, 0)

    # Calculate integral term for energy
    integral_term = abs.(v .* P[:, 1] .* ratio_term)

    # Calculate energy in Joules
    energy_J = trapz(time, integral_term) * (adiabatic_index / (adiabatic_index - 1)) / compressor_efficiency * π * r_in^2

    # Convert to kWh
    return energy_J / 3.6e6
end

"""
Vacuum energy calculation - faithful port of MATLAB VacuumEnergy
"""
function vacuum_energy(time::AbstractVector, state_vars::AbstractMatrix, Patm::Real, ProductEnd::String, Params::AbstractVector, r_in::Real)
    # Extract parameters
    N = Int(Params[1])
    L = Params[18]
    P_0 = Params[17]
    epsilon = Params[6]
    mu = Params[8]
    r_p = Params[7]
    R = Params[9]
    T_0 = Params[5]
    MW_CO2 = Params[19]
    MW_N2 = Params[20]

    # Vacuum parameters
    adiabatic_index = 1.4
    vacuum_efficiency = 0.72
    dz = L / N

    # Dimensionalize variables based on the column end
    P_out = zeros(size(state_vars, 1))
    P = zeros(size(state_vars, 1), 2)
    y = zeros(size(state_vars, 1))
    T = zeros(size(state_vars, 1))
    ro_g = zeros(size(state_vars, 1))

    if ProductEnd == "HPEnd"
        P = state_vars[:, 1:2] .* P_0
        y = state_vars[:, N+3]
        T = state_vars[:, 4*N+9] .* T_0
        ro_g = (y .* MW_CO2 .+ (1 .- y) .* MW_N2) * P[:, 1] ./ R ./ T
        P_out = P[:, 1]
    elseif ProductEnd == "LPEnd"
        P = state_vars[:, N+1:N+2] .* P_0
        y = state_vars[:, 2*N+4]
        T = state_vars[:, 5*N+10] .* T_0
        ro_g = (y .* MW_CO2 .+ (1 .- y) .* MW_N2) * P[:, 2] ./ R ./ T
        P_out = P[:, 2]
    else
        error("CorrectionEnd must be 'HPEnd' or 'LPEnd'")
    end

    # Pressure gradient
    dPdz = 2 .* (P[:, 2] .- P[:, 1]) ./ dz

    # Superficial velocity using Ergun equation
    viscous_term = 150 * mu * (1 - epsilon)^2 / 4 / r_p^2 / epsilon^3
    kinetic_term = (1.75 * (1 - epsilon) / 2 / r_p / epsilon^3) .* ro_g

    safe_kinetic_term = kinetic_term .+ 1e-12
    v = -sign.(dPdz) .* (-viscous_term .+ sqrt.(abs.(viscous_term^2 .+ 4 .* safe_kinetic_term .* abs.(dPdz)))) ./ (2 .* safe_kinetic_term)

    # Calculate compression ratio term for vacuum
    ratio_term = (Patm ./ P_out) .^ ((adiabatic_index - 1) / adiabatic_index) .- 1
    ratio_term = max.(ratio_term, 0)

    # Calculate integral term for energy
    integral_term = abs.(v .* P_out .* ratio_term)

    # Calculate energy in Joules
    energy_J = trapz(time, integral_term) * (adiabatic_index / (adiabatic_index - 1)) / vacuum_efficiency * π * r_in^2

    # Convert to kWh
    return energy_J / 3.6e6
end

# Utility function already available in PSAUtils
# trapz function should be imported from PSAUtils

# =============================================================================
# Jacobian Pattern Functions (faithful MATLAB ports)
# =============================================================================

"""
Jacobian pattern for pressurization step (CoCPressurization)
"""
function JacPressurization(N::Int)
    # Create sparse matrices directly for efficiency
    # Four band Jacobian scheme for advection terms
    A4 = spdiagm(-2 => ones(N), -1 => ones(N + 1), 0 => ones(N + 2), 1 => ones(N + 1))

    # One band Jacobian scheme for adsorption/desorption term
    A1 = spdiagm(0 => ones(N + 2))
    A1[1, 1] = 0
    A1[N+2, N+2] = 0

    # Zero band Jacobian Term
    A0 = spzeros(N + 2, N + 2)

    # Create Overall Jacobian based on individual segments
    J_pres = [A4 A4 A1 A1 A4;
        A4 A4 A1 A1 A4;
        A1 A1 A1 A0 A1;
        A1 A1 A0 A1 A1;
        A4 A1 A1 A1 A4]

    # Modify based on boundary conditions
    # Pressure Inlet
    J_pres[1, :] .= 0
    J_pres[1, 1] = 1

    # Pressure Outlet
    J_pres[N+2, :] = J_pres[N+1, :]
    J_pres[:, N+2] .= 0

    # Mole Fraction Inlet
    J_pres[N+3, :] .= 0
    J_pres[:, N+3] .= 0

    # Mole Fraction Outlet
    J_pres[2*N+4, :] = J_pres[2*N+3, :]
    J_pres[:, 2*N+4] .= 0

    # Temperature Inlet
    J_pres[4*N+9, :] .= 0
    J_pres[:, 4*N+9] .= 0

    # Temperature Outlet
    J_pres[5*N+10, :] = J_pres[5*N+9, :]
    J_pres[:, 5*N+10] .= 0

    return sparse(J_pres)
end

"""
Jacobian pattern for adsorption step
"""
function JacAdsorption(N::Int)
    # Create sparse matrices directly for efficiency
    # Four band Jacobian scheme for advection terms
    A4 = spdiagm(-2 => ones(N), -1 => ones(N + 1), 0 => ones(N + 2), 1 => ones(N + 1))

    # One band Jacobian scheme for adsorption/desorption term
    A1 = spdiagm(0 => ones(N + 2))
    A1[1, 1] = 0
    A1[N+2, N+2] = 0

    # Zero band Jacobian Term
    A0 = spzeros(N + 2, N + 2)

    # Create Overall Jacobian based on individual segments
    J_ads = [A4 A4 A1 A1 A4;
        A4 A4 A1 A1 A4;
        A1 A1 A1 A0 A1;
        A1 A1 A0 A1 A1;
        A4 A1 A1 A1 A4]

    # Modify based on boundary conditions
    # Pressure Inlet
    J_ads[1, :] .= 0
    J_ads[:, 1] .= 0

    # Pressure Outlet
    J_ads[N+2, :] .= 0
    J_ads[:, N+2] .= 0

    # Mole Fraction Inlet
    J_ads[N+3, :] .= 0
    J_ads[:, N+3] .= 0

    # Mole Fraction Outlet
    J_ads[2*N+4, :] = J_ads[2*N+3, :]
    J_ads[:, 2*N+4] .= 0

    # Temperature Inlet
    J_ads[4*N+9, :] .= 0
    J_ads[:, 4*N+9] .= 0

    # Temperature Outlet
    J_ads[5*N+10, :] = J_ads[5*N+9, :]
    J_ads[:, 5*N+10] .= 0

    return sparse(J_ads)
end

"""
Jacobian pattern for CnC depressurization step
"""
function Jac_CnCDepressurization(N::Int)
    # Create sparse matrices directly for efficiency
    # Four band Jacobian scheme for advection terms
    A4 = spdiagm(-1 => ones(N + 1), 0 => ones(N + 2), 1 => ones(N + 1), 2 => ones(N))
    # Fix the first and last rows
    for j in 1:size(A4, 2)
        A4[1, j] = A4[2, j]
        A4[N+2, j] = A4[N+1, j]
    end

    # One band Jacobian scheme for adsorption/desorption term
    A1 = spdiagm(0 => ones(N + 2))
    A1[1, 1] = 0
    A1[N+2, N+2] = 0
    A1[1, 2] = 1
    A1[N+2, N+1] = 1

    # Zero band Jacobian Term
    A0 = spzeros(N + 2, N + 2)

    # Create Overall Jacobian based on individual segments
    J_CnCdepres = [A4 A4 A1 A1 A4;
        A4 A4 A1 A1 A4;
        A1 A1 A1 A0 A1;
        A1 A1 A0 A1 A1;
        A4 A1 A1 A1 A4]

    # Modify based on boundary conditions
    # Pressure Inlet
    J_CnCdepres[1, :] .= 0
    J_CnCdepres[1, 1] = 1

    # Pressure Outlet
    J_CnCdepres[N+2, :] = J_CnCdepres[N+1, :]

    # Mole Fraction Inlet
    J_CnCdepres[N+3, :] = J_CnCdepres[N+4, :]

    # Mole Fraction Outlet
    J_CnCdepres[2*N+4, :] = J_CnCdepres[2*N+3, :]

    # Molar Loading
    J_CnCdepres[2*N+5, :] .= 0
    J_CnCdepres[3*N+6:3*N+7, :] .= 0
    J_CnCdepres[4*N+8, :] .= 0

    # Temperature Inlet
    J_CnCdepres[4*N+9, :] = J_CnCdepres[4*N+10, :]

    # Temperature Outlet
    J_CnCdepres[5*N+10, :] = J_CnCdepres[5*N+9, :]

    return sparse(J_CnCdepres)
end

"""
Jacobian pattern for light reflux step
"""
function Jac_LightReflux(N::Int)
    # Create sparse matrices directly for efficiency
    # Four band Jacobian scheme for advection terms
    A4 = spdiagm(-1 => ones(N + 1), 0 => ones(N + 2), 1 => ones(N + 1), 2 => ones(N))
    # Fix the first and last rows
    for j in 1:size(A4, 2)
        A4[1, j] = A4[2, j]
        A4[N+2, j] = A4[N+1, j]
    end

    # One band Jacobian scheme for adsorption/desorption term
    A1 = spdiagm(0 => ones(N + 2))
    A1[1, 1] = 0
    A1[N+2, N+2] = 0
    A1[1, 2] = 1
    A1[N+2, N+1] = 1

    # Zero band Jacobian Term
    A0 = spzeros(N + 2, N + 2)

    # Create Overall Jacobian based on individual segments
    J_LR = [A4 A1 A1 A1 A4;
        A4 A4 A1 A1 A4;
        A1 A1 A1 A0 A1;
        A1 A1 A0 A1 A1;
        A4 A1 A1 A1 A4]

    # Modify based on boundary conditions
    # Pressure Inlet
    J_LR[1, :] .= 0
    J_LR[1, 1] = 1

    # Pressure Outlet
    J_LR[N+2, :] = J_LR[N+1, :]

    # Mole Fraction Inlet
    J_LR[N+3, :] = J_LR[N+4, :]

    # Mole Fraction Outlet
    J_LR[2*N+4, :] = J_LR[2*N+3, :]

    # Molar Loading
    J_LR[2*N+5, :] .= 0
    J_LR[3*N+6:3*N+7, :] .= 0
    J_LR[4*N+8, :] .= 0

    # Temperature Inlet
    J_LR[4*N+9, :] = J_LR[4*N+10, :]

    # Temperature Outlet
    J_LR[5*N+10, :] = J_LR[5*N+9, :]

    return sparse(J_LR)
end

end # module PSACycleDriver