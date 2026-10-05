module FuncHeavyRefluxModule
# =============================================================================
#  FuncHeavyReflux – Complete Julia translation of MATLAB FuncHeavyReflux.m
#  
#  This is a faithful port of the MATLAB heavy reflux step RHS function,
#  including all mathematical terms: spatial derivatives, velocity calculations,
#  energy balance, component mass balance, and boundary conditions.
#  
#  Key differences from other steps:
#  - Heavy end open with reflux gas feed (y_HR, T_HR, ndot_HR)
#  - Light end pressure maintained at 1.0 (atmospheric)
#  - Standard upwind WENO (no adaptive scheme)
#  
#  The function signature matches DifferentialEquations.jl expectations:
#      dx = FuncHeavyReflux(t, x, Params, IsothermPar)
# =============================================================================

using LinearAlgebra

using ...PSAUtils

export FuncHeavyReflux



function FuncHeavyReflux(t, state_vars, params, isotherm_params)
    """
    Calculate derivatives for heavy reflux step - faithful MATLAB port.

    Args:
        t: Time (not used, required by ODE solver)
        state_vars: State vector [P, y, q1, q2, T] for N+2 volumes
        params: Process parameters from PSAInput
        isotherm_params: Isotherm parameters
        
    Returns:
        derivatives: Time derivatives of state variables
    """

    # Extract process parameters (matching MATLAB exactly)
    N = Int(params[1])
    deltaU_1 = params[2]
    deltaU_2 = params[3]
    ro_s = params[4]
    T_0 = params[5]
    epsilon = params[6]
    r_p = params[7]
    mu = params[8]
    R = params[9]
    v_0 = params[10]
    q_s0 = params[11]
    C_pg = params[12]
    C_pa = params[13]
    C_ps = params[14]
    D_m = params[15]
    K_z = params[16]
    P_0 = params[17]
    L = params[18]
    MW_CO2 = params[19]
    MW_N2 = params[20]
    k_1_LDF = params[21]
    k_2_LDF = params[22]
    P_inlet = params[26]
    y_HR = params[33]    # Heavy reflux composition
    T_HR = params[34]    # Heavy reflux temperature  
    ndot_HR = params[35] # Heavy reflux molar flow rate

    # Initialize state variables
    P = zeros(N + 2)
    y = zeros(N + 2)
    x1 = zeros(N + 2)
    x2 = zeros(N + 2)
    T = zeros(N + 2)

    P[1:N+2] = state_vars[1:N+2]
    y[1:N+2] = max.(state_vars[N+3:2*N+4], 0.0)
    x1[1:N+2] = max.(state_vars[2*N+5:3*N+6], 0.0)
    x2[1:N+2] = state_vars[3*N+7:4*N+8]
    T[1:N+2] = state_vars[4*N+9:5*N+10]

    # Initialize derivatives
    derivatives = zeros(5 * N + 10)
    dPdt = zeros(N + 2)
    dydt = zeros(N + 2)
    dx1dt = zeros(N + 2)
    dx2dt = zeros(N + 2)
    dTdt = zeros(N + 2)

    # Initialize spatial derivatives
    dpdz = zeros(N + 2)
    dpdzh = zeros(N + 1)
    dydz = zeros(N + 2)
    d2ydz2 = zeros(N + 2)
    dTdz = zeros(N + 2)
    d2Tdz2 = zeros(N + 2)

    # Calculate parameters
    dz = 1.0 / N
    D_l = 0.7 * D_m + v_0 * r_p
    Pe = v_0 * L / D_l
    phi = R * T_0 * q_s0 * (1 - epsilon) / epsilon / P_0
    ro_g = P[1:N+2] .* P_0 ./ R ./ T[1:N+2] ./ T_0

    # Boundary Conditions for Heavy Reflux
    y[1] = y_HR
    T[1] = T_HR / T_0

    # Heavy end pressure/velocity handling
    if params[end] == 1
        # Constant pressure inlet
        P[1] = P_inlet
    elseif params[end] == 0
        # Constant velocity inlet - solve for pressure
        MW = MW_N2 + (MW_CO2 - MW_N2) * y[1]

        a_1 = 150 * mu * (1 - epsilon)^2 * dz * L / 2 / 4 / r_p^2 / epsilon^3 / T[1] / T_0 / R
        a_2_1 = 1.75 * (1 - epsilon) / 2 / r_p / epsilon^3 * dz * L / 2
        a_2 = a_2_1 / R / T[1] / T_0 * ndot_HR * MW

        a = a_1 + a_2
        b = P[2] / T[1] * P_0 / R / T_0
        c = -ndot_HR

        vh_inlet = (-b + sqrt(b^2 - 4 * a * c)) / 2 / a / v_0

        a_p = a_1 * T[1] * T_0 * R
        b_p = a_2_1 * MW / R / T[1] / T_0

        P[1] = ((a_p * vh_inlet * v_0 + P[2] * P_0) / (1 - b_p * (vh_inlet * v_0)^2)) / P_0
    else
        error("Please specify whether inlet velocity or pressure is constant")
    end

    # Light end boundary conditions
    y[N+2] = y[N+1]
    T[N+2] = T[N+1]
    if P[N+1] >= 1
        P[N+2] = 1.0  # Atmospheric pressure
    else
        P[N+2] = P[N+1]
    end

    # Velocity Calculations using centralized Ergun equation solver
    vh, Ph, Th, yh, dpdzh = compute_velocity(P, T, y, params; weno_scheme="upwind")

    # Spatial Derivative Calculations using wall values from velocity calculation
    dpdz[2:N+1] = (Ph[2:N+1] - Ph[1:N]) / dz
    dydz[2:N+1] = (yh[2:N+1] - yh[1:N]) / dz
    dTdz[2:N+1] = (Th[2:N+1] - Th[1:N]) / dz

    # Second derivatives
    d2ydz2[3:N] = (y[4:N+1] + y[2:N-1] - 2 * y[3:N]) / dz / dz
    d2ydz2[2] = (y[3] - y[2]) / dz / dz
    d2ydz2[N+1] = (y[N] - y[N+1]) / dz / dz

    d2Tdz2[3:N] = (T[4:N+1] + T[2:N-1] - 2 * T[3:N]) / dz / dz
    d2Tdz2[2] = 4 * (Th[2] + T[1] - 2 * T[2]) / dz / dz
    d2Tdz2[N+1] = 4 * (Th[N] + T[N+2] - 2 * T[N+1]) / dz / dz

    # Temporal Derivatives

    # 1) Adsorbed Mass Balance (LDF) - identical to other steps
    q = Isotherm(y, P .* P_0, T .* T_0, isotherm_params)
    q_1 = q[:, 1] .* ro_s
    q_2 = q[:, 2] .* ro_s

    k_1 = k_1_LDF * L / v_0
    k_2 = k_2_LDF * L / v_0

    dx1dt[2:N+1] = k_1 * (q_1[2:N+1] / q_s0 - x1[2:N+1])
    dx2dt[2:N+1] = k_2 * (q_2[2:N+1] / q_s0 - x2[2:N+1])

    # 2) Column Energy Balance - identical to other steps
    sink_term = (1 - epsilon) * (ro_s * C_ps + q_s0 * C_pa) .+ epsilon .* ro_g[2:N+1] .* C_pg

    # Conduction
    transfer_term = K_z / v_0 / L
    dTdt1 = transfer_term .* d2Tdz2[2:N+1] ./ sink_term

    # Advection
    PvT = Ph[1:N+1] .* vh[1:N+1] ./ Th[1:N+1]
    Pv = Ph[1:N+1] .* vh[1:N+1]
    dTdt2 = -epsilon .* C_pg .* P_0 ./ R ./ T_0 .*
            ((Pv[2:N+1] - Pv[1:N]) - T[2:N+1] .* (PvT[2:N+1] - PvT[1:N])) ./ dz ./ sink_term

    # Heat of adsorption
    generation_term_1 = (1 - epsilon) .* q_s0 .* (-(deltaU_1 .- R .* T[2:N+1] .* T_0)) ./ T_0
    generation_term_2 = (1 - epsilon) .* q_s0 .* (-(deltaU_2 .- R .* T[2:N+1] .* T_0)) ./ T_0
    dTdt3 = (generation_term_1 .* dx1dt[2:N+1] .+ generation_term_2 .* dx2dt[2:N+1]) ./ sink_term

    dTdt[2:N+1] = dTdt1 + dTdt2 + dTdt3

    # 3) Total mass balance - identical to other steps
    dPdt1 = -T[2:N+1] .* (PvT[2:N+1] - PvT[1:N]) ./ dz
    dPdt2 = -phi * T[2:N+1] .* (dx1dt[2:N+1] + dx2dt[2:N+1])
    dPdt3 = P[2:N+1] .* dTdt[2:N+1] ./ T[2:N+1]
    dPdt[2:N+1] = dPdt1 + dPdt2 + dPdt3

    # 4) Component Mass Balance - identical to other steps
    # Diffusion
    dydt1 = (1 / Pe) .* (d2ydz2[2:N+1] + (dydz[2:N+1] .* dpdz[2:N+1] ./ P[2:N+1]) -
                         (dydz[2:N+1] .* dTdz[2:N+1] ./ T[2:N+1]))

    # Advection
    ypvt = yh[1:N+1] .* Ph[1:N+1] .* vh[1:N+1] ./ Th[1:N+1]
    dydt2 = -(T[2:N+1] ./ P[2:N+1]) .* ((ypvt[2:N+1] - ypvt[1:N]) -
                                        y[2:N+1] .* (PvT[2:N+1] - PvT[1:N])) ./ dz

    # Adsorption/desorption
    dydt3 = (phi * T[2:N+1] ./ P[2:N+1]) .* ((y[2:N+1] .- 1) .* dx1dt[2:N+1] .+
                                             y[2:N+1] .* dx2dt[2:N+1])

    dydt[2:N+1] = dydt1 + dydt2 + dydt3

    # Boundary Derivatives - KEY DIFFERENCES for Heavy Reflux
    dPdt[1] = 0.0                     # Fixed pressure at inlet
    dPdt[N+2] = 0.0                   # Fixed pressure at outlet
    dydt[1] = 0.0                     # Fixed composition at inlet
    dydt[N+2] = dydt[N+1]             # No composition gradient at outlet
    dx1dt[1] = 0.0
    dx2dt[1] = 0.0
    dx1dt[N+2] = 0.0
    dx2dt[N+2] = 0.0
    dTdt[1] = 0.0                     # Fixed temperature at inlet
    dTdt[N+2] = dTdt[N+1]             # No temperature gradient at outlet

    # Export derivatives
    derivatives[1:N+2] = dPdt[1:N+2]
    derivatives[N+3:2*N+4] = dydt[1:N+2]
    derivatives[2*N+5:3*N+6] = dx1dt[1:N+2]
    derivatives[3*N+7:4*N+8] = dx2dt[1:N+2]
    derivatives[4*N+9:5*N+10] = dTdt[1:N+2]

    return derivatives
end



end # module FuncHeavyRefluxModule 