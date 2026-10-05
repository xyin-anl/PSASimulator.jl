module PSAInput
# =========================================================================
#  process_input_parameters – faithful port of MATLAB "ProcessInputParameters"
# =========================================================================
export process_input_parameters

"""
    process_input_parameters(process_vars, material, N;
                             feed_gas = "Constant Velocity", y0 = 0.15)

Translate the MATLAB routine *ProcessInputParameters*.

# Arguments
- `process_vars` : length‑8 Vector{<:Real}
      [L, P₀, ṅ₀, t_ads, α, β, P_I, P_l]
- `material`     : Tuple (material_property::Vector, isotherm_par::Vector)
      exactly the same two arrays MATLAB expected as `material{1}` and `material{2}`
- `N`            : integer – number of finite volumes
- `feed_gas`     : `"Constant Pressure"` or `"Constant Velocity"`  
      (defaults to the latter, matching the MATLAB file)
- `y0`           : CO₂ mole fraction in the feed (default 0.15, as in MATLAB)

# Returns
(Params          ::Vector{Float64},   # length 39
 IsothermParams  ::Vector{Float64},   # length 13
 Times           ::Vector{Float64},   # 6‑element vector of step times  [s]
 EconomicParams  ::Vector{Float64})   # length 7

All numeric results are `Float64`, matching MATLAB's default.
"""
function process_input_parameters(process_vars::AbstractVector,
      material::Tuple,
      N::Integer;
      feed_gas::AbstractString="Constant Velocity",
      y0::Real=0.15)

      @assert length(process_vars) == 8 "process_vars must have length 8"
      L, P₀, ṅ₀, t_ads, α, β, P_I, P_l = Float64.(process_vars)

      material_property, isoPar = material
      @assert length(material_property) ≥ 3 "material_property vector should contain (ρ_s, ΔU₁, ΔU₂, …)"
      @assert length(isoPar) ≥ 13 "isotherm_par needs at least 13 elements (matches MATLAB)"

      # ─────────────────────────────────────────────────────────────────────
      #  Constants (identical to MATLAB)
      # ─────────────────────────────────────────────────────────────────────
      R = 8.314             # J mol⁻¹ K⁻¹
      T₀ = 313.15            # K
      y₀ = Float64(y0)
      0 < y₀ < 1 || throw(ArgumentError("y0 must be in (0, 1), got $y0"))
      Ctot₀ = P₀ / R / T₀
      v₀ = ṅ₀ / Ctot₀
      μ = 1.72e-5           # Pa·s
      ε = 0.37
      D_m = 1.2995e-5         # m²/s
      K_z = 0.09              # W m⁻¹ K⁻¹
      C_pg = 30.7              # J mol⁻¹ K⁻¹
      C_pa = 30.7
      MW_CO2 = 0.04402           # kg mol⁻¹
      MW_N2 = 0.02802
      r_p = 1e-3              # m
      C_ps = 1070.0            # J kg⁻¹ K⁻¹
      q_s = 5.84              # mol kg⁻¹
      ρ_s = material_property[1]
      q_s0 = q_s * ρ_s
      k_CO2_LDF = 0.1631
      k_N2_LDF = 0.2044
      ΔU = (material_property[2], material_property[3])

      # Operating‑step default durations (s)
      t_pres = 20.0
      t_CnCdepres = 30.0
      t_CoCdepres = 70.0
      t_LR = t_ads
      t_HR = t_LR
      τ = 0.5               # pressure‑change speed factor
      P_inlet = 1.02              # dimensionless (×P₀ later)

      # ─── Isotherm parameters (dual‑site) ────────────────────────────────
      q_s_b = [isoPar[1], isoPar[7]]
      q_s_d = [isoPar[2], isoPar[8]]
      b = [isoPar[3], isoPar[9]]
      d = [isoPar[4], isoPar[10]]
      ΔU_b = [isoPar[5], isoPar[11]]
      ΔU_d = [isoPar[6], isoPar[12]]
      isotherm_extra = isoPar[13]       # the 13th entry is passed through

      IsothermParams = vcat(q_s_b, q_s_d, b, d, ΔU_b, ΔU_d, isotherm_extra)

      # ─── Assemble Params (length 39, positions exactly match MATLAB) ────
      Params = zeros(Float64, 39)
      Params[1] = N
      Params[2] = ΔU[1]
      Params[3] = ΔU[2]
      Params[4] = ρ_s
      Params[5] = T₀
      Params[6] = ε
      Params[7] = r_p
      Params[8] = μ
      Params[9] = R
      Params[10] = v₀
      Params[11] = q_s0
      Params[12] = C_pg
      Params[13] = C_pa
      Params[14] = C_ps
      Params[15] = D_m
      Params[16] = K_z
      Params[17] = P₀
      Params[18] = L
      Params[19] = MW_CO2
      Params[20] = MW_N2
      Params[21] = k_CO2_LDF
      Params[22] = k_N2_LDF
      Params[23] = y₀
      Params[24] = τ
      Params[25] = P_l
      Params[26] = P_inlet
      Params[27] = 1.0         # y_LP placeholder
      Params[28] = 1.0         # T_LP placeholder
      Params[29] = 1.0         # ṅ_LP placeholder
      Params[30] = α
      Params[31] = β
      Params[32] = P_I
      Params[33] = y₀          # y_HR initialise
      Params[34] = T₀          # T_HR initialise
      Params[35] = ṅ₀ * β
      Params[36] = 0.01        # y_LR guess
      Params[37] = T₀          # T_LR guess
      Params[38] = ṅ₀
      Params[39] = feed_gas == "Constant Pressure" ? 1.0 :
                   feed_gas == "Constant Velocity" ? 0.0 :
                   error("feed_gas must be \"Constant Pressure\" or \"Constant Velocity\"")

      # ─── Times vector (s) – order: press, ads, CnCDepres, LR, CoCDepres, HR
      Times = [t_pres, t_ads, t_CnCdepres, t_LR, t_CoCdepres, t_HR]

      # ─── Economic parameters (unchanged) ────────────────────────────────
      desired_flow = 100.0
      electricity_cost = 0.07
      hour_to_year_conversion = 8000.0
      life_span_equipment = 20.0
      life_span_adsorbent = 5.0
      CEPCI = 536.4

      cycle_time = t_pres + t_ads + t_HR + t_CnCdepres + t_LR
      EconomicParams = [
            desired_flow,
            electricity_cost,
            cycle_time,
            hour_to_year_conversion,
            life_span_equipment,
            life_span_adsorbent,
            CEPCI
      ]

      return (Params=Params,
            IsothermParams=IsothermParams,
            Times=Times,
            EconomicParams=EconomicParams)
end

end # module PSAInput
