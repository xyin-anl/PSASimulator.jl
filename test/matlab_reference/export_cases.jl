# Write the 48 demo cases to cases.csv for run_cases.m
# Columns: scenario (1=S4, 2=S5, 3=S10, 4=S11), material row in Params.mat,
# run type (0=ProcessEvaluation, 1=EconomicEvaluation), then x = [P_0, t_ads, alpha, v_0, beta, P_l]
include(joinpath(@__DIR__, "..", "..", "demo", "demo_data.jl"))

scenarios = [(PURITY_RECOVERY_MATERIALS, OPT_VARS_PURITY, 0),
             (PURITY_RECOVERY_MATERIALS, OPT_VARS_RECOVERY, 0),
             (ECONOMIC_MATERIALS, OPT_VARS_PRODUCTIVITY, 1),
             (ECONOMIC_MATERIALS, OPT_VARS_ENERGY, 1)]

open(joinpath(@__DIR__, "cases.csv"), "w") do io
    for (k, (mats, X, run_type)) in enumerate(scenarios), (i, m) in enumerate(mats)
        println(io, join([k, m.index, run_type, repr.(X[i, 1:6])...], ","))
    end
end
