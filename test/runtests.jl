using PSASimulator
using Test

include(joinpath(@__DIR__, "..", "demo", "demo_data.jl"))

# Same mapping from design variables to process variables as PSACycleSimulation.m
# x = [P_0, t_ads, alpha, v_0, beta, P_l]
process_vars(x) = [1.0, x[1], x[1] * x[4] / 8.314 / 313.15, x[2], x[3], x[5], 1.0e4, x[6]]
material(idx) = (SIMULATION_PARAMETERS[idx, :], ISOTHERM_PARAMETERS[idx, :])

const SCENARIOS = Dict(
    "S4" => (OPT_VARS_PURITY, PURITY_RECOVERY_MATERIALS, :ProcessEvaluation),
    "S5" => (OPT_VARS_RECOVERY, PURITY_RECOVERY_MATERIALS, :ProcessEvaluation),
    "S10" => (OPT_VARS_PRODUCTIVITY, ECONOMIC_MATERIALS, :EconomicEvaluation),
    "S11" => (OPT_VARS_ENERGY, ECONOMIC_MATERIALS, :EconomicEvaluation),
)

function read_reference(path)
    lines = readlines(path)
    header = split(lines[1], ',')
    return [Dict(zip(header, split(l, ','))) for l in lines[2:end]]
end

@testset "PSASimulator" begin
    @testset "process_input_parameters" begin
        x = OPT_VARS_PURITY[4, :]
        ip = process_input_parameters(process_vars(x), material(4), 10)
        @test length(ip.Params) == 39
        @test length(ip.IsothermParams) == 13
        @test length(ip.Times) == 6
        @test length(ip.EconomicParams) == 7
        @test ip.Params[1] == 10
        @test ip.Params[23] == 0.15

        ip = process_input_parameters(process_vars(x), material(4), 10; y0=0.3)
        @test ip.Params[23] == 0.3
        @test ip.Params[33] == 0.3
        @test_throws ArgumentError process_input_parameters(process_vars(x), material(4), 10; y0=1.2)
    end

    @testset "Isotherm" begin
        iso = process_input_parameters(process_vars(OPT_VARS_PURITY[4, :]), material(4), 10).IsothermParams
        P = [1e4, 5e4, 1e5]
        q = Isotherm(fill(0.15, 3), P, fill(313.15, 3), iso)
        @test all(q .>= 0)
        @test issorted(q[:, 1])  # CO2 loading increases with pressure
    end

    # Regression against the original MATLAB code (run under GNU Octave, see
    # matlab_reference/) at the optimal operating points reported by
    # Yancy-Caballero et al. (2020). Observed differences (integrator and CSS
    # tolerance) are below 2.2e-4 except for S4 SIFSIX-3-Ni (1e-3), and below
    # 0.05% relative for the economic objectives.
    @testset "agreement with MATLAB" begin
        for ref in read_reference(joinpath(@__DIR__, "matlab_reference", "reference.csv"))
            X, mats, run_type = SCENARIOS[ref["scenario"]]
            row = findfirst(m -> m.index == parse(Int, ref["idx"]), mats)
            r = psacycle(process_vars(X[row, :]), material(mats[row].index); N=10, run_type=run_type)
            obj = parse.(Float64, (ref["obj1"], ref["obj2"]))
            @testset "$(ref["scenario"]) $(ref["material"])" begin
                if run_type == :ProcessEvaluation  # [-purity, -recovery]
                    @test r.objectives[1] ≈ obj[1] atol = 2e-3
                    @test r.objectives[2] ≈ obj[2] atol = 2e-3
                else  # [-productivity, energy]
                    @test r.objectives[1] ≈ obj[1] rtol = 5e-3
                    @test r.objectives[2] ≈ obj[2] rtol = 5e-3
                end
                @test abs(r.traj[:mass_balance] - 1) < 0.01
            end
        end
    end

    @testset "feed composition" begin
        x = OPT_VARS_PURITY[4, :]
        default = psacycle(process_vars(x), material(4); N=10)
        explicit = psacycle(process_vars(x), material(4); N=10, y0=0.15)
        richer = psacycle(process_vars(x), material(4); N=10, y0=0.30)
        @test explicit.traj[:purity] ≈ default.traj[:purity] atol = 1e-6
        @test richer.traj[:purity] > default.traj[:purity]
        @test 0.99 < richer.traj[:mass_balance] < 1.01
    end
end
