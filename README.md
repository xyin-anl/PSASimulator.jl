<img src="assets/logo.png" alt="PSA Process Diagram" width="300">

#

`PSASimulator.jl` is a Julia-based Pressure Swing Adsorption (PSA) simulation package for modeling gas separation processes, particularly CO₂ capture. This package is a Julia translation of the [original MATLAB-based simulator](https://github.com/PEESEgroup/PSA) at this moment, with additional features and improvements planned in the future.

If you find this package helpful, please cite: [Yin, Xiangyu, and Chrysanthos E. Gounaris. "Computational discovery of Metal–Organic Frameworks for sustainable energy systems: Open challenges." Computers & Chemical Engineering 167 (2022): 108022.](https://www.sciencedirect.com/science/article/pii/S0098135422003568)

[![CI](https://github.com/xyin-anl/PSASimulator.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/xyin-anl/PSASimulator.jl/actions/workflows/CI.yml)

> **Status:** maintenance mode. Bug fixes and dependency updates are welcome; no new features are planned. Related open-source projects: [PSASimulator.py](https://github.com/amvro23/PSASimulator.py) (Python port of the same MATLAB simulator) and the differentiable JAX model of [Glover, Papathanasiou & Pini (2026)](https://arxiv.org/abs/2606.21523).

## Installation

Julia 1.10 or newer is required.

```julia
using Pkg

# Install directly from GitHub
Pkg.add(url="https://github.com/xyin-anl/PSASimulator.jl")

# Or for development (if you want to modify the code):
Pkg.develop(url="https://github.com/xyin-anl/PSASimulator.jl")

# Or if you have cloned the repository locally:
Pkg.develop(path="/path/to/local/PSASimulator.jl")
```

Dependencies: `OrdinaryDiffEqBDF.jl` (stiff ODE solver), `ADTypes.jl`, and the `LinearAlgebra` and `SparseArrays` standard libraries.

## Basic Usage
A high level interface `psacycle` runs the 5-step modified Skarstrom cycle (pressurization, adsorption, heavy reflux, counter-current depressurization, light reflux) until cyclic steady state.

```julia
using PSASimulator

# Define process variables
process_vars = [
    L,      # Bed length [m]
    P_0,    # Feed pressure [Pa]
    n_dot,  # Inlet molar flux [mol/s/m^2]
    t_ads,  # Adsorption time [s]
    alpha,  # Light reflux ratio [-]
    beta,   # Heavy reflux ratio [-]
    P_I,    # Intermediate pressure [Pa]
    P_l     # Light product pressure [Pa]
]

# Material properties: (properties_vector, isotherm_parameters)
material = (material_properties, isotherm_params)

# Run simulation
result = psacycle(process_vars, material;
    N = 10,                           # Number of finite volumes
    run_type = :ProcessEvaluation,    # or :EconomicEvaluation
    y0 = 0.15,                        # CO₂ mole fraction in the feed
    it_disp = true                    # Display iteration progress
)
```

The simulation returns a named tuple with:
- `objectives`: `[-purity, -recovery]` for `:ProcessEvaluation`, or `[-productivity, energy]` for `:EconomicEvaluation`. Productivity is in mol CO₂/kg adsorbent/s and energy in kWh/tonne CO₂, as in the MATLAB code.
- `constraints`: Constraint violations for optimization
- `traj`: Trajectories of the last cycle (`:a` to `:e` for the five steps, `:t1` to `:t5` for their dimensionless times) plus `:purity`, `:recovery` and `:mass_balance`

## Demo

The `demo/` directory evaluates the 16 adsorbents of [Yancy-Caballero et al. (2020)](https://pubs.rsc.org/en/content/articlelanding/2020/me/d0me00060d) at the optimal operating points reported in that paper. It uses its own environment (`demo/Project.toml`) with `DataFrames.jl` and `PrettyTables.jl`:

```bash
cd demo
julia demo_psa_simulator.jl
```

The demo covers four scenarios:
1. Maximized purity with 90% CO₂ recovery constraint
2. Maximized purity with 95% CO₂ recovery constraint
3. Maximized CO₂ productivity
4. Minimized energy consumption [kWh/ton CO₂]

The demo includes 16 different adsorbent materials:
- Metal-Organic Frameworks (MOFs): Co-MOF-74, Cu-BTTri, Mg-MOF-74, MOF-177, etc.
- Zeolites: Zeolite 13X
- Other porous materials: ZIF-8, SIFSIX series

## Validation and Tests

The test suite (`julia --project -e 'using Pkg; Pkg.test()'`) runs all 48 demo cases and compares them with the original MATLAB code run under GNU Octave 10.3 (`test/matlab_reference/`). Purity and recovery agree within 1e-3 (within 2.2e-4 for 47 of the 48 cases; the exception is the slowly converging SIFSIX-3-Ni case), and productivity and energy agree within 0.05% relative. The remaining differences come from the different stiff integrators (`ode15s` vs `QNDF`) and the cyclic steady state tolerance. A demo case takes 1 to 3 s in Julia and 2 to 30 min with the MATLAB code in Octave.

## Mathematical Model

`PSASimulator.jl` solves a system of Partial Differential Equations (PDEs) that describe the dynamics of gas separation in a packed bed column. The model uses dimensionless variables with spatial discretization via finite volumes and WENO schemes.

### Governing Equations

**Component Mass Balance (CO₂ mole fraction)**
```math
\begin{aligned}
\frac{\partial y}{\partial \tau} &= \frac{1}{\mathrm{Pe}}\left(\frac{\partial^2 y}{\partial z^2} + \frac{\partial y}{\partial z}\,\frac{\partial P}{\partial z}\,\frac{1}{P} - \frac{\partial y}{\partial z}\,\frac{\partial T}{\partial z}\,\frac{1}{T}\right)\\
&\quad -\frac{T}{P}\,\frac{\partial}{\partial z}\left(\frac{y P v}{T}\right)
+ \phi\,\frac{T}{P}\Bigl[(y-1)\frac{\partial x_1}{\partial \tau} + y\frac{\partial x_2}{\partial \tau}\Bigr]
\end{aligned}
```

**Total Mass Balance**
```math
\frac{\partial P}{\partial \tau} = -T\,\frac{\partial}{\partial z}\left(\frac{P v}{T}\right) - \phi T\left(\frac{\partial x_1}{\partial \tau} + \frac{\partial x_2}{\partial \tau}\right) + \frac{P}{T}\,\frac{\partial T}{\partial \tau}.
```

**Energy Balance**
```math
\begin{aligned}
\frac{\partial T}{\partial \tau} &= \frac{K_z}{v_0 L}\,\frac{1}{\zeta}\,\frac{\partial^2 T}{\partial z^2}\\
&\quad - \frac{\varepsilon C_{pg} P_0}{R T_0}\,\frac{1}{\zeta}\left[\frac{\partial (P v)}{\partial z} - T\,\frac{\partial}{\partial z}\left(\frac{P v}{T}\right)\right]\\
&\quad + \frac{(1-\varepsilon) q_{s0}}{T_0}\,\frac{1}{\zeta}\Bigl[(-\Delta U_1 + R T_0 T)\frac{\partial x_1}{\partial \tau} + (-\Delta U_2 + R T_0 T)\frac{\partial x_2}{\partial \tau}\Bigr]
\end{aligned}
```
where $\zeta = (1-\varepsilon)(\rho_s C_{ps} + q_{s0} C_{pa}) + \varepsilon\rho_g C_{pg}$ is the effective heat‐capacity term.

**Adsorption Kinetics (Linear Driving Force)**
```math
\frac{\partial x_i}{\partial \tau} = k_i\left(\frac{q_i^{\ast}}{q_{s0}} - x_i\right), \qquad i = 1,2.
```

**Momentum Balance (Ergun Equation)**
```math
-\frac{\partial P}{\partial z}\,\frac{P_0}{L} = \frac{150\,\mu(1-\varepsilon)^2}{4 r_p^{\,2} \varepsilon^{\,2}}\,v + \frac{1.75(1-\varepsilon)\,\rho_g\,\mathrm{MW}}{2 r_p \varepsilon}\,v^2.
```

### Equilibrium Model

The equilibrium loadings $q_1^{\ast}$ and $q_2^{\ast}$ follow a dual–site Langmuir isotherm:
```math
q_i = q_{s,b,i}\,\frac{B_i C_i}{1 + B_1 C_1 + B_2 C_2} + q_{s,d,i}\,\frac{D_i C_i}{1 + D_1 C_1 + D_2 C_2}.
```
with temperature–dependent adsorption constants
```math
B_i = b_i\,e^{-\Delta U_{b,i}/(R T)}, \qquad D_i = d_i\,e^{-\Delta U_{d,i}/(R T)}.
```

### Dimensionless Variables
- Pressure: $P = P_{\text{actual}}/P_0$
- Temperature: $T = T_{\text{actual}}/T_0$
- Velocity: $v = v_{\text{actual}}/v_0$
- Adsorbed amount: $x_i = q_i/q_{s0}$
- Spatial coordinate: $z = z_{\text{actual}}/L$
- Time: $\tau = t v_0 / L$

### Key Parameters
- $\mathrm{Pe} = v_0 L / D_l$ with $D_l = 0.7 D_m + v_0 r_p$ (axial dispersion)
- $\phi = R T_0 q_{s0} (1-\varepsilon) / (\varepsilon P_0)$ (capacity ratio)
- $k_i = k_{i,\text{LDF}} L / v_0$ (dimensionless LDF coefficient)

### Process Steps

The simulator models six process steps, each with specific boundary conditions:

- **Adsorption** - Feed gas at high pressure, CO₂ adsorption
- **Heavy Reflux** - Heavy product recycle for CO₂ recovery enhancement
- **Counter-current Depressurization** - Pressure reduction from feed end
- **Light Reflux** - Light product purge for purity improvement
- **Co-current Depressurization** - Pressure reduction from product end
- **Co-current Pressurization** - Re-pressurization with light product

### Numerical Solution

The PDEs are solved using:
- **Spatial discretization**: Finite volume method with WENO (Weighted Essentially Non-Oscillatory) schemes
- **Temporal integration**: `QNDF` stiff solver from `OrdinaryDiffEqBDF.jl` with sparse finite-difference Jacobians

The simulator iterates through cycles until the state variables at the beginning and end of a cycle converge, indicating cyclic steady state has been reached.