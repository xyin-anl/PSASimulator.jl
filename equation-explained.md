Here is a comprehensive breakdown of the core equations used in the Pressure Swing Adsorption (PSA) simulation. Each equation models a fundamental physical process occurring within the adsorbent bed.

---

### Equation 1: Component Mass Balance (CO₂ mole fraction)

$$
\frac{\partial y}{\partial \tau} = \frac{1}{\mathrm{Pe}}\left(\frac{\partial^2 y}{\partial z^2} + \frac{\partial y}{\partial z}\,\frac{\partial P}{\partial z}\,\frac{1}{P} - \frac{\partial y}{\partial z}\,\frac{\partial T}{\partial z}\,\frac{1}{T}\right) -\frac{T}{P}\,\frac{\partial}{\partial z}\left(\frac{y P v}{T}\right) + \phi\,\frac{T}{P}\Bigl[(y-1)\frac{\partial x_1}{\partial \tau} + y\frac{\partial x_2}{\partial \tau}\Bigr]
$$

* **High-Level Purpose**: This equation tracks how the mole fraction of a specific component (like CO₂) in the gas phase, denoted by **y**, changes over time ($ \tau $) at any point along the length of the column ($ z $).

* **Left-Hand Side**:
    * $ \frac{\partial y}{\partial \tau} $: The **accumulation term**. It represents the rate of change of the CO₂ mole fraction in the gas.

* **Right-Hand Side**:
    1.  **Dispersion Term**: The first term, scaled by the inverse of the **Péclet number** ($ \mathrm{Pe} $), models axial dispersion. This is the "mixing" or "spreading" effect caused by diffusion and varied flow paths, which smears sharp concentration fronts.
    2.  **Advection Term**: The second term, $ -\frac{T}{P}\,\frac{\partial}{\partial z}\left(\frac{y P v}{T}\right) $, represents **advection** (or convection). This is the "flow" term, describing the bulk movement of CO₂ along the column carried by the gas velocity ($ v $).
    3.  **Mass Transfer Term**: The final term is the heart of the separation. It models the rate at which CO₂ is removed from the gas phase and transferred to the solid adsorbent material, driven by the adsorption rates $ \frac{\partial x_1}{\partial \tau} $ and $ \frac{\partial x_2}{\partial \tau} $.

* **Summary**: Rate of CO₂ Change = (Spreading due to Dispersion) - (Movement due to Flow) + (Removal due to Adsorption).

---

### Equation 2: Total Mass Balance

$$
\frac{\partial P}{\partial \tau} = -T\,\frac{\partial}{\partial z}\left(\frac{P v}{T}\right) - \phi T\left(\frac{\partial x_1}{\partial \tau} + \frac{\partial x_2}{\partial \tau}\right) + \frac{P}{T}\,\frac{\partial T}{\partial \tau}
$$

* **High-Level Purpose**: This equation describes how the total gas **pressure (P)** changes over time ($ \tau $) at any point ($ z $) in the bed.

* **Left-Hand Side**:
    * $ \frac{\partial P}{\partial \tau} $: The rate of change of total pressure. A positive value indicates pressurization.

* **Right-Hand Side**:
    1.  **Flow Gradient Term**: $ -T\,\frac{\partial}{\partial z}\left(\frac{P v}{T}\right) $ links pressure change to the gradient in molar gas flow. If more gas flows out of a section than flows in, the pressure in that section drops.
    2.  **Mass Transfer Term**: $ - \phi T\left(\frac{\partial x_1}{\partial \tau} + \frac{\partial x_2}{\partial \tau}\right) $ accounts for pressure change due to molecules leaving the gas phase (adsorption causes pressure to drop) or entering it (desorption causes pressure to rise).
    3.  **Temperature Effect Term**: $ + \frac{P}{T}\,\frac{\partial T}{\partial \tau} $ is a correction based on the Ideal Gas Law. It accounts for pressure changes caused solely by temperature fluctuations.

* **Summary**: Rate of Pressure Change = (Change due to Flow Differences) - (Change due to Adsorption/Desorption) + (Change due to Temperature Swings).

---

### Equation 3: Energy Balance

$$
\frac{\partial T}{\partial \tau} = \frac{K_z}{v_0 L}\,\frac{1}{\zeta}\,\frac{\partial^2 T}{\partial z^2} - \frac{\varepsilon C_{pg} P_0}{R T_0}\,\frac{1}{\zeta}\left[\frac{\partial (P v)}{\partial z} - T\,\frac{\partial }{\partial z}\left(\frac{P v}{T}\right)\right] + \frac{(1-\varepsilon) q_{s0}}{T_0}\,\frac{1}{\zeta}\Bigl[(-\Delta U_1 + R T_0 T)\frac{\partial x_1}{\partial \tau} + (-\Delta U_2 + R T_0 T)\frac{\partial x_2}{\partial \tau}\Bigr]
$$

* **High-Level Purpose**: This equation tracks how the **temperature (T)** changes over time ($ \tau $) by balancing all heat flows and heat generation/consumption. The entire right side is scaled by $ 1/\zeta $, where $ \zeta $ is the effective heat capacity of the bed.

* **Left-Hand Side**:
    * $ \frac{\partial T}{\partial \tau} $: The rate of temperature change.

* **Right-Hand Side**:
    1.  **Heat Conduction Term**: The first term, involving the thermal conductivity ($ K_z $), models the "heat spreading" that smooths out temperature gradients along the bed.
    2.  **Convective Heat Transfer Term**: The second term models the heat carried by the bulk gas flow. This is the "heat flow" term that describes how moving gas transports thermal energy.
    3.  **Heat of Adsorption Term**: The final term is the primary heat source/sink. Adsorption is exothermic and releases heat (driven by $ \Delta U $, the internal energy of adsorption), causing the bed to heat up. Desorption is endothermic and consumes heat, causing the bed to cool down.

* **Summary**: Rate of Temperature Change = (Heat Spreading by Conduction) + (Heat Carried by Gas Flow) + (Heat Generated/Consumed by Adsorption/Desorption).

---

### Equation 4: Adsorption Kinetics (Linear Driving Force)

$$
\frac{\partial x_i}{\partial \tau} = k_i\left(\frac{q_i^{\ast}}{q_{s0}} - x_i\right), \qquad i = 1,2.
$$

* **High-Level Purpose**: This equation describes **how fast** gas molecules are adsorbed onto or desorbed from the solid material.

* **Left-Hand Side**:
    * $ \frac{\partial x_i}{\partial \tau} $: The rate of change of the amount of component `i` adsorbed on the solid.

* **Right-Hand Side**:
    * This side defines the rate based on a simple principle: the rate is proportional to how far the system is from equilibrium.
    * $ k_i $: The **mass transfer coefficient**, a constant representing how quickly molecules can move between the gas and solid phases.
    * $ (q_i^{\ast}/q_{s0} - x_i) $: The **driving force**. This is the difference between the maximum possible adsorbed amount at current conditions ($ q_i^{\ast} $, the equilibrium loading) and the actual current adsorbed amount ($ x_i $). When this difference is large, adsorption or desorption happens quickly.

---

### Equation 5: Momentum Balance (Ergun Equation)

$$
-\frac{\partial P}{\partial z}\,\frac{P_0}{L} = \frac{150\,\mu(1-\varepsilon)^2}{4 r_p^{\,2} \varepsilon^{\,2}}\,v + \frac{1.75(1-\varepsilon)\,\rho_g\,\mathrm{MW}}{2 r_p \varepsilon}\,v^2
$$

* **High-Level Purpose**: This equation calculates the **pressure drop** caused by friction as the gas flows through the packed bed of particles. It relates the gas velocity ($ v $) to the pressure gradient ($ -\frac{\partial P}{\partial z} $).

* **Left-Hand Side**:
    * $ -\frac{\partial P}{\partial z} $: The pressure gradient along the column.

* **Right-Hand Side**: This is the **Ergun equation**, which combines two effects:
    1.  **Viscous Term**: The first term, linear in velocity ($ v $), dominates at low flow rates. It represents pressure loss due to viscous drag (friction) of the gas against the particle surfaces.
    2.  **Inertial Term**: The second term, proportional to velocity squared ($ v^2 $), dominates at high flow rates. It represents pressure loss from the kinetic energy consumed as the gas constantly changes direction around particles.

* **Summary**: Total Pressure Drop = (Pressure Drop from Viscous Friction) + (Pressure Drop from Inertial Effects).

---

### Equations 6 & 7: Equilibrium Model (Dual-Site Langmuir Isotherm)

$$
q_i^* = q_{s,b,i}\,\frac{B_i C_i}{1 + B_1 C_1 + B_2 C_2} + q_{s,d,i}\,\frac{D_i C_i}{1 + D_1 C_1 + D_2 C_2}
$$
$$
B_i = b_i\,e^{-\Delta U_{b,i}/(R T)}, \qquad D_i = d_i\,e^{-\Delta U_{d,i}/(R T)}
$$

* **High-Level Purpose**: These equations together define the **equilibrium loading** ($ q_i^* $), which is the maximum amount of gas that can be adsorbed at a specific pressure, temperature, and composition. This $ q_i^* $ value is the "target" used in the kinetics equation (Eq. 4).

* **Equation 6 (Isotherm)**: This is a **Dual-Site Langmuir isotherm**. It models the surface as having two distinct types of adsorption sites ('b' and 'd').
    * Each term represents the amount adsorbed on one type of site.
    * The model accounts for **competition**, where both gas components ($ C_1 $ and $ C_2 $) compete for the same limited number of sites.
    * The **affinity constants** ($ B_i $ and $ D_i $) represent how strongly each gas sticks to each site type.

* **Equation 7 (Temperature Dependence)**: This equation shows how the affinity constants ($ B_i, D_i $) change with temperature ($ T $).
    * Based on the heat of adsorption ($ \Delta U $), it mathematically describes the fact that adsorption is stronger at lower temperatures. As **T increases**, the affinity constants **decrease**, meaning the material's capacity to hold gas is reduced.
