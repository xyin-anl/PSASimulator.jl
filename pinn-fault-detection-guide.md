# Fault Detection with Physics-Informed Neural Networks (PINNs) in PSASimulator.jl

Physics-Informed Neural Networks (PINNs) offer a powerful method for creating intelligent "digital twins" of physical systems. For a complex, dynamic process like Pressure Swing Adsorption, they can be used to build sophisticated fault detection and diagnosis systems. This guide outlines the high-level strategy and a practical path for implementing such a system.

---

## Section 1: High-Level Strategy - Which Equations to Use?

A crucial strategic question is whether to build one massive model containing all the governing equations or to use a more modular approach.

### The "Holistic" vs. "Decomposed" Approach

#### The "Holistic" Model: A Single, Comprehensive PINN

The idealized approach is to build one large PINN constrained by the entire system of coupled PDEs (Component Mass, Total Mass, Energy, Kinetics, and Momentum). 

*   **Pros**: This model would be incredibly powerful, offering deep diagnostic insights by observing how a fault causes violations across all the physical laws simultaneously.
*   **Cons**: This is a "monster model." In practice, training a single network with many complex, competing physics-based loss terms is extremely difficult. Issues with balancing the influence of each loss term and achieving stable convergence are significant research challenges.

#### The "Decomposed" Hybrid Model: A Practical and Robust Alternative

A more pragmatic and intelligent strategy is to build a hybrid, two-stage system.

**Stage 1: General Anomaly Detection**
This stage answers the simple question: *"Is the system operating normally?"*

*   **Model**: An **Autoencoder** trained exclusively on data from a wide range of normal operating conditions.
*   **Function**: It acts as a fast, robust, and data-driven first line of defense. If the model fails to reconstruct the incoming sensor data accurately, it flags an anomaly and triggers Stage 2.

**Stage 2: Specialist Fault Diagnosis Ensemble**
This stage answers the follow-up question: *"What is wrong?"*

*   **Model**: An **ensemble of small, specialized PINNs**. Each PINN is an expert on only one part of the physics.
*   **Example Ensemble Members**:
    *   **"Leak Detector" PINN**: Constrained only by the **Total Mass Balance** equation.
    *   **"Degradation Detector" PINN**: Constrained only by the **Component Mass Balance** and **Kinetics** equations.
    *   **"Blockage Detector" PINN**: Constrained only by the **Momentum Balance (Ergun Equation)**.
*   **Function**: When an anomaly is detected, the data is passed to all specialists. The diagnosis is made by observing which specialist's physics loss "complains the loudest." If the Leak Detector's loss spikes, the fault is a leak.

**Conclusion**: The hybrid approach is superior. It avoids the training nightmare of a single complex model while creating a modular, interpretable, and scalable system that is robust for both fault detection (Stage 1) and diagnosis (Stage 2).

---

## Section 2: Practical Guide - Training a Specialist PINN

Here is a concrete path for training one of the specialist models from the ensemble. We will use the **"Leak Detector" PINN** as our example.

### Step 1: Data Preparation and Formatting

This is the most critical step.

**1. Generate Raw Data:**
*   **Normal Dataset**: Use `PSASimulator.jl` to generate simulations under various normal operating conditions.
*   **Faulty Dataset**: Generate a separate set of simulations where a leak fault has been intentionally injected. This is for testing.

**2. Extract Variables:**
*   From the simulation `traj` object, extract the variables appearing in the **Total Mass Balance** equation: `P(z, τ)`, `T(z, τ)`, `x₁(z, τ)`, `x₂(z, τ)`, and the calculated velocity `v(z, τ)`.

**3. Structure the Data for the PINN:**
*   A PINN is trained on a "point cloud." You must create three sets of these collocation points:
    *   **Initial Condition (IC) Points**: The state of the system at `τ = 0`. A list of `[z, 0]` points with their known `[P, T, ...]` values.
    *   **Boundary Condition (BC) Points**: The state at the column inlet (`z = 0`) and outlet (`z = 1`) for all time `τ`. A list of `[0, τ]` and `[1, τ]` points with their known values.
    *   **Physics Collocation Points**: A large number of random `[z, τ]` coordinate points scattered throughout the entire simulation domain. These points do not need measured outcomes; their purpose is to enforce the PDE.

### Step 2: Building and Training the PINN

**1. Network Architecture:**
*   A standard feed-forward neural network (Multi-Layer Perceptron).
*   **Input**: 2 neurons for a `(z, τ)` coordinate.
*   **Hidden Layers**: Several layers with `tanh` activation.
*   **Output**: Neurons to predict the state variables: `P_pred, T_pred, x₁_pred, x₂_pred`.

**2. The Physics-Informed Loss Function:**
*   The total loss is a weighted sum: `Loss = w_ic*Loss_ic + w_bc*Loss_bc + w_phys*Loss_phys`.
*   `Loss_ic` & `Loss_bc`: The Mean Squared Error (MSE) between the PINN's predictions and the known values at the initial and boundary collocation points.
*   `Loss_phys`: The core of the PINN.
    1.  For all Physics Collocation Points, get the network's output `(P_pred, ...)`.
    2.  Use **automatic differentiation** to compute the necessary derivatives (`∂P_pred/∂τ`, etc.) from the network's output.
    3.  Plug these predictions and derivatives into the **Total Mass Balance equation**.
    4.  The result is the "residual." This loss term is the MSE of the residual, forcing it towards zero.

**3. Training Process:**
*   Train the network using the collocation points from the **Normal Dataset** only.
*   The optimizer adjusts the network's weights to minimize the total combined loss, resulting in a network that has learned a solution to the PDE.

### Step 3: Using the Trained PINN for Fault Detection

1.  Take data from the **Faulty (Leaky) Dataset**.
2.  Feed the `(z, τ)` coordinates into your trained PINN.
3.  Calculate **only the physics loss** (the residual of the Total Mass Balance equation).
4.  **Result**: Because the data comes from a leaky system, it will be inconsistent with the physical law the PINN has mastered. The equation residual will be **large and non-zero**. This spike in the physics loss is your fault detection signal.
