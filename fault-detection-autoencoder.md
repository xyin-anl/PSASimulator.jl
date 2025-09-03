### The Goal

The core idea is to train a model that learns the "signature" of your PSA system's **normal behavior** so well that it can immediately tell when something deviates from that normal signature, even if it has never seen that specific deviation before.

---

### The Recipe: Fault Detection with an Autoencoder

#### **Step 1: Prepare Your Ingredients (The Data)**

You need two main batches of data, generated using the simulator.

**Ingredient 1: A Large Dataset of *Normal* Operations**
*   **What it is**: This is the most important ingredient. It's a comprehensive dataset representing the system running correctly under a wide range of acceptable conditions.
*   **How to create it**:
    *   Use the `psacycle` function to run hundreds or thousands of simulations.
    *   For each simulation, slightly vary the `process_vars` (like feed pressure `P_0`, adsorption time `t_ads`, etc.) within their expected *normal operating ranges*.
    *   The `simulation_data` in the `demo` folder is a small-scale example of this. You need to generate a much larger and more diverse set.
*   **Why you need it**: The autoencoder will be trained **exclusively** on this data. It needs to learn the full spectrum of what "normal" looks like.

**Ingredient 2: A Labeled Dataset of *Faulty* Operations**
*   **What it is**: This is your "testing set." It's a collection of simulations where you have intentionally introduced specific faults.
*   **How to create it**:
    *   Use the fault injection techniques we discussed (and that are documented in `fault-injection.md`).
    *   For each type of fault you want to be able to detect (e.g., "Stuck Valve," "Column Leak," "Adsorbent Degradation"), run a number of simulations.
    *   Keep this data separate and label it. For example, all simulations with a stuck valve are labeled "Fault Type 1," all simulations with a leak are "Fault Type 2," and so on.
*   **Why you need it**: This data will be used to **validate** your system. You will use it to check if your trained model can successfully spot the anomalies you created. **Do not use this data for training the autoencoder.**

#### **Step 2: The Cooking Method (Training the Autoencoder)**

1.  **Choose Your Model**: An autoencoder is a type of neural network with two parts: an "encoder" that compresses the input data into a smaller representation, and a "decoder" that tries to reconstruct the original data from that compressed version.

2.  **Train the Model**:
    *   Feed your entire **Normal Operations Dataset** (Ingredient 1) to the autoencoder.
    *   The model's goal is to minimize the **reconstruction error**, which is the difference between the original input data and the reconstructed output.
    *   Because it only ever sees normal data, the model becomes an expert at reconstructing normal PSA cycles. Conversely, it will be very bad at reconstructing faulty data that doesn't fit the patterns it has learned.

#### **Step 3: The Taste Test (Detecting Faults)**

1.  **Set a Threshold**:
    *   After training, pass some of your normal data (a validation subset of Ingredient 1 that the model wasn't trained on) through the autoencoder and calculate the reconstruction error for each sample.
    *   This will give you a distribution of errors for "normal" data. You can then pick a threshold (e.g., the 99th percentile of these errors).
    *   **Any reconstruction error above this threshold will be flagged as an anomaly.**

2.  **Deploy and Detect**:
    *   Now, you can feed new, unseen data into the model in real-time (or from your **Faulty Operations Dataset** for testing).
    *   Calculate the reconstruction error for this new data.
    *   **If the error is *above* your threshold, you have detected a fault!** If it's below, the system is running normally.

### Summary of the Workflow

**Generate Data** -> **Train Autoencoder on Normal Data** -> **Determine Anomaly Threshold** -> **Detect Faults in New Data**.

This approach is powerful because the autoencoder can detect faults it has never been explicitly trained on. As long as the fault produces a signal that deviates from the learned normal behavior, the reconstruction error will spike, and your system will raise an alarm.