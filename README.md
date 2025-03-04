## **Modeling Synaptic GABA-B Currents in a Hippocampal Microcircuit**

### **Authors**
- Andrea F. Campos Pérez
- Regina A. Mejia Ortiz
- Allan I. Rico Becerra

*July, 2020.*

---

### **Overview**
This repository contains the **Julia-based computational model and simulation scripts** used to study the role of **GABAergic synaptic currents** in modulating neuronal firing patterns in hippocampal microcircuits. The study focuses on how **GABA<sub>B</sub>** currents influence the generation of **spike-and-wave oscillations (SW)** associated with absence epilepsy.

The repository includes:
- A Julia script implementing intrinsic and synaptic currents for a pyramidal neuron and an interneuron.
- A Jupyter Notebook that runs simulations with different synaptic current configurations.
- The full research study (PDF), detailing the methodology, results, and discussion. (In Spanish)

---

### **Background**
Epilepsy, one of the most common severe brain disorders, manifests through convulsive and absence seizures. Electroencephalogram (EEG) recordings from epilepsy patients and animal models reveal **spike-and-wave (SW) patterns** linked to thalamocortical circuits. However, recent findings suggest a synchronized interaction between the **hippocampus and thalamocortical networks** during absence seizures.

This study investigates the role of **GABA<sub>B</sub> receptor-mediated inhibition** in generating **SW patterns** in a **hippocampal microcircuit model**, comprising:
- A **pyramidal neuron** (excitatory)
- An **interneuron** (inhibitory)

We examine how different configurations of **GABA<sub>A</sub>** and **GABA<sub>B</sub>** currents influence neuronal firing dynamics.

---

### **Model Implementation**
The computational model is implemented in **Julia** and follows **McCormick & Huguenard (1992)** for intrinsic neuronal currents and **Destexhe (1998)** for synaptic currents.

**Key Features:**
- **Voltage-dependent intrinsic currents**:  
  - Sodium currents: **I<sub>Na</sub>, I<sub>Nap</sub>**  
  - Calcium currents: **I<sub>T</sub>, I<sub>L</sub>**  
  - Potassium currents: **I<sub>K</sub>, I<sub>A</sub>, I<sub>K2</sub>, I<sub>C</sub>**  
  - Hyperpolarization-activated current: **I<sub>H</sub>**
- Synaptic currents: **AMPA**, **GABA<sub>A</sub>**, **GABA<sub>B</sub>**
- Presynaptic neurotransmitter release dynamics
- Euler integration method with **0.01 ms timestep**

### **Files in This Repository**
| File | Description |
|------|------------|
| `Fx_model.jl` | Julia script defining neuron dynamics, intrinsic and synaptic currents |
| `Sim_ModEpi.ipynb` | Jupyter Notebook running the simulation with different current conditions |
| `GABA_B_Synaptic_Current_Hippocampal_Model_Study.pdf` | Full research study detailing methodology, results, and analysis (in Spanish) |

---

### **Simulation Results**
- **Baseline Activity:** Without GABAergic currents, neurons exhibit **tonic firing**.
- **GABA<sub>A</sub> & GABA<sub>B</sub> Present:** Pyramidal neurons show **irregular burst firing**, synchronized with interneuron activity.
- **GABA<sub>A</sub> Absent, GABA<sub>B</sub> Present:** **Regular, periodic bursts**, characteristic of SW oscillations in absence seizures.
- **GABA<sub>A</sub> Present, GABA<sub>B</sub> Absent:** **Irregular burst patterns**, lacking SW-like periodicity.

These findings support the hypothesis that **GABA<sub>B</sub> currents play a key role in generating periodic burst firing, a hallmark of SW oscillations**.

---

### **How to Run the Simulation**
1. Install **Julia** (version 1.3.1 or later).
2. Clone this repository:
   ```sh
   git clone https://github.com/andreaf-campos/GABAB_Hippocampal_Model.git
   cd GABAB_Hippocampal_Model
   ```
3. Open the **Jupyter Notebook**:
   ```sh
   jupyter notebook Sim_ModEpi.ipynb
   ```
4. Run the notebook to explore different synaptic configurations.

---

### **Future Directions**
- Expand the model to include **hippocampal networks** instead of a two-neuron microcircuit.
- Integrate the hippocampus into a **thalamocortical loop model** for SW oscillation analysis.
- Simulate the effect of **GABAergic drugs** to explore **therapeutic interventions** for absence epilepsy.

---

### **Main References**
- Destexhe, A. (1998). *Spike-and-wave oscillations based on the properties of GABA(B) receptors*. Journal of Neuroscience, 18(21), 9099–9111.
- McCormick, D. A., & Huguenard, J. R. (1992). *A model of the electrophysiological properties of thalamocortical relay neurons*. Journal of Neurophysiology, 68(4), 1384–1400.

---

### **License**
This project is licensed under the MIT License.

If you use this work in academic research or derivative studies, please cite the original authors:

*Campos Pérez, A. F., Mejia Ortiz, R. A., & Rico Becerra, A. I. (2020).*
