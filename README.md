# Multi-type Spike Model

The **Multi-type Spike Model** is a BEAST 2 clock model for modelling **punctuated evolution** under structured birth–death processes.

It extends the Gamma Spike model of Douglas et al. (2025) by allowing:

- **State‑dependent speciation and extinction** 
- **Time varying (skyline) phylodynamic parameters**

> **Note:** Unlike the original Gamma Spike model, the Multi-type Spike Model **does not have the option to estimate the number or timing of hidden events** (“stubs”) explicitly. It integrates over the number of hidden speciation events numerically. However, all other features of the Gamma Spike model are retained.

---

## Model overview

The evolutionary distance along branch $e$ is:

$$
d^e = r^e \ \tau^e + s^e
$$

where:

- $r^e$ is the **branch rate**  
- $\tau^e$ is the **branch length** (time)  
- $s^e$ is the **spike** on branch $e$
  
The term $r^e \tau^e$ represents gradual change along the branch, while $s^e$ represents the contribution of punctuated evolution (“spikes”) along the branch.

---

## Spike Distribution

Each branching event contributes an instantaneous burst of substitutions.
If branch $e$ belongs to type $i$ and has $N_i^e$ branching events
(observed + hidden), then its total spike magnitude is:

$$
s_i^e \sim S^\mu_i \cdot \mathrm{Gamma}\\left(
\text{shape} = S^\alpha \cdot N_i^e\ , 
\text{scale} = \frac{1}{S^\alpha}
\right).
$$

---

## Model Averaging
As in the original Gamma Spike model, the posterior support for punctuated evolution can be assessed directly from the data. The clock indicator (`useSpikeModel`) is assigned the prior

$\mathbb{I} ∼ Bernoulli(0.5)$,

such that the posterior probability of $\mathbb{I} = 1$ provides evidence for punctuated evolution relative to a relaxed (or strict) clock model.

---

### Key Parameters

| Parameter | Meaning | Notes |
|----------|---------|-------|
| $S^\mu$ or $S^\mu_i$ | Spike mean | Expected substitutions per site **per branching event**. Optionally type‑specific in multi‑type mode. |
| $S^\alpha$ or $S^\alpha_i$ | Spike shape | Controls variance in spike sizes. Larger $S^\alpha$ $\rightarrow$ more uniform spikes. Optionally type‑specific in multi‑type mode. |
| $\mathbb{I}$ or $\mathbb{I}_i$ | Spike model indicator | Bernoulli prior; posterior support quantifies evidence for punctuated evolution. Optionally type‑specific in multi‑type mode. |


> In single‑type models, all lineages share the same $S^\mu$ and $S^\alpha$.  
> In multi‑type models, each type $i$ can be associated with its own $S^\mu_i$ and/or $S^\alpha_i$ and/or $\mathbb{I}_i$ parameters, enabling testing of type‑specific differences in punctuated evolution.

---

## Testing for Differences in Punctuated Evolution Between Types

Type‑specific spike parameters allow for hypothesis testing:

- **Trait‑dependent evolution:**  
  Are certain traits (e.g., body size, ecological niche) associated with larger spikes?

- **Geography‑dependent evolution:**  
  Do lineages in region A exhibit larger spikes than those in region B?

---

## Dependencies

- **BDMM-Prime**  

- **BEAST 2.7**

---
## Installation

### Building from source

To build **MultiTypeSpikeModel** from source, you will need the following installed and in your execution path:

- OpenJDK version 17 or greater
- A recent version of OpenJFX (JavaFX SDK)
- The Apache Ant build system

The build expects `beast2`, `BeastFX`, and `BDMM-Prime` checked out as sibling directories relative to this repository:

```bash
git clone [https://github.com/CompEvol/beast2.git](https://github.com/CompEvol/beast2.git)
git clone [https://github.com/CompEvol/BeastFX.git](https://github.com/CompEvol/BeastFX.git)
git clone [https://github.com/tgvaughan/BDMM-Prime.git](https://github.com/tgvaughan/BDMM-Prime.git)
git clone [https://github.com/EwanCiuffi/MultiTypeSpikeModel.git](https://github.com/EwanCiuffi/MultiTypeSpikeModel.git)
```

Once installed, issue the following command from the root directory of this repository:

```bash
JAVA_FX_HOME=/path/to/openjfx/lib ant

JAVA_FX_HOME=/path/to/openjfx/lib ant install
```

---

## Citation

If you use this model, please cite:

- Manuscript in preparation (TBC)  

- **Gamma Spike Model**
  Douglas, J., Bouckaert, R., Harris, S., Carter, C., & Wills, P. (2025).  
  *Evolution is coupled with branching across many granularities of life.*  
  Proc. R. Soc. B. 292:20250182  
  http://doi.org/10.1098/rspb.2025.0182
  
- **BDMM-Prime**
  Vaughan, T. G., & Stadler, T. (2025).  
  *Bayesian phylodynamic inference of multi-type population trajectories using genomic data.*  
  *Molecular Biology and Evolution* 42: msaf130  
  https://doi.org/10.1093/molbev/msaf130

- **BEAST 2**
  Bouckaert, R., Vaughan, T. G., Barido-Sottani, J., Duchêne, S., Fourment, M., Gavryushkina, A., ... & Drummond, A. J. (2019).  
  *BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis.*  
  *PLoS Computational Biology* 15(4): e1006650  
  https://doi.org/10.1371/journal.pcbi.1006650
