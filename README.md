# README — Dual-Polarized Varactor-Loaded RIS (S2P-driven)

This package models a dual-polarized RIS **unit cell** using only the **measured varactor S-parameter files** and produces two figures that match the publication’s plotting semantics:

1) **Magnitude (dB)**: **copolar** (solid) and **cross-polar** (dashed)  
2) **Phase (degrees)**: **copolar** (solid) and **cross-polar** (dashed)

The varactor behavior is taken **exclusively** from the S2P files (x-branch uses the **10 V** file; y-branch uses the **0 V** file). Nothing else (datasheet \(L_s\), \(Q\), \(C(V)\), etc.) is injected into the varactor branch. **Cross-polarization is never neglected** (a non-zero \(C_m\) is used).

---

## File layout (one function per file)

```
main_ris_cocross_from_s2p.m          % Entry point (edit two S2P paths, then run)
readS2P.m                             % Touchstone reader (dB/MA/RI; HZ/K/M/G units)
S2Z_2port.m                           % 2-port S → Z (robust to 2x2xN or Nx2x2)
Z2SeriesTwoTerminal.m                 % Z11+Z22−Z12−Z21 reduction (series 2-terminal)
loadVaractorZseries.m                 % Uses the three above to return Zvar(f)
hammerEffectivePermittivity.m         % Hammerstad ε_eff(W,h,ε_r)
hammerDeltaL.m                        % Hammerstad ΔL(W,h,ε_eff)
cpwCprime.m                           % CPW per-unit-length capacitance 4ε₀ε_eff K/K’
computeZaxis.m                        % Builds Zxx or Zyy from Zvar + board physics
reflectivityFromZs.m                  % Closed-form 2x2 reflectivity entries
computeCoCross.m                      % Reduce R to two traces (co, cross) + phases
plotMagnitudeCoCross.m                % Figure 1 (dB): co (solid), cross (dashed)
plotPhaseCoCross.m                    % Figure 2 (deg): co (solid), cross (dashed)
```

---

## Quick start

1) Place all `.m` files in one folder on your MATLAB path.  
2) Open **`main_ris_cocross_from_s2p.m`** and set:
   - `file_x_10V` → your **10 V** S2P (x-branch)  
   - `file_y_0V`  → your **0 V** S2P (y-branch)
3) Run `main_ris_cocross_from_s2p`.  
   You will get two figures: **Magnitude (dB)** and **Phase (deg)** versus **3.2–3.8 GHz**.

---

## Chronological call flow (what runs, why it exists, and how it contributes)

### 1) `main_ris_cocross_from_s2p.m` (entry point)
**What it does**
- Defines **physical constants** \((\mu_0,\ \varepsilon_0,\ \eta_0)\), **substrate** (ε_r, tanδ, h), **conductor** (σ_Cu), and **geometry** for both polarizations (patch dims, gap geometry, return-path via/loop).
- Computes **parasitic shunt capacitors** per polarization:
  - Patch capacitance \(C_{p,\{x,y\}}\) using **Hammerstad** fringing (via `hammerEffectivePermittivity`, `hammerDeltaL`).
  - Gap/grid capacitance \(C_{g,\{x,y\}}\) using **CPW K/K’** (via `cpwCprime`).
- Builds a **non-zero cross-polar mutual** \(C_m\) (also via CPW form) — ensuring **cross-pol is included**.
- Defines the **frequency sweep** (default 3.2–3.8 GHz → matches the template bounds).
- **Loads measured varactor behavior** for each polarization:
  - Calls `loadVaractorZseries` for **x** with the **10 V** file and for **y** with the **0 V** file → returns **measured** series **Zvar(f)**.
- Builds **surface impedances** \(Z_{xx}(f)\) and \(Z_{yy}(f)\) using `computeZaxis` (details below).
- Builds **cross terms** \(Z_{xy}(f)=Z_{yx}(f)=\frac{1}{j\omega C_m}\).
- Calls `reflectivityFromZs` to obtain **\(R_{xx},R_{xy},R_{yx},R_{yy}\)**.
- Reduces to **two traces** via `computeCoCross` (copolar & cross-polar; both **magnitude in dB** and **phase in degrees**).
- Calls `plotMagnitudeCoCross` and `plotPhaseCoCross` to produce the two figures.

**Why it matters**
- This script is the **scenario definition** and **orchestrator**. It ensures the model uses the S2P files **as the source of truth** for varactor behavior and that **cross-pol** is included. All other blocks are purely **unit-cell physics** (board return, copper loss, grid/patch shunts).

---

### 2) `loadVaractorZseries.m`  →  **measured series \(Z_{\text{var}}(f)\)**

**What it does**
1. **Reads** the Touchstone file via `readS2P` (handles `# HZ|kHz|MHz|GHz S DB|MA|RI R <Z0>`).
2. Converts S→Z with `S2Z_2port` for each frequency sample.
3. Reduces the 2-port network to a **two-terminal series impedance**:
   \[
   Z_\mathrm{series}(f) = Z_{11}(f) + Z_{22}(f) - Z_{12}(f) - Z_{21}(f).
   \]
4. **Interpolates** this measured \(Z_\mathrm{series}(f)\) onto the analysis grid `fvec` using PCHIP.

**Why it matters**
- This is the **only place** where varactor behavior enters the model. By using the measured \(Z_\mathrm{series}\), we **avoid double-counting** package/bond \(L_s\) or ESR. Everything thereafter treats the varactor as a **black-box series impedance**.

---

### 3) `readS2P.m`  →  **Touchstone parser**

**What it does**
- Parses the header to determine **frequency units**, **format** (DB/MA/RI), and **reference impedance** \(Z_0\).
- Parses each numeric row into **\(S_{11}, S_{21}, S_{12}, S_{22}\)** as **complex** numbers.
- Returns \(f\) (Hz) and \(\mathbf S(f)\) in a consistent shape.

**Why it matters**
- Ensures we faithfully interpret the manufacturer data. The entire chain relies on these values being correct.

---

### 4) `S2Z_2port.m`  →  **convert S to Z**

**What it does**
- Uses the standard relation:
  \[
  \mathbf Z = Z_0(\mathbf I+\mathbf S)(\mathbf I-\mathbf S)^{-1}
  \]
  per frequency sample.
- Accepts either **2×2×N** or **N×2×2** and returns **2×2×N**.

**Why it matters**
- Many vendors and tools differ in how they shape arrays; this file is a robust adapter.

---

### 5) `Z2SeriesTwoTerminal.m`  →  **series reduction**

**What it does**
- Computes \(Z_{11}+Z_{22}-Z_{12}-Z_{21}\), which is the **series two-terminal** equivalent seen between the device pins.

**Why it matters**
- The RIS branch is a **series** varactor path in the unit cell; this number is the correct series element to embed into the board model.

---

### 6) `hammerEffectivePermittivity.m` & `hammerDeltaL.m`  →  **patch parasitics**

**What they do**
- Hammerstad/Wheeler formulas to compute **effective permittivity** and **fringe length extension** for a microstrip patch of width \(W\) on thickness \(h\).
- From these, `main` computes **effective patch area** and **\(C_p\)** per axis:
  \[
  C_{p} = \varepsilon_0\,\varepsilon_\text{eff}\,\frac{A_\text{eff}}{h}.
  \]

**Why it matters**
- These are the **shunt** capacitive paths that load the series branch and strongly shape \(R(f)\).

---

### 7) `cpwCprime.m`  →  **gap/grid shunt capacitance**

**What it does**
- Computes CPW **per-unit-length** capacitance:
  \[
  C' = 4\varepsilon_0 \varepsilon_\text{eff} \frac{K(k)}{K(k')},\quad
  k=\frac{g}{g+2w_\text{edge}}
  \]
- Multiplies by engaged length \(L_\text{edge}\) to get **\(C_g\)** per axis.

**Why it matters**
- This captures the **grid gap** loading (another shunt path). It sets bandwidth and the dispersion of the surface impedance.

---

### 8) `computeZaxis.m`  →  **surface impedance per polarization** \(Z_{xx}\) or \(Z_{yy}\)

**What it does**
- Takes the **measured** \(Z_{\text{var}}(f)\) for one axis and adds **only the board contributions**:
  - **Return inductance** \(L_\text{via}+L_\text{loop}\).
  - **Copper loop loss** via skin effect \(R_\text{loop}(f)\).
- Builds the **series branch**:  
  \(Z_\text{series}(f) = Z_{\text{var}}(f) + j\omega L_\text{ret} + R_\text{loop}(f)\)
- Builds **shunt admittances** with dielectric loss:
  \[
  Y_p = j\omega C_p(1 - j\tan\delta),\quad
  Y_g = j\omega C_g\!\left(1 - j\frac{\tan\delta}{2}\right)
  \]
- Returns:
  \[
  Z_{uu}(f) = \frac{1}{\,1/Z_\text{series} + Y_p + Y_g\,}\quad (u\!\in\!\{x,y\})
  \]

**Why it matters**
- This is the **unit-cell physics** layer. It converts the device’s measured series impedance into a **sheet impedance** per polarization for the RIS boundary model.

---

### 9) (In `main`) cross terms \(Z_{xy}=Z_{yx}=\frac{1}{j\omega C_m}\)

**What it does**
- Models polarization coupling by a **mutual capacitance \(C_m\)** (non-zero) computed from gap geometry using the same CPW approach.

**Why it matters**
- Guarantees **cross-polarization is included** (no matter how small), as required.

---

### 10) `reflectivityFromZs.m`  →  **R-matrix entries**

**What it does**
- Uses the 2×2 sheet impedance \(\mathbf Z_s\) with free-space impedances \(\eta_1,\eta_2\) (broadside → both \(\eta_0\)) to compute:
  \[
  \mathbf R = (\mathbf Z_s - \operatorname{diag}(\eta_1,\eta_2))\,
              (\mathbf Z_s + \operatorname{diag}(\eta_1,\eta_2))^{-1}.
  \]
- Returns the four entries \(R_{xx},R_{xy},R_{yx},R_{yy}\) as vectors over frequency.

**Why it matters**
- This is the **observable** quantity the plots use (co/cross magnitudes and phases).

---

### 11) `computeCoCross.m`  →  **two traces from R** (magnitude dB + phase deg)

**What it does**
- Reduces the 2×2 R-matrix into **copolar** and **cross-polar** traces. Default `'avg'` mode yields:
  - **Co-pol magnitude**: RMS across x/y incidence,
    \[
    |R|_\text{co} = \sqrt{\tfrac{1}{2}(|R_{xx}|^2 + |R_{yy}|^2)}
    \]
  - **Cross-pol magnitude**: RMS across x→y and y→x,
    \[
    |R|_\text{xpol} = \sqrt{\tfrac{1}{2}(|R_{xy}|^2 + |R_{yx}|^2)}
    \]
  - Phases: phasor-averaged `angle(0.5*(Rxx+Ryy))` and `angle(0.5*(Rxy+Ryx))`, then `unwrap` and convert to degrees.  
- Alternatives:
  - `'xinc'`: co=\(|R_{xx}|\), cross=\(|R_{yx}|\)  
  - `'yinc'`: co=\(|R_{yy}|\), cross=\(|R_{xy}|\)  
  - `'max'`: worst-case across axes (element-wise max of magnitudes)

**Why it matters**
- Produces exactly **two** curves per bias (co and cross), ready for the requested plot style.

---

### 12) `plotMagnitudeCoCross.m`  →  **Figure 1 (Magnitude, dB)**

**What it does**
- Plots **co-pol** magnitude (solid) and **cross-pol** magnitude (dashed) vs **GHz**.
- Axes, labels, line styles follow the provided template semantics. Optional `xlim`/`ylim` arguments.

**Why it matters**
- This is the primary figure for **reflectivity magnitude**, aligned with your stylistic/semantic expectations.

---

### 13) `plotPhaseCoCross.m`  →  **Figure 2 (Phase, degrees)**

**What it does**
- Plots **co-pol** phase (solid) and **cross-pol** phase (dashed) vs **GHz**, using `unwrap` to keep curves continuous. Optional `xlim`/`ylim`.

**Why it matters**
- Complements the magnitude plot, matching the template’s appearance and units.

---

## Modeling guarantees & assumptions

- **Varactor behavior**: derived **only** from the S2P files (x=10 V, y=0 V).  
  No datasheet \(L_s\), \(Q\), or junction models are added (prevents double-counting).
- **Board physics**: return-path inductance and copper loss are **outside** the packaged device and therefore **are** added in `computeZaxis`.
- **Cross-pol**: \(C_m\neq 0\) by construction. You can tune `gd, wd, Ld` to your layout; or replace \(Z_{xy}\) with a more detailed coupling if you have one.
- **Interpolation**: PCHIP on frequency for the S2P → analysis grid transfer; safe and shape-preserving.

---

## Extending & generalizing

- **Different biases**: swap the filenames in `main` to any S2P pair; the pipeline stays identical.  
- **Multiple overlays**: loop over several (x,y) bias pairs, call the same functions, and overlay traces with different colors/linestyles.
- **Non-broadside**: replace \(\eta_1,\eta_2\) as needed and/or extend \(\mathbf Z_s\) to include angle-dependent effects.
- **Alternative co/cross definitions**: call `computeCoCross` with `'xinc'`, `'yinc'`, or `'max'` as needed.

---

## Troubleshooting

- **“File not found”**: verify `file_x_10V`, `file_y_0V` paths.  
- **Size mismatch in plots**: ensure all returned vectors are the same length; the plotting helpers force column vectors and keep sizes aligned.  
- **Out-of-band extrapolation**: we use `'extrap'` in `interp1`. If you want to **avoid extrapolation**, clamp `fvec` to the S2P’s frequency span.

---

## What each function contributes (one line each)

- `main_ris_cocross_from_s2p` — Orchestrates the entire workflow; defines scenario; produces figures.  
- `readS2P` — Faithful decoding of Touchstone S-parameters.  
- `S2Z_2port` — Converts S→Z per frequency sample.  
- `Z2SeriesTwoTerminal` — Extracts the correct series 2-terminal varactor impedance from the 2-port network.  
- `loadVaractorZseries` — Returns **measured** varactor series \(Z(f)\) on the analysis grid (source of truth).  
- `hammerEffectivePermittivity`, `hammerDeltaL` — Patch fringing model → \(C_p\).  
- `cpwCprime` — CPW gap model → per-unit-length C; becomes \(C_g\) and \(C_m\).  
- `computeZaxis` — Forms the **surface impedance** \(Z_{xx}\) or \(Z_{yy}\) by embedding measured \(Z_\text{var}\) into the unit-cell network.  
- `reflectivityFromZs` — Computes observable **reflectivity matrix** entries \(R_{ij}(f)\).  
- `computeCoCross` — Reduces \(R\) to **co/cross** traces (dB & degrees) using a clear physical convention.  
- `plotMagnitudeCoCross`, `plotPhaseCoCross` — Render the two required figures with solid/dashed semantics and correct units.