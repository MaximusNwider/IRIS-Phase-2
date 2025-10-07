# RIS Unit Cell — Simulink LTI Project (S2P-driven)

This project programmatically builds a **Simulink model** of a dual-polarized RIS **unit cell** under the **LTI** assumption. It uses ONLY the measured varactor **S2P** files (x-branch = **10 V**, y-branch = **0 V**) to form the 2×2 frequency response **R(jω)**, fits a continuous-time **MIMO state-space** **R(s)**, and draws a Simulink diagram:

```
Ex_in ──┐
        ├─[Mux]──► [ State-Space (R(s)) ] ──► [Demux] ──► Ex_ref
Ey_in ──┘                                      │
                                               └────────► Ey_ref
```

## Requirements
- MATLAB with Control System Toolbox and Simulink
- Functions used: `frd`, `fitfrd`, `ss`, `balreal`, `minreal`, `ellipke`
- Two Touchstone files:
  - `SMV1408-219-S-PAR-Vr-10V.s2p` (x-branch @ 10 V)
  - `SMV1408-219-S-PAR-Vr-0V.s2p`  (y-branch @ 0 V)  ← **Provide this** in the project folder

## How to run
1. Place both S2P files in this folder (next to the `.m` files).
2. Open MATLAB and `cd` into this folder:
   ```matlab
   cd <path_to>/ris_simulink_project
   addpath(pwd);
   ```
3. Run the orchestrator:
   ```matlab
   ris_draw_simulink_top;
   ```
4. The script creates and saves **`RIS_UnitCell_LTI.slx`** (model name `RIS_UnitCell_LTI`).
   Open it with:
   ```matlab
   open_system RIS_UnitCell_LTI
   ```

## File map
- `ris_draw_simulink_top.m` — Orchestrates the build (read S2P → FRD → fit ss → draw Simulink)
- `ris_build_reflectivity_frd.m` — Builds 2×2 FRD `R_frd` (Rxx,Rxy,Ryx,Ryy) across `fvec`
- `ris_fit_ss.m` — Fits a continuous-time MIMO state-space `R(s)` to the FRD
- `ris_make_simulink.m` — Programmatically draws the Simulink diagram with a State-Space block
- `loadVaractorZseries.m` — Reads S2P, converts S→Z, reduces to series Z, interpolates to analysis grid
- `readS2P.m`, `S2Z_2port.m`, `Z2SeriesTwoTerminal.m` — Touchstone & network utilities
- `hammerEffectivePermittivity.m`, `hammerDeltaL.m`, `cpwCprime.m` — EM helpers (patch & CPW)
- `computeZaxis.m` — Builds per-axis surface impedance from measured Zvar + board contributions
- `reflectivityFromZs.m` — Closed-form reflectivity entries from the sheet impedance matrix

## Notes
- The varactor branches are taken **exclusively** from the S2P data (no datasheet Ls/Q/C(V) added).
- Cross-polarization **is included** via a mutual capacitance `Cm` derived from gap geometry.
- You can adjust the fitting **order** in `ris_draw_simulink_top.m` (8–16 is typical). Check the fit report.
- To avoid extrapolating outside S2P ranges, clamp the sweep `fvec` to the common frequency span.
