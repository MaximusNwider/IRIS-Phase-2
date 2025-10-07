% File: main_ris_cocross_from_s2p.m
% Dual-pol RIS: build Zs from geometry + measured varactor S2P,
% compute reflectivity matrix R, and plot two traces (co, cross)
% for magnitude (dB) and phase (deg).

clear; clc;

%% ---- USER: point to the varactor S2P files ----
file_x_10V = 'SMV1408-219-S-PAR-Vr-10V.s2p';   % x-branch (reverse 10 V)
file_y_0V  = 'SMV1408-219-S-PAR-Vr-0V.s2p';    % y-branch (reverse 0 V)

assert(exist(file_x_10V,'file')==2, 'Missing S2P: %s', file_x_10V);
assert(exist(file_y_0V, 'file')==2, 'Missing S2P: %s', file_y_0V);

%% ---- Physical constants, materials, geometry (as in your topology code) ----
mu0   = 4*pi*1e-7;
eps0  = 8.854187817e-12;
eta0  = sqrt(mu0/eps0);

er    = 2.2;          % F4B220
tand  = 1e-3;
h     = 2e-3;
sigma_cu = 5.8e7;

% Parasitic patches (per axis)
Lp_x  = 4.0e-3;  Wp_x  = 3.0e-3;
Lp_y  = 4.0e-3;  Wp_y  = 3.0e-3;

% Gap / grid geometry
g        = 0.50e-3;
wedge_x  = 6.0e-3;  wedge_y  = 6.0e-3;
Ledge_x  = 6.0e-3;  Ledge_y  = 6.0e-3;

% Via / loop return
d_via   = 0.60e-3;
lcur_x  = 6.0e-3;  wcur_x  = 1.0e-3;
lcur_y  = 6.0e-3;  wcur_y  = 1.0e-3;

% ---- Cross-pol coupling: DO NOT neglect ----
% Use CPW form to estimate mutual capacitance Cm from geometry (gd,wd,Ld)
gd   = 1.0e-3; wd = 2.0e-3; Ld = 6.0e-3;   % set from your layout
eeff_d = hammerEffectivePermittivity(er, h, wd);
Kfun = @(k) ellipke(k.^2);
kd  = gd./(gd + 2*wd);  kpd = sqrt(1 - kd.^2);
Ccpd = 4*eps0*eeff_d * (Kfun(kd)/Kfun(kpd));   % per-unit-length
Cm   = Ccpd * Ld;                               % mutual C (non-zero)

% Parasitic & grid shunt caps per axis (Hammerstad + CPW)
eeff_x = hammerEffectivePermittivity(er, h, Wp_x);
eeff_y = hammerEffectivePermittivity(er, h, Wp_y);
dL_x   = hammerDeltaL(eeff_x, h, Wp_x);
dL_y   = hammerDeltaL(eeff_y, h, Wp_y);
Ap_eff_x = (Lp_x + 2*dL_x) * (Wp_x + 2*dL_x);
Ap_eff_y = (Lp_y + 2*dL_y) * (Wp_y + 2*dL_y);
Cp_x     = eps0*eeff_x * Ap_eff_x / h;
Cp_y     = eps0*eeff_y * Ap_eff_y / h;

Ccp_x = cpwCprime(eps0, eeff_x, wedge_x, g);
Ccp_y = cpwCprime(eps0, eeff_y, wedge_y, g);
Cg_x  = Ccp_x * Ledge_x;
Cg_y  = Ccp_y * Ledge_y;

%% ---- Frequency sweep (match your plotting template semantics) ----
fvec = linspace(3.2e9, 3.8e9, 401);   % Hz
wvec = 2*pi*fvec;

%% ---- Varactor SERIES impedance from S2P ONLY (no datasheet Ls/ESR) ----
Zvar_x = loadVaractorZseries(file_x_10V, fvec);   % x-branch (10 V)
Zvar_y = loadVaractorZseries(file_y_0V,  fvec);   % y-branch ( 0 V)

%% ---- Build per-axis surface impedances and cross terms ----
Zxx = computeZaxis(fvec, Zvar_x, struct('er',er,'tand',tand,'h',h,'sigma',sigma_cu), ...
    struct('Lp',Lp_x,'Wp',Wp_x,'Ledge',Ledge_x,'wedge',wedge_x,'lcur',lcur_x,'wcur',wcur_x,'d_via',d_via,'g',g,'Cp',Cp_x,'Cg',Cg_x));
Zyy = computeZaxis(fvec, Zvar_y, struct('er',er,'tand',tand,'h',h,'sigma',sigma_cu), ...
    struct('Lp',Lp_y,'Wp',Wp_y,'Ledge',Ledge_y,'wedge',wedge_y,'lcur',lcur_y,'wcur',wcur_y,'d_via',d_via,'g',g,'Cp',Cp_y,'Cg',Cg_y));

Zxy = 1 ./ (1i*wvec*Cm);     % cross-polar branch from mutual C
Zyx = Zxy;                   % reciprocity at broadside

%% ---- Reflectivity matrix entries (broadside) ----
eta1 = eta0; eta2 = eta0;
[Rxx,Rxy,Ryx,Ryy] = reflectivityFromZs(fvec, Zxx, Zyy, Zxy, Zyx, eta1, eta2);

%% ---- Reduce to two traces: co & cross (use 'avg' across x/y) ----
[co_mag_db, xpol_mag_db, co_phase_deg, xpol_phase_deg] = ...
    computeCoCross(fvec, Rxx, Rxy, Ryx, Ryy, 'avg');

%% ---- Figure 1: Magnitude (dB) ----
plotMagnitudeCoCross(fvec, co_mag_db, xpol_mag_db, ...
    'Co/Cross Reflectivity (V_x=10V S2P, V_y=0V S2P)', ...
    [3.2 3.8], [-10 1]);

%% ---- Figure 2: Phase (degrees) ----
plotPhaseCoCross(fvec, co_phase_deg, xpol_phase_deg, ...
    'Co/Cross Phase (V_x=10V S2P, V_y=0V S2P)', ...
    [3.2 3.8], [-360 0]);
