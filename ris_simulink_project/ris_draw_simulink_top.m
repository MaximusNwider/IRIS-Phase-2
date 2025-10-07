% Builds R(jw) from S2P → fits MIMO ss → draws Simulink model.
clear; clc;

% --- USER: set your S2P files (place in this folder) ---
file_x_10V = fullfile(pwd,'SMV1408-219-S-PAR-Vr-10V.s2p');   % x-branch (10 V)
file_y_0V  = fullfile(pwd,'SMV1408-219-S-PAR-Vr-0V.s2p');    % y-branch ( 0 V)

assert(exist(file_x_10V,'file')==2, 'Missing S2P: %s', file_x_10V);
assert(exist(file_y_0V, 'file')==2, 'Missing S2P: %s', file_y_0V);

% --- Materials / geometry (as in the topology) ---
P.mu0 = 4*pi*1e-7;  P.eps0 = 8.854187817e-12;  P.eta0 = sqrt(P.mu0/P.eps0);
P.er=2.2; P.tand=1e-3; P.h=2e-3; P.sigma=5.8e7;

G.Lp_x=4.0e-3; G.Wp_x=3.0e-3; G.Lp_y=4.0e-3; G.Wp_y=3.0e-3;
G.g=0.50e-3; G.wedge_x=6.0e-3; G.wedge_y=6.0e-3; G.Ledge_x=6.0e-3; G.Ledge_y=6.0e-3;
G.d_via=0.60e-3; G.lcur_x=6.0e-3; G.wcur_x=1.0e-3; G.lcur_y=6.0e-3; G.wcur_y=1.0e-3;
G.gd=1.0e-3; G.wd=2.0e-3; G.Ld=6.0e-3;  % cross-pol geometry for Cm

% Sweep (GHz semantics from your plots)
fvec = linspace(3.2e9, 3.8e9, 401);

% --- Build R(jw) as FRD (2x2) using S2P-only varactor branches ---
[R_frd, Rxx,Rxy,Ryx,Ryy] = ris_build_reflectivity_frd(fvec, file_x_10V, file_y_0V, P, G);

% --- Fit a MIMO state-space model R(s) ---
order = 10;  % try 8–16
sysR  = ris_fit_ss(R_frd, order);

% --- Draw Simulink model ---
mdlName = 'RIS_UnitCell_LTI';
ris_make_simulink(sysR, mdlName);
open_system(mdlName);
