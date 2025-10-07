%% ================================================================
%  Dual-Polarized Varactor-Loaded RIS — Eigenvalues of Reflectivity Matrix
%  1) Paper constants + canonical geometry/electromagnetics (fill Table I)
%  2) Symbolic R(f) with explicit physical forms (f is symbolic)
%  3) Numeric sweep (3.3–3.7 GHz) and |eigs(R)| in dB, Vx=10 V, Vy=2.5 V
%  ================================================================

%% 1) Initialization — materials, device, geometry
% Physical constants
mu0   = 4*pi*1e-7;
eps0  = 8.854187817e-12;
eta0  = sqrt(mu0/eps0);

% Substrate: F4B220 (publication)
er    = 2.2;
tand  = 1e-3;
h     = 2e-3;            % 2 mm

% Conductor (copper)
sigma_cu = 5.8e7;        % S/m

% Varactor: Skyworks SMV1408 (canonical junction law; replace with S2P if available)
Cj0   = 2.0e-12;         % F (datasheet/extracted)
phi   = 0.7;             % V
m     = 0.5;             % grading coefficient
Ls    = 0.5e-9;          % H  (pkg/bond; refine from S2P)
Qx_nom = 60;             % placeholder Q(Vx,f); replace with measured Q
Qy_nom = 60;             % placeholder Q(Vy,f)

% 1-bit biases (publication)
Vy    = 10.0;            % V (x-pol branch)
Vx    =  2.5;            % V (y-pol branch)

% -------- Geometry (SET to Table I for exact replication) --------
% Parasitic patches (per axis)
Lp_x  = 4.0e-3;          % m  (placeholder)
Wp_x  = 3.0e-3;          % m  (placeholder)
Lp_y  = 4.0e-3;          % m  (placeholder)
Wp_y  = 3.0e-3;          % m  (placeholder)

% Gap / grid geometry (engaged edge length, edge width, inter-cell gap)
g        = 0.50e-3;      % m  (placeholder)
wedge_x  = 6.0e-3;       % m  (placeholder)
wedge_y  = 6.0e-3;       % m  (placeholder)
Ledge_x  = 6.0e-3;       % m  (placeholder)
Ledge_y  = 6.0e-3;       % m  (placeholder)

% Via / loop path (return)
d_via   = 0.60e-3;       % m  (placeholder)
lcur_x  = 6.0e-3;        % m  (placeholder)
wcur_x  = 1.0e-3;        % m  (placeholder)
lcur_y  = 6.0e-3;        % m  (placeholder)
wcur_y  = 1.0e-3;        % m  (placeholder)

% Cross-pol mutual (optional). Set to 0 to neglect.
gd   = 1.0e-3;           % m  (placeholder)
wd   = 2.0e-3;           % m  (placeholder)
Ld   = 6.0e-3;           % m  (placeholder)
useCm = false;            % set false to force zero cross-pol mutual

% Helper: Hammerstad fringing for parasitic patches
hammer_eeff = @(Er,W) (Er+1)/2 + (Er-1)/2 ./ sqrt(1 + 12*h./W);
hammer_dL   = @(eeff,W) 0.412*h .* ((eeff+0.3).*(W/h+0.264)) ./ ((eeff-0.258).*(W/h+0.8));

eeff_x = hammer_eeff(er, Wp_x);
eeff_y = hammer_eeff(er, Wp_y);
dL_x   = hammer_dL(eeff_x, Wp_x);
dL_y   = hammer_dL(eeff_y, Wp_y);

Ap_eff_x = (Lp_x + 2*dL_x) * (Wp_x + 2*dL_x);
Ap_eff_y = (Lp_y + 2*dL_y) * (Wp_y + 2*dL_y);
Cp_x     = eps0*eeff_x * Ap_eff_x / h;
Cp_y     = eps0*eeff_y * Ap_eff_y / h;

% Coplanar-gap per-unit-length capacitance (elliptic integral closed form)
Kfun = @(k) ellipke(k.^2);  % MATLAB's K(m) with m=k^2
k_x   = g./(g + 2*wedge_x);  kp_x = sqrt(1 - k_x.^2);
k_y   = g./(g + 2*wedge_y);  kp_y = sqrt(1 - k_y.^2);
Ccp_x = 4*eps0*eeff_x * (Kfun(k_x)/Kfun(kp_x));   % per-unit-length
Ccp_y = 4*eps0*eeff_y * (Kfun(k_y)/Kfun(kp_y));
Cg_x  = Ccp_x * Ledge_x;
Cg_y  = Ccp_y * Ledge_y;

% Optional mutual (cross-pol) capacitance using same CPW form
if useCm
    kd   = gd./(gd + 2*wd);  kpd = sqrt(1 - kd.^2);
    eeff_d = hammer_eeff(er, wd);
    Ccpd = 4*eps0*eeff_d * (Kfun(kd)/Kfun(kpd));
    Cm   = Ccpd * Ld;
else
    Cm   = 0;
end

%% 2) Symbolic reflectivity matrix R(f) (broadside; frequency symbolic)
syms f positive
w = 2*pi*f;

% Wave impedances at broadside
eta1 = eta0;  eta2 = eta0;

% Varactor junction caps (explicit junction law)
Cvar_x = Cj0 / (1 + Vx/phi)^m;
Cvar_y = Cj0 / (1 + Vy/phi)^m;

% ESR from nominal Q (replace with measured Q(V,f) as needed)
Qx = sym(Qx_nom);
Qy = sym(Qy_nom);
R_ESR_x = 1./(w*Cvar_x*Qx);
R_ESR_y = 1./(w*Cvar_y*Qy);

% Varactor impedances
Zvar_x = R_ESR_x + 1i*w*Ls + 1./(1i*w*Cvar_x);
Zvar_y = R_ESR_y + 1i*w*Ls + 1./(1i*w*Cvar_y);

% Return inductances (via + loop)
Lvia    = mu0*h*(log(4*h/d_via) - 1);
Lloop_x = mu0*lcur_x/pi;
Lloop_y = mu0*lcur_y/pi;
Lret_x  = Lvia + Lloop_x;
Lret_y  = Lvia + Lloop_y;

% Copper skin-effect loop loss
Rs      = sqrt(pi*mu0*w/sigma_cu);
Rloop_x = Rs * (lcur_x/max(wcur_x,1e-12));
Rloop_y = Rs * (lcur_y/max(wcur_y,1e-12));

% Series (through) branch impedances
Zseries_x = Zvar_x + 1i*w*Lret_x + Rloop_x;
Zseries_y = Zvar_y + 1i*w*Lret_y + Rloop_y;

% Shunt admittances (with dielectric loss)
Yp_x = 1i*w*Cp_x * (1 - 1i*tand);
Yp_y = 1i*w*Cp_y * (1 - 1i*tand);
Yg_x = 1i*w*Cg_x * (1 - 1i*(tand/2));
Yg_y = 1i*w*Cg_y * (1 - 1i*(tand/2));

% Per-axis shunt admittances and impedances
Yx = 1./Zseries_x + Yp_x + Yg_x;
Yy = 1./Zseries_y + Yp_y + Yg_y;
Zxx = 1./Yx;
Zyy = 1./Yy;

% Cross term
if Cm==0
    Zxy = sym(0);
else
    Zxy = 1./(1i*w*Cm);
end
Zyx = Zxy;  % reciprocity

% Reflectivity matrix R = (Zs - diag(eta1,eta2)) * (Zs + diag(eta1,eta2))^{-1}
Delta = (Zxx+eta1).*(Zyy+eta2) - Zxy.*Zyx;
Rxx   = ((Zxx-eta1).*(Zyy+eta2) - Zxy.*Zyx) ./ Delta;
Rxy   = (2*eta1*Zxy) ./ Delta;
Ryx   = (2*eta2*Zyx) ./ Delta;
Ryy   = ((Zyy-eta2).*(Zxx+eta1) - Zyx.*Zxy) ./ Delta;
Rmat  = [Rxx, Rxy; Ryx, Ryy];

%% 3) Numeric sweep and eigenvalue plot (3.3–3.7 GHz)
fvec = linspace(3.3e9, 3.7e9, 401);
lam1 = zeros(size(fvec));
lam2 = zeros(size(fvec));

for k = 1:numel(fvec)
    Rk = double(subs(Rmat, f, fvec(k)));
    ev = eig(Rk);
    % Sort by magnitude descending to stabilize color/ordering
    [~, idx] = sort(abs(ev), 'descend');
    ev = ev(idx);
    lam1(k) = ev(1);
    lam2(k) = ev(2);
end

toDB = @(x) 20*log10(abs(x) + eps);
lam1_dB = toDB(lam1);
lam2_dB = toDB(lam2);

figure('Color','w');
plot(fvec/1e9, lam1_dB, 'LineWidth',1.8); hold on;
plot(fvec/1e9, lam2_dB, 'LineWidth',1.8);
grid on; box on;
xlabel('Frequency (GHz)');
ylabel('Eigenvalue magnitude (dB)');
ylim([-0.7 0.1]);
title('Eigenvalues of Reflectivity Matrix, broadside (V_x=10 V, V_y=2.5 V)');
legend('|\lambda_1(\mathbf{R})|','|\lambda_2(\mathbf{R})|','Location','best');
