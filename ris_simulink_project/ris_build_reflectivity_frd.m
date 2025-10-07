function [R_frd, Rxx,Rxy,Ryx,Ryy] = ris_build_reflectivity_frd(fvec, file_x_10V, file_y_0V, P, G)
% Returns FRD of 2x2 reflectivity matrix and individual entries over fvec (Hz).

% Load varactor branches strictly from S2P
Zvar_x = loadVaractorZseries(file_x_10V, fvec);
Zvar_y = loadVaractorZseries(file_y_0V,  fvec);

% Parasitic patch & grid caps (Hammerstad + CPW)
eeff_x = hammerEffectivePermittivity(P.er, P.h, G.Wp_x);
eeff_y = hammerEffectivePermittivity(P.er, P.h, G.Wp_y);
dL_x   = hammerDeltaL(eeff_x, P.h, G.Wp_x);
dL_y   = hammerDeltaL(eeff_y, P.h, G.Wp_y);
Ap_eff_x = (G.Lp_x + 2*dL_x) * (G.Wp_x + 2*dL_x);
Ap_eff_y = (G.Lp_y + 2*dL_y) * (G.Wp_y + 2*dL_y);
Cp_x = P.eps0*eeff_x * Ap_eff_x / P.h;
Cp_y = P.eps0*eeff_y * Ap_eff_y / P.h;
Ccp_x = cpwCprime(P.eps0, eeff_x, G.wedge_x, G.g);
Ccp_y = cpwCprime(P.eps0, eeff_y, G.wedge_y, G.g);
Cg_x  = Ccp_x * G.Ledge_x;
Cg_y  = Ccp_y * G.Ledge_y;

% Cross-pol Cm via CPW
eeff_d = hammerEffectivePermittivity(P.er, P.h, G.wd);
k  = G.gd./(G.gd + 2*G.wd);  kp = sqrt(1-k.^2);
K  = ellipke(k.^2); Kp = ellipke(kp.^2);
Cprime_d = 4*P.eps0*eeff_d*(K./Kp);
Cm = Cprime_d * G.Ld;

% Axis impedances (use ONLY measured Zvar)
Zxx = computeZaxis(fvec, Zvar_x, P, struct('Lp',G.Lp_x,'Wp',G.Wp_x,'Ledge',G.Ledge_x, ...
    'wedge',G.wedge_x,'lcur',G.lcur_x,'wcur',G.wcur_x,'d_via',G.d_via,'g',G.g,'Cp',Cp_x,'Cg',Cg_x));
Zyy = computeZaxis(fvec, Zvar_y, P, struct('Lp',G.Lp_y,'Wp',G.Wp_y,'Ledge',G.Ledge_y, ...
    'wedge',G.wedge_y,'lcur',G.lcur_y,'wcur',G.wcur_y,'d_via',G.d_via,'g',G.g,'Cp',Cp_y,'Cg',Cg_y));

w = 2*pi*fvec;
Zxy = 1./(1i*w*Cm);  Zyx = Zxy;           % reciprocity
eta1 = P.eta0; eta2 = P.eta0;

[Rxx,Rxy,Ryx,Ryy] = reflectivityFromZs(fvec, Zxx, Zyy, Zxy, Zyx, eta1, eta2);

% Pack into FRD (2x2xN)
N = numel(fvec);
R = zeros(2,2,N);
for k=1:N, R(:,:,k) = [Rxx(k) Rxy(k); Ryx(k) Ryy(k)]; end
R_frd = frd(R, fvec, 'FrequencyUnit','Hz');   % Control System Toolbox
end
