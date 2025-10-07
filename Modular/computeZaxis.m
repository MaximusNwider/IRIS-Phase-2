% File: computeZaxis.m
function Zuu = computeZaxis(fvec, Zvar, M, A)
% M: struct('er','tand','h','sigma'), A: struct('Lp','Wp','Ledge','wedge','lcur','wcur','d_via','g','Cp','Cg')
mu0=4*pi*1e-7; w=2*pi*fvec(:); Zvar=Zvar(:);
% Return path (via+loop) and copper loop loss
Lvia  = mu0*M.h*(log(4*M.h/A.d_via) - 1);
Lloop = mu0*A.lcur/pi;
Lret  = Lvia + Lloop;
Rs    = sqrt(pi*mu0*w/M.sigma);
Rloop = Rs .* (A.lcur/max(A.wcur,1e-12));
% Series branch: measured varactor + board return + copper loss
Zseries = Zvar + 1i*w*Lret + Rloop;
% Shunt admittances (dielectric loss)
Yp = 1i*w*A.Cp .* (1 - 1i*M.tand);
Yg = 1i*w*A.Cg .* (1 - 1i*(M.tand/2));
Ytot = 1./Zseries + Yp + Yg;
Zuu  = 1 ./ Ytot;
Zuu  = Zuu(:).';
end
