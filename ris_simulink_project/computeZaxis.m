function Zuu = computeZaxis(fvec, Zvar, P, A)
mu0 = 4*pi*1e-7; w = 2*pi*fvec(:); Zvar = Zvar(:);
Lvia  = mu0*P.h*(log(4*P.h/A.d_via)-1);
Lloop = mu0*A.lcur/pi;  Lret = Lvia + Lloop;
Rs    = sqrt(pi*mu0*w/P.sigma);
Rloop = Rs .* (A.lcur/max(A.wcur,1e-12));
Zseries = Zvar + 1i*w*Lret + Rloop;
Yp = 1i*w*A.Cp .* (1 - 1i*P.tand);
Yg = 1i*w*A.Cg .* (1 - 1i*(P.tand/2));
Ytot = 1./Zseries + Yp + Yg;
Zuu = (1 ./ Ytot).';
end
