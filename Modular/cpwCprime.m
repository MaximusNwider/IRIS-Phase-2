% File: cpwCprime.m
function Cprime = cpwCprime(eps0, eeff, wedge, g)
k  = g./(g + 2*wedge); kp = sqrt(1 - k.^2);
K  = ellipke(k.^2);    Kp = ellipke(kp.^2);
Cprime = 4*eps0*eeff * (K./Kp);
end