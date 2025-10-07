function Zvar = loadVaractorZseries(pathS2P, fvec)
[f,S,~,~,Z0,~] = readS2P(pathS2P,50);
Z2p  = S2Z_2port(S,Z0);
Zser = Z2SeriesTwoTerminal(Z2p);                 % series two-terminal
Zvar = interp1(f(:), Zser(:), fvec(:), 'pchip','extrap').';
end
