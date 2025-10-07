% File: loadVaractorZseries.m
function Zvar = loadVaractorZseries(pathS2P, fvec)
[f,S,~,~,Z0,~] = readS2P(pathS2P,50);
Z2p  = S2Z_2port(S,Z0);
Zser = Z2SeriesTwoTerminal(Z2p);            % measured series branch
Zvar = interp1(f(:), Zser(:), fvec(:), 'pchip', 'extrap');
Zvar = Zvar(:).';                            % row vector
end
