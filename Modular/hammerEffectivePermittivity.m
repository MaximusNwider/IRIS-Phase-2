% File: hammerEffectivePermittivity.m
function eeff = hammerEffectivePermittivity(er, h, W)
eeff = (er+1)/2 + (er-1)/2 ./ sqrt(1 + 12*h./W);
end
