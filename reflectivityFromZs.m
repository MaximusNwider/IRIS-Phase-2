% File: reflectivityFromZs.m
function [Rxx,Rxy,Ryx,Ryy] = reflectivityFromZs(fvec, Zxx, Zyy, Zxy, Zyx, eta1, eta2)
Delta = (Zxx+eta1).*(Zyy+eta2) - Zxy.*Zyx;
Rxx = ((Zxx-eta1).*(Zyy+eta2) - Zxy.*Zyx) ./ Delta;
Rxy = (2*eta1).*Zxy ./ Delta;
Ryx = (2*eta2).*Zyx ./ Delta;
Ryy = ((Zyy-eta2).*(Zxx+eta1) - Zyx.*Zxy) ./ Delta;
end
