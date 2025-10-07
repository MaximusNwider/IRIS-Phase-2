% File: Z2SeriesTwoTerminal.m
function Zser = Z2SeriesTwoTerminal(Z2p)

if(size(Z2p,1)==2 && size(Z2p,2)==2)
    Z=Z2p;


elseif(size(Z2p,2)==2 && size(Z2p,3)==2)
    Z=permute(Z2p,[2 3 1]);
 
else
    msg = 'Error. \BadShape';
    error(msg);
end

Zser = squeeze(Z(1,1,:)+Z(2,2,:)-Z(1,2,:)-Z(2,1,:));
end
