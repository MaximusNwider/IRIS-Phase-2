function Z = S2Z_2port(S,Z0)
if size(S,1)==2 && size(S,2)==2
    S2=S;
elseif size(S,2)==2 && size(S,3)==2
    S2=permute(S,[2 3 1]);
else
    error('S2Z_2port:BadShape','S must be 2x2xN or Nx2x2');
end
N=size(S2,3); Z=zeros(2,2,N); I=eye(2);
for k=1:N
    Sk=S2(:,:,k); Z(:,:,k)=Z0*(I+Sk)/(I-Sk);
end
end
