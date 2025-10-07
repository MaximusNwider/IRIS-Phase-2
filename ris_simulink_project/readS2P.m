function [f, S, fmt, fUnit, Z0, info] = readS2P(fname, fallbackZ0)
if nargin<2, fallbackZ0=50; end
fid=fopen(fname,'r'); assert(fid>0); c=onCleanup(@() fclose(fid));
hdrFound=false; comments={}; data=[];
while true
    ln=fgetl(fid); if ~ischar(ln), break; end; ln=strtrim(ln); if isempty(ln), continue; end
    if ln(1)=='!', comments{end+1}=ln; continue; end %#ok<AGROW>
    if ln(1)=='#'
        toks=regexp(upper(ln),'\s+','split'); fUnit=toks{2}; fmt=toks{4}; Z0=fallbackZ0;
        idx=find(strcmpi(toks,'R'),1); if ~isempty(idx)&&idx<numel(toks), z=str2double(toks{idx+1}); if ~isnan(z), Z0=z; end, end
        hdrFound=true; continue
    end
    v=sscanf(ln,'%f'); if ~isempty(v), data=[data; v(:).']; end %#ok<AGROW>
end
assert(hdrFound,'No Touchstone header'); info=struct('comments',{comments});
switch upper(fUnit), case 'HZ', fs=1; case 'KHZ', fs=1e3; case 'MHZ', fs=1e6; case 'GHZ', fs=1e9; otherwise, error('bad unit'); end
f=data(:,1)*fs; %#ok<*NBRAK>
rp=data(:,2:end); S=zeros(2,2,numel(f));
switch upper(fmt)
    case 'RI', s11=rp(:,1)+1i*rp(:,2); s21=rp(:,3)+1i*rp(:,4); s12=rp(:,5)+1i*rp(:,6); s22=rp(:,7)+1i*rp(:,8);
    case 'MA', toC=@(m,a)m.*exp(1i*pi/180*a); s11=toC(rp(:,1),rp(:,2)); s21=toC(rp(:,3),rp(:,4)); s12=toC(rp(:,5),rp(:,6)); s22=toC(rp(:,7),rp(:,8));
    case 'DB', toC=@(d,a)10.^(d/20).*exp(1i*pi/180*a); s11=toC(rp(:,1),rp(:,2)); s21=toC(rp(:,3),rp(:,4)); s12=toC(rp(:,5),rp(:,6)); s22=toC(rp(:,7),rp(:,8));
    otherwise, error('fmt');
end
S(1,1,:)=s11; S(2,1,:)=s21; S(1,2,:)=s12; S(2,2,:)=s22;
end
