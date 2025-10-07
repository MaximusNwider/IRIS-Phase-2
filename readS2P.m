% File: readS2P.m
function [f, S, fmt, fUnit, Z0, info] = readS2P(fname, fallbackZ0)
if nargin<2 || isempty(fallbackZ0), fallbackZ0 = 50; end
fid = fopen(fname,'r'); assert(fid>0); c = onCleanup(@() fclose(fid));
hdrFound=false; comments={}; data=[];
while true
    ln = fgetl(fid); if ~ischar(ln), break; end
    ln = strtrim(ln); if isempty(ln), continue; end
    if strncmp(ln,'!',1), comments{end+1}=ln; continue; end %#ok<AGROW>
    if strncmp(ln,'#',1)
        toks = regexp(upper(strtrim(ln)),'\s+','split');
        assert(numel(toks)>=4,'Malformed header in %s',fname);
        fUnit = toks{2}; assert(strcmpi(toks{3},'S'),'Only S supported');
        fmt   = toks{4}; Z0 = fallbackZ0;
        idxR  = find(strcmpi(toks,'R'),1,'first');
        if ~isempty(idxR) && idxR<numel(toks)
            ztmp = str2double(toks{idxR+1}); if ~isnan(ztmp), Z0 = ztmp; end
        end
        hdrFound=true; continue;
    end
    vals = sscanf(ln,'%f'); if ~isempty(vals), data=[data; vals(:).']; end %#ok<AGROW>
end
assert(hdrFound,'No Touchstone header in %s',fname);
info = struct('comments',{comments});
switch upper(fUnit)
    case 'HZ',  fscale=1; case 'KHZ', fscale=1e3; case 'MHZ', fscale=1e6; case 'GHZ', fscale=1e9;
    otherwise, error('Unsupported frequency unit: %s', fUnit);
end
assert(size(data,2)>=9,'Unexpected column count in %s',fname);
f = data(:,1)*fscale; N=numel(f); S=zeros(2,2,N);
rp = data(:,2:end);
switch upper(fmt)
    case 'RI'
        s11=rp(:,1)+1i*rp(:,2); s21=rp(:,3)+1i*rp(:,4);
        s12=rp(:,5)+1i*rp(:,6); s22=rp(:,7)+1i*rp(:,8);
    case 'MA'
        toC=@(m,a) m.*exp(1i*(pi/180).*a);
        s11=toC(rp(:,1),rp(:,2)); s21=toC(rp(:,3),rp(:,4));
        s12=toC(rp(:,5),rp(:,6)); s22=toC(rp(:,7),rp(:,8));
    case 'DB'
        toC=@(db,a) 10.^(db/20).*exp(1i*(pi/180).*a);
        s11=toC(rp(:,1),rp(:,2)); s21=toC(rp(:,3),rp(:,4));
        s12=toC(rp(:,5),rp(:,6)); s22=toC(rp(:,7),rp(:,8));
    otherwise, error('Unsupported S format: %s',fmt);
end
S(1,1,:)=s11; S(2,1,:)=s21; S(1,2,:)=s12; S(2,2,:)=s22;
end
