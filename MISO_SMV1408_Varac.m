%% ========================================================================
%  build_varactor_model_from_s2p.m  (Script with local functions)
%  - Reads a directory of .s2p files (one per bias)
%  - Converts S -> Z (two-port) -> Z_series (two-terminal)
%  - Aligns onto a common frequency grid (overlap range)
%  - Builds a callable model.fnZ(f,V) to return Z_series via interpolation
%  - Optionally decomposes Z into ESR(f,V), Ceq(f,V), and Ls_est(V)
%  - Saves model to MAT (optional)
%  ========================================================================

%% USER INPUTS
dataDir   = "path/to/your/s2p/folder";     % <<<<<< set me
pattern   = "*.s2p";                       % e.g., 'SMV1408_bias_*.s2p'
saveMat   = false;                         % true to save model to MAT file
matOut    = "varactor_model_from_s2p.mat"; % MAT filename if saveMat=true

% (Optional) Known Z0 if file header lacks it (usually 50 ohm)
fallbackZ0 = 50;

% (Optional) If your filenames don't contain bias, provide a map here:
% manualBias = containers.Map( ...
%    ["fileA.s2p","fileB.s2p"], [10.0, 2.5]);  % volts
manualBias = containers.Map('KeyType','char','ValueType','double');

% (Optional) Robust fitting of Ls & C (set true to estimate Ls per bias)
doFitRLC = true;   % recommended: true

%% GATHER FILES
files = dir(fullfile(dataDir, pattern));
assert(~isempty(files), 'No .s2p files found with pattern "%s" in %s', pattern, dataDir);

allBias = [];
raw = struct('fname',{},'bias',{},'f',{},'S',{},'Z0',{});

for k = 1:numel(files)
    fname = fullfile(files(k).folder, files(k).name);
    [f, S, fmt, fUnit, Z0, info] = readS2P(fname, fallbackZ0);
    Vbias = detectBias(files(k).name, info, manualBias);  % try filename first, then comments, then map
    fprintf('Loaded %-40s  Bias=%g V  Z0=%g Ohm  Format=%s  Unit=%s\n', ...
        files(k).name, Vbias, Z0, fmt, fUnit);

    % Convert S->Z, then two-terminal series impedance
    Z2p = S2Z(S, Z0);                          % [2x2xN] Z-parameters
    Zseries = squeeze(Z2p(1,1,:) + Z2p(2,2,:) - Z2p(1,2,:) - Z2p(2,1,:));  % Nx1

    raw(end+1) = struct('fname',files(k).name, 'bias',Vbias, ...
                        'f',f(:), 'S',S, 'Z0',Z0, 'Zseries',Zseries(:)); %#ok<SAGROW>
    allBias(end+1) = Vbias; %#ok<SAGROW>
end

% Sort by bias ascending
[allBias, order] = sort(allBias);
raw = raw(order);

%% BUILD COMMON FREQUENCY GRID (overlap range)
fmins = arrayfun(@(r) min(r.f), raw);
fmaxs = arrayfun(@(r) max(r.f), raw);
fmin  = max(fmins);
fmax  = min(fmaxs);
assert(fmax>fmin, 'Frequency ranges across files do not overlap.');

% Choose a grid density ~ min native resolution
df_list = arrayfun(@(r) median(diff(unique(r.f))), raw, 'UniformOutput', true);
df      = max(min(df_list), 1); % at least 1 Hz to avoid degenerate step
fGrid   = (ceil(fmin/df)*df : df : floor(fmax/df)*df).';

%% INTERPOLATE EACH BIAS ONTO COMMON GRID
nF = numel(fGrid);
nV = numel(allBias);
Zgrid = complex(nF, nV);

for j = 1:nV
    fj  = raw(j).f;
    Zj  = raw(j).Zseries;
    % Interpolate real & imag separately (monotonic f)
    Zr  = interp1(fj, real(Zj), fGrid, 'pchip', 'extrap');
    Zi  = interp1(fj, imag(Zj), fGrid, 'pchip', 'extrap');
    Zgrid(:,j) = Zr + 1i*Zi;
end

%% OPTIONAL R-L-C DECOMPOSITION (per bias)
% Im(Z) = w*Ls - 1/(w*C). Multiply by w:  y = w*Im(Z) = (w^2)*Ls - 1/C  (linear in x=w^2)
% Fit y vs x to get Ls (slope) and C (from intercept).
ESR = real(Zgrid);
Ls_est = zeros(1,nV);
Ceq = zeros(nF,nV);
if doFitRLC
    w = 2*pi*fGrid;
    x = (w.^2);
    for j = 1:nV
        y = w .* imag(Zgrid(:,j));
        % Robust linear fit y = a*x + b  (a=Ls, b=-1/C)
        p = robustfit(x, y);   % p(1)=intercept b, p(2)=slope a
        a = p(2);  b = p(1);
        Ls_est(j) = max(real(a), 0);           % enforce non-negative
        Cest = -1./max(real(b), eps);          % nominal C from intercept
        % Frequency-dependent C using fitted Ls_est and measured X:
        X   = imag(Zgrid(:,j));
        Ceq(:,j) = 1 ./ (w .* max(w*Ls_est(j) - X, 1e-30)); % keep positive denom
        % Where Ceq is unreasonable, fall back to Cest
        bad = ~isfinite(Ceq(:,j)) | (Ceq(:,j)<=0);
        Ceq(bad,j) = Cest;
    end
else
    Ls_est(:) = 0;
    Ceq(:)    = NaN;
end

%% BUILD INTERPOLANTS & MODEL HANDLE
% Interpolants over (f,V) for Re(Z), Im(Z)
Fr = griddedInterpolant({fGrid, allBias}, real(Zgrid), 'pchip', 'nearest');
Fi = griddedInterpolant({fGrid, allBias}, imag(Zgrid), 'pchip', 'nearest');

model = struct();
model.f         = fGrid;            % Hz
model.V         = allBias;          % V
model.Zseries   = Zgrid;            % [nF x nV] complex
model.ESR       = ESR;              % [nF x nV]
model.Ceq       = Ceq;              % [nF x nV]
model.Ls_est    = Ls_est;           % [1 x nV]
model.Z0        = raw(1).Z0;        % assume common
model.files     = string({raw.fname});
model.fnZ       = @(f,V) interpZ_series(f, V, Fr, Fi);       % complex Z(f,V)
model.fnRCLeq   = @(f,V) rcl_from_model(f, V, model);        % optional: ESR, Ceq, Ls_est

if saveMat
    save(matOut, 'model', '-v7.3');
    fprintf('Saved model -> %s\n', matOut);
end

%% (Optional) sanity plot of |Z| vs f for a couple of biases
%{
figure('Color','w'); hold on; grid on; box on;
fplot = linspace(min(model.f), max(model.f), 400);
bsel = [model.V(1), model.V(end)];
for bv = bsel
    Zp = model.fnZ(fplot, bv);
    plot(fplot/1e9, abs(Zp), 'DisplayName', sprintf('|Z| @ V=%.3g V', bv), 'LineWidth',1.6);
end
xlabel('Frequency (GHz)'); ylabel('|Z_{series}| (\Omega)'); legend show;
%}

%% ======================= LOCAL FUNCTIONS =========================
function [f, S, fmt, fUnit, Z0, info] = readS2P(fname, fallbackZ0)
    % Robust .s2p reader (2-port Touchstone); supports RI/MA/DB formats.
    % Returns:
    %   f    : Nx1 frequency (Hz)
    %   S    : 2x2xN complex S-matrix
    %   fmt  : 'RI'|'MA'|'DB'
    %   fUnit: 'HZ'|'KHZ'|'MHZ'|'GHZ'
    %   Z0   : reference impedance (Ohm)
    %   info : struct with header/comments (for bias parsing)
    arguments
        fname (1,:) char
        fallbackZ0 (1,1) double = 50
    end
    fid = fopen(fname,'r');
    assert(fid>0, 'Cannot open %s', fname);
    c = onCleanup(@() fclose(fid));

    hdrFound = false;
    comments = strings(0,1);
    data = [];
    while true
        ln = fgetl(fid);
        if ~ischar(ln); break; end
        ln = strtrim(ln);
        if isempty(ln); continue; end
        if startsWith(ln,'!') || startsWith(ln,'#') && ~hdrFound
            if startsWith(ln,'!')
                comments(end+1) = string(ln); %#ok<AGROW>
            elseif startsWith(ln,'#')
                % Header example: "# GHZ S RI R 50"
                toks = split(upper(strtrim(ln)));
                % tokens: #  FUNIT  S  FORMAT  R  Z0
                assert(numel(toks)>=5, 'Malformed header in %s', fname);
                fUnit = toks(2);
                assert(toks(3)=="S", 'Only S-parameters supported');
                fmt   = toks(4);   % 'RI'|'MA'|'DB'
                Z0    = fallbackZ0;
                idxR  = find(toks=="R", 1);
                if ~isempty(idxR) && idxR< numel(toks)
                    Z0 = str2double(toks(idxR+1));
                    if isnan(Z0); Z0 = fallbackZ0; end
                end
                hdrFound = true;
            end
            continue;
        end
        if startsWith(ln,'!'); comments(end+1) = string(ln); continue; end %#ok<AGROW>
        % Data line: f  S11  S21  S12  S22 (pairs)
        vals = sscanf(ln, '%f');
        data = [data; vals(:).']; %#ok<AGROW>
    end
    assert(hdrFound, 'No Touchstone header (# ...) found in %s', fname);
    info = struct('comments',{comments});

    % Convert frequency unit to Hz
    switch fUnit
        case "HZ",   fscale = 1;
        case "KHZ",  fscale = 1e3;
        case "MHZ",  fscale = 1e6;
        case "GHZ",  fscale = 1e9;
        otherwise, error('Unsupported frequency unit: %s', fUnit);
    end

    % Parse data into S
    % Touchstone 2-port: cols = 1 + 8 numbers (4 complex pairs)
    % Order: f, (S11 pair), (S21 pair), (S12 pair), (S22 pair)
    assert(size(data,2) >= 9, 'Unexpected .s2p column count in %s', fname);
    f  = data(:,1) * fscale;
    N  = numel(f);
    S  = zeros(2,2,N);
    rawPairs = data(:,2:end);

    switch fmt
        case "RI"
            % real/imag pairs
            s11 = rawPairs(:,1)+1i*rawPairs(:,2);
            s21 = rawPairs(:,3)+1i*rawPairs(:,4);
            s12 = rawPairs(:,5)+1i*rawPairs(:,6);
            s22 = rawPairs(:,7)+1i*rawPairs(:,8);
        case "MA"
            % magnitude/angle(deg)
            toC = @(m,a) m .* exp(1i*(pi/180).*a);
            s11 = toC(rawPairs(:,1), rawPairs(:,2));
            s21 = toC(rawPairs(:,3), rawPairs(:,4));
            s12 = toC(rawPairs(:,5), rawPairs(:,6));
            s22 = toC(rawPairs(:,7), rawPairs(:,8));
        case "DB"
            % dB/angle(deg)
            toC = @(db,a) 10.^(db/20) .* exp(1i*(pi/180).*a);
            s11 = toC(rawPairs(:,1), rawPairs(:,2));
            s21 = toC(rawPairs(:,3), rawPairs(:,4));
            s12 = toC(rawPairs(:,5), rawPairs(:,6));
            s22 = toC(rawPairs(:,7), rawPairs(:,8));
        otherwise
            error('Unsupported S-parameter format: %s', fmt);
    end

    S(1,1,:) = s11;
    S(2,1,:) = s21;
    S(1,2,:) = s12;
    S(2,2,:) = s22;
end

function Vbias = detectBias(filename, info, manualBias)
    % Try filename like *_10V.s2p, *_2.5v.s2p, *V=10.s2p, *Bias10V*, etc.
    Vbias = NaN;
    fn = char(filename);
    % 1) Common patterns
    pat = '(?<num>\d+(\.\d+)?)\s*(V|v)\b';
    m = regexp(fn, pat, 'names');
    if ~isempty(m)
        Vbias = str2double(m(1).num);
    end
    % 2) Look into comments for 'Bias' or 'V='
    if isnan(Vbias) && ~isempty(info) && isfield(info,'comments')
        for c = info.comments(:).'
            ln = char(c);
            m = regexp(ln, '(Bias|V)\s*[:=]\s*(?<num>\d+(\.\d+)?)', 'names', 'ignorecase');
            if ~isempty(m); Vbias = str2double(m(1).num); break; end
        end
    end
    % 3) Manual map fallback
    if isnan(Vbias) && isKey(manualBias, fn)
        Vbias = manualBias(fn);
    end
    assert(~isnan(Vbias), 'Could not determine bias for file: %s', filename);
end

function Z = S2Z(S, Z0)
    % Convert 2-port S to Z:  Z = Z0 * (I+S) * inv(I-S)
    % S: 2x2xN
    N = size(S,3);
    Z = zeros(2,2,N);
    I = eye(2);
    for k = 1:N
        Sk = S(:,:,k);
        Z(:,:,k) = Z0 * (I + Sk) / (I - Sk);
    end
end

function Z = interpZ_series(fq, V, Fr, Fi)
    % Complex interpolation of Z_series over (f,V)
    fq = fq(:);
    V  = V(:).';
    [FF,VV] = ndgrid(fq, V);
    Zr = Fr(FF, VV);
    Zi = Fi(FF, VV);
    Z  = Zr + 1i*Zi;
    if isvector(fq) && isscalar(V)
        Z = Z(:); % N x 1
    end
end

function rcl = rcl_from_model(fq, V, model)
    % Return struct with ESR, Ceq, and Ls_est interpolated at (f,V)
    fq = fq(:);
    V  = V(:).';
    % GriddedInterpolant for ESR and Ceq (may contain NaNs if not fitted)
    FrE = griddedInterpolant({model.f, model.V}, model.ESR, 'pchip', 'nearest');
    FrC = griddedInterpolant({model.f, model.V}, model.Ceq, 'pchip', 'nearest');
    ESR = FrE(fq, V);
    Ceq = FrC(fq, V);
    % Ls_est is per-bias scalar; interp across bias
    LsI = griddedInterpolant(model.V, model.Ls_est, 'pchip', 'nearest');
    Ls  = LsI(V) * ones(numel(fq), numel(V));
    rcl = struct('ESR',ESR, 'Ceq',Ceq, 'Ls',Ls);
end
