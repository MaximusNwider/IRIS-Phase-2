% File: computeCoCross.m
function [co_db, xpol_db, co_phase_deg, xpol_phase_deg] = ...
    computeCoCross(fvec, Rxx, Rxy, Ryx, Ryy, mode)
% mode: 'avg' (default), 'xinc', 'yinc', or 'max'
if nargin<6 || isempty(mode), mode='avg'; end
Rxx=Rxx(:); Rxy=Rxy(:); Ryx=Ryx(:); Ryy=Ryy(:); fvec=fvec(:);

switch lower(mode)
    case 'xinc'   % x-polarized incidence
        co   = Rxx;  xpol = Ryx;
    case 'yinc'   % y-polarized incidence
        co   = Ryy;  xpol = Rxy;
    case 'max'    % worst case across axes
        % pick element with larger magnitude at each f
        pickCo   = abs(Rxx)>=abs(Ryy);
        pickXpol = abs(Rxy)>=abs(Ryx);
        co   = Rxx;  co(~pickCo)   = Ryy(~pickCo);
        xpol = Rxy;  xpol(~pickXpol)= Ryx(~pickXpol);
    otherwise     % 'avg' — equal mix of x/y excitations
        % magnitude via RMS power, phase via phasor average
        co_mag   = sqrt(0.5*(abs(Rxx).^2 + abs(Ryy).^2));
        xpol_mag = sqrt(0.5*(abs(Rxy).^2 + abs(Ryx).^2));
        co_ph    = angle(0.5*(Rxx + Ryy));
        xpol_ph  = angle(0.5*(Rxy + Ryx));
        co_db    = 20*log10(max(co_mag, realmin('double')));
        xpol_db  = 20*log10(max(xpol_mag, realmin('double')));
        co_phase_deg   = rad2deg(unwrap(co_ph));
        xpol_phase_deg = rad2deg(unwrap(xpol_ph));
        return;
end

% dB and phase for the chosen pair
co_db    = 20*log10(max(abs(co),   realmin('double')));
xpol_db  = 20*log10(max(abs(xpol), realmin('double')));
co_phase_deg   = rad2deg(unwrap(angle(co)));
xpol_phase_deg = rad2deg(unwrap(angle(xpol)));
end
