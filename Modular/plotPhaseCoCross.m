% File: plotPhaseCoCross.m
function plotPhaseCoCross(fvec, co_deg, xpol_deg, ttl, xlim_GHz, ylim_deg)
figure('Color','w');
plot(fvec/1e9, co_deg,   'LineWidth',1.8); hold on;
plot(fvec/1e9, xpol_deg, 'LineWidth',1.8, 'LineStyle','--');
grid on; box on;
xlabel('Frequency (GHz)'); ylabel('Phase (degrees)');
if nargin>=4 && ~isempty(ttl), title(ttl); end
if nargin>=5, xlim(xlim_GHz); end
if nargin>=6, ylim(ylim_deg);  end
legend('Copolar (solid)','Crosspolar (dashed)','Location','best');
end


%[co_db,xpol_db,co_ph,xpol_ph] = computeCoCross(fvec,Rxx,Rxy,Ryx,Ryy,'xinc'); % or 'yinc'
