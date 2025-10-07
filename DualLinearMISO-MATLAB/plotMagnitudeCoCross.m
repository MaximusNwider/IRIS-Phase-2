% File: plotMagnitudeCoCross.m
function plotMagnitudeCoCross(fvec, co_db, xpol_db, ttl, xlim_GHz, ylim_dB)
figure('Color','w');
plot(fvec/1e9, co_db,   'LineWidth',1.8); hold on;
plot(fvec/1e9, xpol_db, 'LineWidth',1.8, 'LineStyle','--');
grid on; box on;
xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)');
if nargin>=4 && ~isempty(ttl), title(ttl); end
if nargin>=5, xlim(xlim_GHz); end
if nargin>=6, ylim(ylim_dB);  end
legend('Copolar (solid)','Crosspolar (dashed)','Location','best');
end
