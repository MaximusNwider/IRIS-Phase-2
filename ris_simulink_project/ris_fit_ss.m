function sysR = ris_fit_ss(R_frd, order)
% Fit MIMO FRD with a stable, minimal state-space model for time simulation.
opts = fitfrdOptions('Weighting','rel', 'Display','on');  % relative error across band
sysR = fitfrd(R_frd, order, opts);    % returns tf/ss; cast to ss for Simulink
sysR = ss(sysR);  sysR = minreal(balreal(sysR), 1e-8);
end
