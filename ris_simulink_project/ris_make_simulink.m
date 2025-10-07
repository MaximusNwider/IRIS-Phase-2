function ris_make_simulink(sysR, mdl)
% Creates: two inputs -> MIMO State-Space -> two outputs
if bdIsLoaded(mdl), close_system(mdl,0); end
new_system(mdl); open_system(mdl);

% Put matrices in base workspace with unique names
assignin('base',[mdl '_A'],sysR.A);
assignin('base',[mdl '_B'],sysR.B);
assignin('base',[mdl '_C'],sysR.C);
assignin('base',[mdl '_D'],sysR.D);

% Blocks
add_block('simulink/Sources/In1',[mdl '/Ex_in'],'Position',[40 60 70 80]);
add_block('simulink/Sources/In1',[mdl '/Ey_in'],'Position',[40 130 70 150]);
add_block('simulink/Signal Routing/Mux',[mdl '/Mux'],'Inputs','2','Position',[110 60 140 150]);

add_block('simulink/Continuous/State-Space',[mdl '/R_ss'], ...
    'A',[mdl '_A'],'B',[mdl '_B'],'C',[mdl '_C'],'D',[mdl '_D'], ...
    'Position',[200 50 360 160]);

add_block('simulink/Signal Routing/Demux',[mdl '/Demux'], ...
    'Outputs','2','Position',[410 60 440 150]);

add_block('simulink/Sinks/Out1',[mdl '/Ex_ref'],'Position',[480 70 510 90]);
add_block('simulink/Sinks/Out1',[mdl '/Ey_ref'],'Position',[480 130 510 150]);

% Wires
add_line(mdl,'Ex_in/1','Mux/1');  add_line(mdl,'Ey_in/1','Mux/2');
add_line(mdl,'Mux/1','R_ss/1');   add_line(mdl,'R_ss/1','Demux/1');
add_line(mdl,'Demux/1','Ex_ref/1'); add_line(mdl,'Demux/2','Ey_ref/1');

% Annotation
add_block('simulink/Notes/Note',[mdl '/Note'], ...
    'Position',[200 10 520 40], 'FontSize','12', 'UserDataPersistent','on');
set_param([mdl '/Note'],'Text',sprintf('RIS Reflectivity R(s) — 2x2 MIMO LTI\nInputs: [Ex; Ey]  Outputs: [Ex_ref; Ey_ref]'));

save_system(mdl);
end
