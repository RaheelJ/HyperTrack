% Target = testTruthGenerator();
% Project_Details = Target.create_project_detail_struct();
% Config = 'Configuration/Output_ABH.xml';
% Target.initialize(Project_Details, Config);

FID             = 3;
config_file     = 'Output_ABH.xml';
solution_in     = [];
current_time    = 0;

[status, message, solution_out, start_states, targets] = Matlab_dll(FID, config_file, solution_in, current_time);

% FID             = 3;
% iTarget_in      = iTarget_out;
% config_file     = 'Configuration/Output_ABH.xml';
% current_time    = 0;
% 
% [iTarget_out, start_time, targets, num_targets, status, message] = Matlab_dll(iTarget_in, FID, config_file, current_time);