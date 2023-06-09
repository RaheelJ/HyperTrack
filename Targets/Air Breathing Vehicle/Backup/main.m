clear
clc
close all

% Target = testTruthGenerator();
% Project_Details = Target.create_project_detail_struct();
% Config = 'Configuration/Output_ABH.xml';
% Target.initialize(Project_Details, Config);


truthgen = testTruthGenerator();
project_detail = iTruthGenerator.create_project_detail_struct();      
truthgen.initialize(project_detail, 'Configuration/AB_Hypersonic_Aircraft-1_d1_0.xml');

num_targets = truthgen.get_expected_num_of_targets();
targets = truthgen.get_targets(0);
start_time= truthgen.get_start_time();

changes = iTruthGenerator.create_target_change_struct();      
truthgen.apply_changes(changes);

