function [status, message, solution_out, start_states, targets] = Matlab_dll(FID, config_file, solution_in, current_time) %#codegen

    % FID = Function ID for methods of the iTarget class 
    
    status          = 0;
    message         = 'Nil';
    iTarget         = testTruthGenerator();
    start_states    = iTarget.my_target_start_states;
    targets         = iTarget.my_targets;
    solution_out    = [];
    
    switch (FID)
        case 1
            project_details     = iTarget.create_project_detail_struct();
            [status, message]   = iTarget.initialize(project_details, config_file);
            start_states        = iTarget.my_target_start_states;
            targets             = iTarget.my_targets;
            solution_out        = iTarget.my_solution;
            %save('workspace.mat') 
        case 2
            iTarget.my_solution     = solution_in;
            targets                 = iTarget.get_targets(current_time);
        case 3 % For debugging and testing
            start_states(1).a       = 'start_states_1';
            start_states(2).a       = 'start_states_2';
            targets(1).field1       = 99.99;
            targets(1).field2       = 20.02;
            targets(2).field1       = 1;
            targets(2).field2       = 2;
    end
end

