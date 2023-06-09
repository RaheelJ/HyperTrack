classdef iTruthGenerator < handle

    properties (Access = protected)
    end
    
    methods (Access = public, Static)    
        function s = create_project_detail_struct(empty)
            s=struct(...
                'reference_time', 0, ...         % Reference time in seconds from 1970/01/01
                'lrf_latitude', 0, ...           % Local reference frame latitude
                'lrf_longitude', 0, ....         % Local reference frame longitude
                'lrf_altitude', 0);              % Local reference frame altitude                
            if exist('empty','var') && (empty)
                s(1) = [];
            end
        end
        
        function s = create_target_transmitters_struct(empty)
            s=struct(...                
                'ais', 0, ...                  % 0: no AIS transmiiter, 1: has AIS transmitter
                'ads_b', 0, ...                % 0: no ADS-B transmiiter, 1: has ADS-B transmitter
                'esm', 0);                     % 0: no ESM transmiiter, 1: has ESM transmitter
            if exist('empty','var') && (empty)
                s(1) = [];
            end
        end
        function s = create_target_struct(empty)
            s=struct(...
                'id', 0, ...                    % target ID
                'name', '', ...                  % target name
                'type', 'other', ...                  % target type; options: aircraft,ship,boat,satellite,submarine,helicopter,car,truck,person,missile,other
                'is_alive', 1,...               % 1: alive, 0: dead 
                'time', 0, ...                  % time from reference time
                'pos_in_geo', 0, ...            % specifies if the target state in GEO or LRF (local reference frame). "pos_in_geo = 0" means LRF, "pos_in_geo = 1" means GEO
                'pos_x_or_lat', 0, ...          % target x (in m) if "pos_in_geo = 0, otherwise target latitude (in deg)
                'pos_y_or_long', 0, ...         % target y (in m) if "pos_in_geo = 0, otherwise target longitude (in deg)
                'pos_z_or_alt', 0, ...          % target z (in m) if "pos_in_geo = 0, otherwise target altitude (in m)
                'vel_x', 0, ...                 % target velocity in the x/east direction (in m/s)
                'vel_y', 0, ...                 % target velocity in the y/north direction (in m/s)
                'vel_z', 0,...                  % target velocity in the z/up direction (in m/s)
                'transmitters', iTruthGenerator.create_target_transmitters_struct(0));
            if exist('empty','var') && (empty)
                s(1) = [];
            end
        end
        function s = create_object_state_struct(empty)
            s=struct(...
                'time', 0, ...                        % time from reference time
                'pos_in_geo', 0, ...                  % specifies if the object position in GEO or LRF (local reference frame). "pos_in_geo = 0" means LRF, "pos_in_geo = 1" means GEO
                'pos_x_or_lat', 0, ...                % position x (in m) if "pos_in_geo = 0, otherwise object latitude (in deg)
                'pos_y_or_long', 0, ...               % position y (in m) if "pos_in_geo = 0, otherwise object longitude (in deg)
                'pos_z_or_alt', 0, ...                % position z (in m) if "pos_in_geo = 0, otherwise object altitude (in m)
                'vel_x', 0, ...                       % object velocity in the x/east direction (in m/s)
                'vel_y', 0, ...                       % object velocity in the y/north direction (in m/s)
                'vel_z', 0);                          % object velocity in the z/up direction (in m/s)
            if exist('empty','var') && (empty)
                s(1) = [];
            end
        end
        function s = create_waypoint_struct(empty)
            s=struct(...
				'latitude', 0.0, ...
	            'longitude', 0.0, ...
	            'altitude', 0.0, ...
                'rest_time', 0.0, ...
                'leaving_speed', 0.0 ...
			);
			if exist('empty','var') && (empty)
                s(1) = [];
            end
        end
		function s = create_target_change_struct(empty)
            s=struct(...
				'id', 0 , ...
				'internal_target', 0, ...  % 0: target comes fom outside, 1: target created here
	            'name', '' , ...
	            'start', 0, ...
	            'terminate', 0, ...
	            'change_start_time', 0, ...
                'change_end_time', 0, ...
	            'start_time', 0, ...
	            'end_time', 0, ...
	            'change_heading', 0, ...
	            'new_heading', 0.0, ...
	            'new_vertical_heading', 0.0, ...
	            'change_speed', 0, ...
	            'new_speed', 0, ...
	            'change_currrent_state', 0, ...
	            'current_state', iTruthGenerator.create_object_state_struct(), ...
	            'change_waypoints', 0, ...
	            'waypoints', iTruthGenerator.create_waypoint_struct(1) ...			
			);
			if exist('empty','var') && (empty)
                s(1) = [];
            end
		end
    end

    methods (Access = private)        
    end   
    
    methods (Access = public) 
    end
     
    methods (Access=public, Abstract)
        %%% Initialize the truth-generator with given configuration file
        %%% project_details is a struct defined under
        %%% create_project_detail_struct()        
        %%% config file is the project file (xml, m, or any user specified file). It can be empty too.        
        %%% return status(1: success, 0: false), and error_message (if status = 0) 
        [status, message] = initialize(this, project_details, config);

        %%% Reinitialize the truth-generator
        %%% Reinitialize is same as running the project again. So, except parameter setting, all the data may need to be cleared.
        reinitialize(this);

        %%% Get the expected (or exact) number of targets, used for memory allocation
        [num_targets] = get_expected_num_of_targets(this);
		
        %%% Get the start time of the first target
        [start_time] = get_start_time(this);
		
        %%% Get the target states of all the targets at the given time return array of targets as an array of struct defined under create_target_struct()
        [targets] = get_targets(this, current_time);
		%%% Apply the changes
		%%% changes is an array of struct defined under
        %%% create_target_change_struct()        
        apply_changes(changes);
    end

end 