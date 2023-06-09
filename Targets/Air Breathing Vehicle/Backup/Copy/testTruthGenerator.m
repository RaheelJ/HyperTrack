classdef testTruthGenerator  < iTruthGenerator
    %% Public Variables
    properties (Access = public)        
        my_last_update_time;        
        my_project_details; 
        %my_target_states;
        my_num_targets;
        my_targets;
        my_target_start_states;
        my_solution;
        my_status=0;
    end
    
    %% Private Methods
    methods (Access = private)
    end
    
    %% Public Methods
    methods (Access = public)
        
        %Constructor
        function this = testTruthGenerator()
            this = this@iTruthGenerator();
            this.my_last_update_time = -1;
            disp('TruthGenerator constructed');
        end
        
        %Initialization and Trajectory Optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [status, message] = initialize(this, project_details, config)
            this.my_status = 0;
            this.my_project_details = project_details;
            disp(project_details);
                     
            fprintf('config: %s\n', config);
            
            restoredefaultpath;
            mfile_name          = mfilename('fullpath');
            [pathstr,name,ext]  = fileparts(mfile_name);
            cd(pathstr);
            cd("Optimization/Optimal Control Tools/OPTI-master/");
            opti_Install(0, 0, 0);
            cd("../../../")
            
            paths_added = Add_Paths();
            if (paths_added == 0)
                status=0;
                message='Initialization Error! (Could not add required paths)';
                disp(message);
                return;
            end
            
            %Read the config file to extract trajectory endpoints, waypoints and no fly zones
            [Alt0,Lon0,Lat0,Altf,Lonf,Latf,error1]   = Read_Endpoints(config);
            [Alt1,Lon1,Lat1,Speed,n_wp,error2]       = Read_Waypoints(config);
            [Alt2,Lon2,Lat2,Radius,n_nfz,error3]     = Read_NFZ(config);
            if ((error1+error2+error3)>0)
                status=0;
                message='Initialization Error! (Inputs Parameters Violate Limits)';
                disp(message);
                return;
            end
            
            %Initialize values of the target and store as the starting states
            this.my_num_targets = 1;
            this.my_targets = iTruthGenerator.create_target_struct(0);
            this.my_targets.pos_x_or_lat = deg2rad(Lat0);
            this.my_targets.pos_y_or_long = deg2rad(Lon0);
            this.my_targets.pos_in_geo = 1;
            this.my_target_start_states = this.my_targets;
            disp('TruthGenerator initialized');
            
            %Setting up the trajectory endpoints and waypoints for the optimization problem
            Lat=[Lat0];
            Lon=[Lon0];
            for i=1:n_wp
                Lat=[Lat Lat1(i)];
                Lon=[Lon Lon1(i)];
            end
            Lat=[Lat Latf];
            Lon=[Lon Lonf];
            
            % Fetch the problem definition
            [problem,guess]=myProblem;          
            %Add no fly zones to the problem
            if (n_nfz>0)
                if (n_nfz>10)
                    problem.data.auxdata.Lat_NFZ(1:10) = deg2rad(Lat2(1:10));
                    problem.data.auxdata.Lon_NFZ(1:10) = deg2rad(Lon2(1:10));
                    problem.data.auxdata.n_NFZ = 10;
                    problem.constraints.gl(3:10+2) = Radius(1:10);
                else
                    problem.data.auxdata.Lat_NFZ(1:n_nfz) = deg2rad(Lat2);
                    problem.data.auxdata.Lon_NFZ(1:n_nfz) = deg2rad(Lon2);
                    problem.data.auxdata.n_NFZ = n_nfz;
                    problem.constraints.gl(3:n_nfz+2) = Radius;
                end
            end
            % Get options and solver settings (h method)
            options = problem.settings(78);
            %Solve the problem for optimal trajectory and control inputs
            [solution,status] = Solve(problem,guess,options,Alt1,Lat,Lon,Speed,n_wp);
            if(status == 1)
                message = 'Optimal solution found!';
                this.my_solution=solution;
				this.my_status=1;
            else 
                message = 'Optimal solution not found!';
                disp(message);
                return;
            end
            disp(message);
            %assignin('base','Trajectory',solution);
            %Plot_Solution(solution, n_wp+1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Reinitialze to the initial states
        function reinitialize(this)
            this.my_last_update_time = -1;  
            this.my_targets = this.my_target_start_states;
            disp('TruthGenerator reinitialized');
        end
        
        %Get the number of targets
        function [num_targets] = get_expected_num_of_targets(this)
            num_targets = this.my_num_targets;
        end
        
        %Get the initial time at which target takesoff 
        function [start_time] = get_start_time(this)
            disp('Here is the start time: ');
            start_time = 0.0
        end
        
        %Get the current state of the target
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [targets] = get_targets(this, current_time)
                if (this.my_status==1 && this.my_targets.is_alive == 1)
                    this.my_targets.time = current_time;
                    temp = size(this.my_solution);
                    n_seg=temp(2);
                    for i=1:n_seg
                        if (current_time >= this.my_solution(i).t0 && current_time <= this.my_solution(i).tf)
                            this.my_targets.pos_x_or_lat = rad2deg(speval(this.my_solution(i),'X',3,current_time));
                            this.my_targets.pos_y_or_long = rad2deg(speval(this.my_solution(i),'X',2,current_time));
                            this.my_targets.pos_z_or_alt = speval(this.my_solution(i),'X',1,current_time);
                            
                            V = speval(this.my_solution(i),'X',4,current_time);
                            path_ang = speval(this.my_solution(i),'X',5,current_time);
                            head_ang = speval(this.my_solution(i),'X',6,current_time);
                            
                            this.my_targets.vel_x = V*cos(path_ang)*sin(head_ang);
                            this.my_targets.vel_y = V*cos(path_ang)*cos(head_ang);
                            this.my_targets.vel_z = V*sin(path_ang);
                            break;
                        end
                    end
                    if(current_time > this.my_solution(n_seg).tf)
                    this.my_targets.is_alive = 0;
                    end
                end
            %disp('Here is the target state: ');
            targets = this.my_targets;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Make changes to the target
		function apply_changes(this, changes)
			for i=1:length(changes)
				if (changes(i).terminate == 1)
					%terminate();
					this.my_targets.is_alive = 0;
                    this.my_last_update_time = this.my_targets.time;
                elseif (changes(i).change_end_time == 1)
                    % need to modify
                end
			end
		end
		
    end
    
    %% Static methods for coordinate conversions
    methods (Static)
        function [e, n, u] = geo2enu(latitude, longitude, altitude, lrf_lat, lrf_lon, lrf_alt,reference_ellipsoid, angleUnit)
            % default angle unit: degrees
            
            if (nargin < 8)
                angleUnit = 'Degrees';
            end
            
            if (nargin < 7)
                reference_ellipsoid = referenceEllipsoid('wgs84', 'm');
            end
            
            if (strcmpi(angleUnit, 'Degrees'))
                latitude1 = latitude*pi/180;
                longitude1 = longitude*pi/180;
                
            else
                latitude1 = latitude;
                longitude1 = longitude;
            end
            
            [x, y, z] = testTruthGenerator.geo2ecef(latitude1, longitude1, altitude, reference_ellipsoid, 'Radians');
            
            [e, n, u] = testTruthGenerator.ecef2enu(x, y, z, lrf_lat, lrf_lon, lrf_alt, reference_ellipsoid, angleUnit);
            
        end
        
        function [x, y, z] = geo2ecef(phi, lambda, h, ellipsoid, angleUnit)
            % default angle unit: radians
            
            if nargin <= 4 || strcmpi(angleUnit, 'Radians')
                [x, y, z] = testTruthGenerator.geod2ecef(phi*180/pi, lambda*180/pi, h);
            else
                [x, y, z] = testTruthGenerator.geod2ecef(phi, lambda, h);
            end
        end
        
        function [e, n, u] = ecef2enu(X, Y, Z, LAT0, LON0, H0, SPHEROID, angleUnit)
            % default angle unit: degrees
            
            if nargin >= 8 &&  strcmpi(angleUnit, 'Radians')
                lla = [LAT0, LON0, H0]';
            else
                lla = [LAT0*pi/180, LON0*pi/180, H0]';
            end
            Lambda = lla(2);
            Phi = lla(1);
            
            [x0, y0, z0] = testTruthGenerator.geo2ecef(lla(1), lla(2), H0);
            
            XYZ2ENU = [-sin(Lambda)           cos(Lambda)           0;
                -sin(Phi)*cos(Lambda) -sin(Phi)*sin(Lambda)  cos(Phi);
                cos(Phi)*cos(Lambda)  cos(Phi)*sin(Lambda)  sin(Phi)];
            
            enu    = XYZ2ENU*[X-x0,Y-y0,Z-z0]';
            
            e = enu(1);
            n = enu(2);
            u = enu(3);
        end
        
        function [x, y, z] = geod2ecef(latitude, longitude, altitude)
            % default angle unit: degrees
            
            % Input checking/conversion.
            error(nargchk(1, 3, nargin));
            if nargin == 1
                sizelatitude = size(latitude);
                first3 = find(sizelatitude == 3, 1, 'first');
                latitude = reshape(permute(latitude, [first3, 1:(first3 - 1), ...
                    (first3 + 1):ndims(latitude)]), 3, []);
                sizelatitude(first3) = 1;
                longitude = reshape(latitude(2, :), sizelatitude);
                altitude = reshape(latitude(3, :), sizelatitude);
                latitude = reshape(latitude(1, :), sizelatitude);
            end
            latitude = latitude*pi/180; longitude = longitude*pi/180;
            
            % WGS84 parameters.
            a = 6378137; f = 1/298.257223563; b = a*(1 - f); e2 = 1 - (b/a)^2;
            
            % Conversion from:
            % en.wikipedia.org/wiki/Geodetic_system#Conversion_calculations
            Nphi = a ./ sqrt(1 - e2*sin(latitude).^2);
            x = (Nphi + altitude).*cos(latitude).*cos(longitude);
            y = (Nphi + altitude).*cos(latitude).*sin(longitude);
            z = (Nphi.*(1 - e2) + altitude).*sin(latitude);
            
            % Shape output according to number of arguments.
            if nargout <= 1
                if nargin == 1
                    x = cat(first3, x, y, z);
                else
                    dims = ndims(x);
                    if dims == 2
                        if size(x, 2) == 1
                            x = [x, y, z];
                        else
                            x = [x; y; x];
                        end
                    else
                        x = cat(dims + 1, x, y, z);
                    end
                end
            end
        end
    end
end