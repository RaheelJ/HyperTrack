function [solution,status] = Solve(xproblem,xguess,xoptions,Alt,Lat,Lon,Speed,n_wp)
    
    problem=xproblem;
	guess=xguess;
	options=xoptions;
	i=1;
    
    % Modify the problem to add user entered takeoff point 
	problem.states.x0(2:3) = [deg2rad(Lon(i)) deg2rad(Lat(i))];
	problem.states.x0l(2:3) = [deg2rad(Lon(i)) deg2rad(Lat(i))]; 
	problem.states.x0u(2:3) = [deg2rad(Lon(i)) deg2rad(Lat(i))];
    if(n_wp>0)
        problem.states.xfl(1:4) = [Alt(i) deg2rad(Lon(i+1)) deg2rad(Lat(i+1)) Speed(i)]; 
        problem.states.xfu(1:4) = [Alt(i) deg2rad(Lon(i+1)) deg2rad(Lat(i+1)) Speed(i)];
        guess.states(:,1)=[80 Alt(i)];
        guess.states(:,4)=[120 Speed(i)];
    else
        problem.states.xfl(2:3) = [deg2rad(Lon(i+1)) deg2rad(Lat(i+1))]; 
        problem.states.xfu(2:3) = [deg2rad(Lon(i+1)) deg2rad(Lat(i+1))];
    end
	
    % Change the trajectory guess accordingly
	guess.states(:,2)=[deg2rad(Lon(i)) deg2rad(Lon(i+1))];
	guess.states(:,3)=[deg2rad(Lat(i)) deg2rad(Lat(i+1))];

	% Solve the problem for Target trajectory
	[solution(i),MRHistory]=solveMyProblem(problem,guess,options);
	
	% Check if optimal solution found
	if(MRHistory.status ~= 0)
        status = 0;
		return;
    end
	
    status = 1;
    % Save the state variables at the first waypoint
	tf          = solution(1).tf;
    alt0        = speval(solution(1),'X',1,tf);
    lon0        = speval(solution(1),'X',2,tf);
    lat0        = speval(solution(1),'X',3,tf);
    speed0      = speval(solution(1),'X',4,tf);
    path_ang0   = speval(solution(1),'X',5,tf);
    head_ang0   = speval(solution(1),'X',6,tf);
    m0          = 120e3;

    alpha0      = speval(solution(1),'U',1,tf);
    bank_ang0   = speval(solution(1),'U',2,tf);
    u_mf0       = speval(solution(1),'U',3,tf);
	
	for i=2:n_wp
		problem=xproblem;
		guess=xguess;
		options=xoptions;
        
		%Initial Time. t0<tf
		problem.time.t0_min=tf;
		problem.time.t0_max=tf;
		guess.t0=tf;

		% Final time. Let tf_min=tf_max if tf is fixed.
		problem.time.tf_min=tf+100;     
		problem.time.tf_max=tf+4000; 
		guess.tf=tf+2000;
        
        % Modify the problem to add user entered waypoint and starting states 
        problem.states.x0 = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
        problem.states.x0l = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
        problem.states.x0u = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
        problem.states.xfl(1:4) = [Alt(i) deg2rad(Lon(i+1)) deg2rad(Lat(i+1)) Speed(i)]; 
        problem.states.xfu(1:4) = [Alt(i) deg2rad(Lon(i+1)) deg2rad(Lat(i+1)) Speed(i)];
		
        % Change the guess trajectory accordingly
		guess.states(:,1)=[alt0 Alt(i)];
		guess.states(:,2)=[lon0 deg2rad(Lon(i+1))];
		guess.states(:,3)=[lat0 deg2rad(Lat(i+1))];
		guess.states(:,4)=[speed0 Speed(i)];
		guess.states(:,5)=[path_ang0 0];
		guess.states(:,6)=[head_ang0 0];
		guess.states(:,7)=[m0 54431];

        % Add the starting states of inputs
		problem.inputs.u0l = [alpha0 bank_ang0 u_mf0];
		problem.inputs.u0u = [alpha0 bank_ang0 u_mf0];

        % Modify the control input guesses accordingly
		guess.inputs(:,1) = [alpha0     0];
		guess.inputs(:,2) = [bank_ang0  0];
		guess.inputs(:,3) = [u_mf0    163];

        %Solve the problem for Target trajectory
        [solution(i),MRHistory]=solveMyProblem(problem,guess,options);
		
        %Check if optimal solution found
        if(MRHistory.status ~= 0)
            status = 0;
            return;
        end
	
        status = 1;
        % Save the state variables at the next waypoint
		tf          = solution(i).tf;
		alt0        = speval(solution(i),'X',1,tf);
		lon0        = speval(solution(i),'X',2,tf);
		lat0        = speval(solution(i),'X',3,tf);
		speed0      = speval(solution(i),'X',4,tf);
		path_ang0   = speval(solution(i),'X',5,tf);
		head_ang0   = speval(solution(i),'X',6,tf);
		m0          = 120e3;

		alpha0      = speval(solution(i),'U',1,tf);
		bank_ang0   = speval(solution(i),'U',2,tf);
		u_mf0       = speval(solution(i),'U',3,tf);
    end
    
    if (n_wp>0)
        i=n_wp+1;
        problem=xproblem;
        guess=xguess;
        options=xoptions;

        %Initial Time. t0<tf
        problem.time.t0_min=tf;
        problem.time.t0_max=tf;
        guess.t0=tf;

        % Final time. Let tf_min=tf_max if tf is fixed.
        problem.time.tf_min=tf+100;     
        problem.time.tf_max=tf+4000; 
        guess.tf=tf+2000;
        
        % Modify the problem to add user entered landing point and starting states 
        problem.states.x0 = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
        problem.states.x0l = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
        problem.states.x0u = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
        problem.states.xfl(2:3) = [deg2rad(Lon(i+1)) deg2rad(Lat(i+1))]; 
        problem.states.xfu(2:3) = [deg2rad(Lon(i+1)) deg2rad(Lat(i+1))];
        
        % Change the guess trajectory accordingly
        guess.states(:,1)=[alt0 80];
        guess.states(:,2)=[lon0 deg2rad(Lon(i+1))];
        guess.states(:,3)=[lat0 deg2rad(Lat(i+1))];
        guess.states(:,4)=[speed0 120];
        guess.states(:,5)=[path_ang0 0];
        guess.states(:,6)=[head_ang0 0];
        guess.states(:,7)=[m0 54431];
        
        % Add the starting states of inputs
        problem.inputs.u0l = [alpha0 bank_ang0 u_mf0];
        problem.inputs.u0u = [alpha0 bank_ang0 u_mf0];
    
        % Modify the control input guesses accordingly
        guess.inputs(:,1) = [alpha0     0];
        guess.inputs(:,2) = [bank_ang0  0];
        guess.inputs(:,3) = [u_mf0   163];

        %Solve the problem for Target trajectory
        [solution(i),MRHistory]=solveMyProblem(problem,guess,options);

        %Check if optimal solution found
        if(MRHistory.status ~= 0)
            status = 0;
            return;
        end
	    status = 1;
    end
end



