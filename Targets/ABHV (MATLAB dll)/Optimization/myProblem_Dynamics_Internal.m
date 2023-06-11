
function [dx, g_neq] = myProblem_Dynamics_Internal(x,u,p,t,data)
% Dynamics for internal model 

%------------- BEGIN CODE --------------

    %Stored data
    auxdata         = data.auxdata;         
    
    %Constants
    Re              = auxdata.Re;           %Radius of Earth
    mu              = auxdata.mu;           %Gravitational Parameter
    S               = auxdata.S;            %Reference Area
    omega           = auxdata.omega;        %Earth's Rotation Rate
    g0              = auxdata.g0;           %Gravitational Acceleration
    scale_height    = auxdata.scale_height; %Scale Height
    atm_den_sea     = auxdata.atm_den_sea;  %Atmospheric Density at Sea-Level
    A_ic            = auxdata.A_ic;         %Inlet Capture Area
    phi_st          = auxdata.phi_st;       %Stoichiometric Fuel-Air Ratio       
    Lat_NFZ         = auxdata.Lat_NFZ;
    Lon_NFZ         = auxdata.Lon_NFZ;
    
    %Atmospheric Model
    table_speed     = auxdata.data_speed;
    table_altitude  = auxdata.data_altitude;
    
    %Aerodynamic Model
    table_C_L       = auxdata.data_C_L;          
    table_C_D       = auxdata.data_C_D;
    table_alpha     = auxdata.data_alpha;
    table_mach      = auxdata.data_mach;
    
    %Propulsion Model 
    table_CAR       = auxdata.data_CAR;     %Inlet Capture Area Ratio
    table_throttle  = auxdata.data_throttle;
    table_I_sp      = auxdata.data_I_sp;
    
    %Define states
    alt         = x(:, 1);					%Altitude
    lon         = x(:, 2);					%Longitude
    lat         = x(:, 3);					%Latitude
    speed       = x(:, 4);                  %Speed
    path_ang    = x(:, 5);  				%Flight Path Angle
    head_ang    = x(:, 6);  				%Velocity Heading Angle
    m           = x(:, 7);                  %Mass of Vehicle

    %Define inputs
    alpha       = u(:, 1);                  %Angle of Attack
    bank_ang    = u(:, 2);                  %Bank Angle
    u_mf        = u(:, 3);                  %Mass Flow Rate    
    
    r           = alt+Re;                               %Radial Distance from Centre of the Earth
	alpha_deg   = rad2deg(alpha);                       
    amb_den     = atm_den_sea*exp(-alt./scale_height);  %Ambient Atmospheric Density                     
    q           = power(speed, 2).*amb_den/2;           %Dynamic Pressure
    
    c           = interp1(table_altitude, table_speed, alt, 'spline');  %Speed of Sound
    mach        = speed./c;                                             %Mach Number
    
    %Calculation of Thrust and Specific Impulse
    CAR         = interp2(table_mach, table_alpha, table_CAR, mach, alpha_deg, 'spline');
    throttle    = u_mf./(amb_den.*speed.*CAR*A_ic*phi_st);
    I_sp        = interp2(table_mach, table_throttle, table_I_sp, mach, throttle, 'spline');
    T           = g0*u_mf.*I_sp;                        %Thrust
    
    %Calculation of Lift and Drag
    C_L         = interp2(table_mach, table_alpha, table_C_L, mach, alpha_deg, 'spline');
    C_D         = interp2(table_mach, table_alpha, table_C_D, mach, alpha_deg, 'spline');
	L           = S*q.*C_L;                                                 
	D           = S*q.*C_D;                             
    
    %Define ODE right-hand side
    dx(:, 1)    = speed.*sin(path_ang);
	dx(:, 2)    = speed.*cos(path_ang).*sin(head_ang)./(r.*cos(lat));
	dx(:, 3)    = speed.*cos(path_ang).*cos(head_ang)./r;
	
	term1       = (1./m).*(T.*cos(alpha)-D);
	term2       = mu*sin(path_ang)./power(r, 2);
	term3       = r*power(omega, 2).*cos(lat).*(sin(path_ang).*cos(lat)-cos(path_ang).*sin(lat).*cos(head_ang));
	dx(:, 4)    = term1-term2+term3;
	
	term1       = cos(bank_ang).*(T.*sin(alpha)+L)./(m.*speed);
	term2       = cos(path_ang).*((speed./r)-(mu./(power(r, 2).*speed)));
	term3       = 2*omega*cos(lat).*sin(head_ang);
	term4       = (r./speed).*power(omega, 2).*cos(lat).*(cos(path_ang).*cos(lat)+sin(path_ang).*sin(lat).*cos(head_ang));
	dx(:, 5)    = term1+term2+term3+term4;
	
	term1       = (sin(bank_ang)./(m.*speed.*cos(path_ang))).*(T.*sin(alpha)+L);
	term2       = (speed./r).*cos(path_ang).*sin(head_ang).*tan(lat);
	term3       = 2*omega*(tan(path_ang).*cos(lat).*cos(head_ang)-sin(lat));
	term4       = (r.*power(omega, 2)./(speed.*cos(path_ang))).*sin(lat).*cos(lat).*sin(head_ang);
	dx(:, 6)    = term1+term2-term3+term4;
	
	dx(:, 7)    = -T./(I_sp*g0);
    
    %Define Path constraints
    g_neq(:, 1)  = q;           % Constraint on dynamic pressure
    g_neq(:, 2)  = throttle;    % Constraint on maximum throttle
    
    % No fly zone constraints
    g_neq(:, 3)  = sqrt((lon-Lon_NFZ(1)).^2 + (lat-Lat_NFZ(1)).^2);
    g_neq(:, 4)  = sqrt((lon-Lon_NFZ(2)).^2 + (lat-Lat_NFZ(2)).^2);
    g_neq(:, 5)  = sqrt((lon-Lon_NFZ(3)).^2 + (lat-Lat_NFZ(3)).^2);
    g_neq(:, 6)  = sqrt((lon-Lon_NFZ(4)).^2 + (lat-Lat_NFZ(4)).^2);
    g_neq(:, 7)  = sqrt((lon-Lon_NFZ(5)).^2 + (lat-Lat_NFZ(5)).^2);
    g_neq(:, 8)  = sqrt((lon-Lon_NFZ(6)).^2 + (lat-Lat_NFZ(6)).^2);
    g_neq(:, 9)  = sqrt((lon-Lon_NFZ(7)).^2 + (lat-Lat_NFZ(7)).^2);
    g_neq(:, 10) = sqrt((lon-Lon_NFZ(8)).^2 + (lat-Lat_NFZ(8)).^2);
    g_neq(:, 11) = sqrt((lon-Lon_NFZ(9)).^2 + (lat-Lat_NFZ(9)).^2);
    g_neq(:, 12) = sqrt((lon-Lon_NFZ(10)).^2 + (lat-Lat_NFZ(10)).^2);
    
   
%------------- END OF CODE --------------