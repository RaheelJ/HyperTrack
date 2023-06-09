function [problem,guess] = myProblem
%File for optimal control problem definition

%------------- BEGIN CODE --------------

% Plant model name, provide in the format of function handle
InternalDynamics=@myProblem_Dynamics_Internal; 
SimDynamics=@myProblem_Dynamics_Sim;

% Analytic derivative files (optional), provide in the format of function handle
problem.analyticDeriv.gradCost=[];
problem.analyticDeriv.hessianLagrangian=[];
problem.analyticDeriv.jacConst=[];

% Settings file
problem.settings=@settings_myProblem;

%Initial Time. t0<tf
problem.time.t0_min=0;
problem.time.t0_max=0;
guess.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
problem.time.tf_min=100;     
problem.time.tf_max=4000; 
guess.tf=2000;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Constants
auxdata.Re              = 6378.166*1000;            %Radius of Earth
auxdata.mu              = 398600.4405e9;            %Gravitational Parameter
auxdata.S               = 557;                      %Reference Area
auxdata.omega           = 0.00007292115856;         %Earth's Rotation Rate
auxdata.g0              = 9.8066498;                %Gravitational Acceleration
auxdata.scale_height    = 10.4e3;                   %Scale Height
auxdata.atm_den_sea     = 1.225;                    %Atmospheric Density at Sea-Level
auxdata.A_ic            = 28;                       %Inlet Capture Area
auxdata.phi_st          = 0.0292;                   %Stoichiometric Fuel-Air Ratio
auxdata.ee              = 0.9999;                   %Convex Combination Coefficient for Objective
auxdata.Lon_NFZ 		= [0 0 0 0 0 0 0 0 0 0];
auxdata.Lat_NFZ 		= [0 0 0 0 0 0 0 0 0 0];
auxdata.n_NFZ           = 0;

% Look-up Tables

% Speed vs Altitude
auxdata.data_speed		= [1116,	968,	968,	995,	1082,	1082,	964,	899]*0.3048;
auxdata.data_altitude 	= [0,		36152,	65824,	105520,	155349,	168677,	235571,	282152]*0.3048;

% Lift Coefficient 
auxdata.data_C_L =    [-0.1330,    -0.1357,	   -0.1414,	   -0.1384,	   -0.1430,	   -0.1399,	   -0.1459,	   -0.1301,	   -0.1050,	   -0.0781,	   -0.0584,	   -0.0471,	   -0.0392;
                       -0.0418,	   -0.0424,	   -0.0454,	   -0.0412,	   -0.0371,	   -0.0313,	   -0.0409,	   -0.0401,	   -0.0327,	   -0.0283,	   -0.0257,	   -0.0240,	   -0.0200;
                        0.0494,	    0.0509,	    0.0506,	    0.0560,	    0.0688,	    0.0773,	    0.0641,	    0.0499,	    0.0396,	    0.0215,	    0.0070,	   -0.0009,	   -0.0008;
                        0.1406,	    0.1442,	    0.1466,	    0.1532,	    0.1747,	    0.1859,	    0.1691,	    0.1399,    	0.1119,	    0.0713,	    0.0397,	    0.0222,	    0.0184;
                        0.2318,	    0.2375,	    0.2426,	    0.2504,	    0.2806,	    0.2945,	    0.2741,	    0.2299,	    0.1842,	    0.1211,	    0.0724,	    0.0453,	    0.0376;
                        0.3230,	    0.3308,	    0.3386,	    0.3476,	    0.3865,	    0.4031,	    0.3791,	    0.3199,	    0.2565,	    0.1709,	    0.1051,	    0.0684,	    0.0568;
                        0.4142,	    0.4241,	    0.4346,	    0.4448,	    0.4924,	    0.5117,	    0.4841,	    0.4099,	    0.3288,	    0.2207,	    0.1378,	    0.0915,    	0.0760;
                        0.5054, 	0.5174,	    0.5306,	    0.5420,	    0.5983,	    0.6203,	    0.5891,	    0.4999,	    0.4011,	    0.2705,	    0.1705,	    0.1146,	    0.0952;
                        0.5966,	    0.6107,	    0.6266,	    0.6392,	    0.7042,	    0.7289,	    0.6941,	    0.5899,	    0.4734,	    0.3203,	    0.2032,	    0.1377,	    0.1144];

% Drag Coefficients                    
auxdata.data_C_D =     [0.0421,	    0.0459,	    0.0492,	    0.0569,	    0.0711,	    0.0933,	    0.0903,	    0.0875,	    0.0769,	    0.0569,	    0.0400,	    0.0294,	    0.0244;
                        0.0353,  	0.0363,    	0.0399,	    0.0487,	    0.0613,	    0.0838,	    0.0824,	    0.0794,	    0.0702,	    0.0516,	    0.0360,	    0.0264,	    0.0219;
                        0.0343,	    0.0338,	    0.0387,	    0.0488,	    0.0610,	    0.0841,	    0.0830,	    0.0786,	    0.0692,	    0.0504,	    0.0347,	    0.0252,	    0.0209;
                        0.0391,	    0.0385,	    0.0454,	    0.0573,	    0.0703,	    0.0944,	    0.0922,	    0.0852,	    0.0742,	    0.0534,	    0.0362,	    0.0259,	    0.0215;
                        0.0497,	    0.0504,	    0.0602,	    0.0740,	    0.0889,	    0.1145,	    0.1101,	    0.0993,	    0.0849,	    0.0605,	    0.0403,	    0.0285,	    0.0236;
                        0.0662,	    0.0696,	    0.0831,	    0.0992,	    0.1171,	    0.1445,	    0.1364,	    0.1207,	    0.1015,	    0.0718,	    0.0472,	    0.0329,	    0.0272;
                        0.0884,	    0.0959,	    0.1139,	    0.1326,	    0.1548,	    0.1844,	    0.1714,	    0.1496,	    0.1239,	    0.0873,	    0.0569,	    0.0392,	    0.0323;
                        0.1164,	    0.1294,	    0.1528, 	0.1744,	    0.2019,	    0.2341,	    0.2150,	    0.1859,	    0.1522,	    0.1069,	    0.0692,	    0.0474,	    0.0390;
                        0.1503,	    0.1702,	    0.1997,	    0.2245,	    0.2586,	    0.2938,	    0.2672,	    0.2295,	    0.1863,	    0.1307,	    0.0843,	    0.0575,	    0.0472];

% Angle of Attack and Mach Number                    
auxdata.data_alpha =  [-3,			0,			3, 			6,			9,			12,			15,			18,			21]';     
auxdata.data_mach  =  [0.4,         0.6, 		0.8, 		0.9,		0.95,		1.05,		1.2,		1.5,		2,			3,			6,			12,			24];

% Inlet Capture Area Ratio
auxdata.data_CAR    =  [1.09449,   	0.53018,	0.31459,   	0.26226,   	0.24452,	0.22157,   	0.20981,   	0.23464,   	0.34159,   	0.62377,   	1.45141,   	2.76052,   	4.68122;
						1.18766,   	0.58304,   	0.35204,   	0.29597,   	0.27698,   	0.25238,   	0.23978,   	0.26638,   	0.38098,   	0.68332,   	1.57007,   	2.97269,   	5.03057;
						1.28082,   	0.63590,   	0.38950,   	0.32969,   	0.30943,   	0.28319,   	0.26975,   	0.29813,   	0.42037,   	0.74286,   	1.68873,   	3.18486,   	5.37993;
						1.37399,   	0.68875,   	0.42696,   	0.36341,   	0.34188,   	0.31401,   	0.29973,   	0.32987,   	0.45975,   	0.80240,   	1.80739,   	3.39702,   	5.72929;
						1.46715,   	0.74161,   	0.46441,   	0.39713,   	0.37433,   	0.34482,   	0.32970,   	0.36162,   	0.49914,   	0.86194,   	1.92605,   	3.60919,   	6.07865;
						1.56032,   	0.79447,   	0.50187,   	0.43085,   	0.40679,   	0.37563,   	0.35967,   	0.39337,   	0.53852,   	0.92148,   	2.04471,   	3.82136,   	6.42801;
						1.65348,   	0.84732,   	0.53933,   	0.46457,   	0.43924,   	0.40644,   	0.38964,  	0.42511,   	0.57791,   	0.98102,   	2.16336,   	4.03352,   	6.77737;
						1.74664,   	0.90018,   	0.57678,   	0.49829,   	0.47169,   	0.43726,   	0.41962,   	0.45686,   	0.61729,   	1.04057,   	2.28202,   	4.24569,   	7.12673;	
						1.83981,   	0.95304,   	0.61424,   	0.53201,   	0.50414,   	0.46807,   	0.44959,   	0.48860,   	0.65668,   	1.10011,   	2.40068,   	4.45786,   	7.47609];

% Throttle and Specific Impulse
auxdata.data_throttle =[0,			0.25,		0.5,		0.75,		1,			1.25,		1.5,		1.75,		2]';		
auxdata.data_I_sp   =  [0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000;
						1693.1500, 	1693.1500, 	1686.9000, 	1682.2125,	1679.0875, 	1661.9000, 	1636.9000, 	1568.1500, 	1443.1500, 	1568.1500, 	1318.1500,  755.6500,  	505.6500;
						2262.3999, 	2262.3999, 	2253.5625, 	2246.9343, 	2242.5156, 	2218.2126, 	2182.8625, 	2085.6499, 	1908.9000, 	2085.6499, 	1732.1499,  936.7750,  	583.2750;
						2699.6499, 	2699.6499, 	2688.8250, 	2680.7063, 	2675.2937, 	2645.5249, 	2602.2251, 	2483.1499, 	2266.6501, 	2483.1499, 	2050.1499, 	1075.9000,  642.9000;
						3068.1499, 	3068.1499, 	3055.6499, 	3046.2749, 	3040.0249, 	3005.6499, 	2955.6499, 	2818.1499, 	2568.1499, 	2818.1499, 	2318.1499, 	1193.1500,  693.1500;
						3392.6501, 	3392.6501, 	3378.6748, 	3368.1936, 	3361.2063, 	3322.7749, 	3266.8750, 	3113.1499, 	2833.6499, 	3113.1499, 	2554.1501, 	1296.4000,  737.4000;
						3686.8999, 	3686.8999, 	3671.5874, 	3660.1030, 	3652.4468, 	3610.3374, 	3549.0874, 	3380.6499, 	3074.3999, 	3380.6499, 	2768.1499,	1390.0250,  777.5250;
						3956.4001, 	3956.4001, 	3939.8625, 	3927.4595, 	3919.1907, 	3873.7124, 	3807.5625, 	3625.6499, 	3294.8999, 	3625.6499, 	2964.1501, 	1475.7750,  814.2750;
						4206.6499, 	4206.6499, 	4188.9751, 	4175.7188, 	4166.8813, 	4118.2749, 	4047.5750, 	3853.1499, 	3499.6499, 	3853.1499, 	3146.1499, 	1555.4000,  848.4000];
 
% Initial conditions for system.
alt0        = 80; 
lat0        = deg2rad(30); 
lon0        = deg2rad(30);
speed0      = 120; 
path_ang0   = deg2rad(4); 
head_ang0   = 0; 
m0          = 136078;

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0 = [alt0 lon0 lat0 speed0 path_ang0 head_ang0 m0];
problem.states.x0l = [alt0 lon0 lat0 speed0 path_ang0 -deg2rad(180) m0]; 
problem.states.x0u = [alt0 lon0 lat0 speed0 path_ang0 deg2rad(180) m0];


TOL = 1e-3;
% State bounds. xl=< x <=xu
alt_min         = 0;                    alt_max         = 84e3;
lon_min         = -deg2rad(180);        lon_max         = deg2rad(180);
lat_min         = -deg2rad(89);         lat_max         = deg2rad(89);
speed_min       = 115;                  speed_max       = 2800;
head_ang_min    = -deg2rad(180);        head_ang_max    = deg2rad(180);
path_ang_min    = -deg2rad(89);         path_ang_max    = deg2rad(89);
m_min           = 54431;                m_max           = m0;

% Input Bounds
alpha_min       = -deg2rad(3);          alpha_max       = deg2rad(21);
bank_ang_min    = -deg2rad(45);         bank_ang_max    = deg2rad(45);
u_alpha_min     = -deg2rad(6);          u_alpha_max     = deg2rad(6);
u_bank_min      = -deg2rad(6);          u_bank_max      = deg2rad(6);
u_mf_min        = 18;                   u_mf_max        = 163;


problem.states.xl = [alt_min, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, m_min];
problem.states.xu = [alt_max, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, m_max];

% State rate bounds. xrl=< x_dot <=xru
problem.states.xrl = [-inf -inf -inf -inf -inf -inf -inf]; 
problem.states.xru = [ inf  inf  inf  inf  inf  inf  inf]; 

% State error bounds
problem.states.xErrorTol_local      = [2 deg2rad(2e-5) deg2rad(2e-5) 2 deg2rad(0.8) deg2rad(0.8) 2]; 
problem.states.xErrorTol_integral   = [2 deg2rad(2e-5) deg2rad(2e-5) 2 deg2rad(0.8) deg2rad(0.8) 2]; 

% State constraint error bounds
problem.states.xConstraintTol       = [2 deg2rad(2e-5) deg2rad(2e-5) 2 deg2rad(0.8) deg2rad(0.8) 2]; 
problem.states.xrConstraintTol      = [2 deg2rad(2e-5) deg2rad(2e-5) 2 deg2rad(0.8) deg2rad(0.8) 2]; 

% Terminal state bounds. xfl=< xf <=xfu
altf        = 80;
lonf        = deg2rad(0);
latf        = deg2rad(0);
path_angf   = deg2rad(0);
head_angf   = deg2rad(0);

problem.states.xfl = [altf, lonf, latf, speed_min, path_ang_min, head_ang_min, m_min]; 
problem.states.xfu = [altf, lonf, latf, speed_max, path_ang_max, head_ang_max, m_max];

% Guess the state trajectories with [x0 ... xf]
guess.states(:,1)=[alt0 altf];
guess.states(:,2)=[lon0 lonf];
guess.states(:,3)=[lat0 latf];
guess.states(:,4)=[speed0 speed_max];
guess.states(:,5)=[path_ang0 path_angf];
guess.states(:,6)=[head_ang0 head_angf];
guess.states(:,7)=[m0 m_min];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul = [alpha_min bank_ang_min u_mf_min];
problem.inputs.uu = [alpha_max bank_ang_max u_mf_max];

% Bounds on the first control action
alpha0      = deg2rad(4);
bank_ang0   = 0;

problem.inputs.u0l = [alpha0 bank_ang0 u_mf_min];
problem.inputs.u0u = [alpha0 bank_ang0 u_mf_max];

% Input rate bounds
problem.inputs.url = [u_alpha_min u_bank_min -u_mf_max/6]; 
problem.inputs.uru = [u_alpha_max u_bank_max  u_mf_max/6]; 

% Input constraint error bounds
problem.inputs.uConstraintTol  = [deg2rad(0.8) deg2rad(0.8) 2];
problem.inputs.urConstraintTol = [deg2rad(0.8) deg2rad(0.8) 2];

% Guess the input sequences with [u0 ... uf]
guess.inputs(:,1) = [alpha0     0];
guess.inputs(:,2) = [bank_ang0  bank_ang0];
guess.inputs(:,3) = [u_mf_max   u_mf_max];

% Path constraint function 
problem.constraints.ng_eq       = 0; % number of quality constraints in format of g(x,u,p,t) == 0
problem.constraints.gTol_eq     = []; % equality cosntraint error bounds

problem.constraints.gl          = [0 0 -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf]; % Lower ounds for inequality constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gu          = [95.76e3 2 inf inf inf inf inf inf inf inf inf inf]; % Upper ounds for inequality constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gTol_neq    = [10 0.1 2*ones(1, 10)]; % inequality constraint error bounds

% Bounds for boundary constraints bl =< b(x0,xf,u0,uf,p,t0,tf) =< bu
problem.constraints.bl=0;
problem.constraints.bu=0;
problem.constraints.bTol=0; 

% store the necessary problem parameters used in the functions
problem.data.auxdata=auxdata;

% Get function handles and return to Main.m
problem.data.InternalDynamics=InternalDynamics;
problem.data.functionfg=@fg;
problem.data.plantmodel = func2str(InternalDynamics);
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.sim.functions=SimDynamics;
problem.sim.inputX=[];
problem.sim.inputU=1:length(problem.inputs.ul);
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc,@b_unscaled};
problem.data.functions_unscaled=problem.functions_unscaled;
problem.data.ng_eq=problem.constraints.ng_eq;
problem.constraintErrorTol=[problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.constraints.gTol_eq,problem.constraints.gTol_neq,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,xr,u,ur,p,t,data)
%------------- BEGIN CODE --------------
auxdata         = data.auxdata;
ee              = auxdata.ee;

u1=u(:, 1);
u2=u(:, 2);
u3=u(:, 3);

stageCost = (1-ee)*(u1.*u1+u2.*u2+u3.*u3);
%------------- END OF CODE --------------

function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 
%------------- BEGIN CODE --------------
auxdata         = data.auxdata;
ee              = auxdata.ee;

boundaryCost = ee*tf;
%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,data,varargin)
% Leave it here
varargin=varargin{1};
%------------- BEGIN CODE --------------
bc(1,:)=0;
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if strcmp(options.transcription,'hpLGR') && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

