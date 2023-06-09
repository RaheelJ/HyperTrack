#include<psopt.h>
using namespace PSOPT;

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

struct Constants{
	double A_ref;								//Aerodyanmic Reference Area
	double Re;									//Radius of Earth
	double omega, mu, g0;						//Earth's rotation rate, gravitational parameter, gravitational acceleration
	double phi_st, A_ic;						//Stoichiometric fuel-air ratio and inlet capture area
	double ee;									//Convex combination coefficient
	double R[3];								//Weight Matrix
	MatrixXd *C_L;								//Data points for calculation of the lift coefficient
	MatrixXd *C_D;								//Data points for calculation of the drag coefficient
	MatrixXd *mach;								//Data points for mach number
	MatrixXd *alpha, *I_sp, *CAR;				//Data points for calculation of thrust
	MatrixXd *throttle;							//Fuel-air equivalence ratio
};
typedef struct Constants Constants_;

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

//Auxiliary Functions//
adouble rad_to_deg(adouble in) 
{
    return (180 * in / M_PI);
}
adouble deg_to_rad(adouble in) 
{
    return (M_PI * in / 180);
}
double rad_to_deg(double in) 
{
    return (180 * fmod(in, 2 * M_PI) / M_PI);
}
double deg_to_rad(double in) 
{
    return (M_PI * fmod(in, 360) / 180);
}
adouble Limit(adouble in, double upper_lim, double lower_lim)
{
	adouble out;
    if (in < lower_lim) { out = lower_lim; }
    else if (in > upper_lim) { out = upper_lim; }
	else { out = in; }
	return out;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

//Density conversion from (Kg/m^3) to (lbs/ft^3)//
double conv_den(double den_metric)
{
	double den_imperial;
	den_imperial=den_metric*0.062427960576145;
	return den_imperial;
}

//Distance conversion from (meters) to (foot)//
double conv_dist(double meters)
{
	double feet=meters*3.28084;
	return feet;
}

//Calculate Ambient Density (lbs/ft^3)
double amb_den(adouble alt)
{
	double scale_height=conv_dist(10.4e3);
	double atm_den_sea=conv_den(1.225);
	double density=atm_den_sea*exp(-alt.value()/scale_height);
	return density;
}

//Calculate speed of sound based on curve fitted over 1976 US standard atmospheric data
double speed_sound(adouble alt)
{
	double p1, p2, p3, p4, p5, p6, p7, p8;
	double h, speed;
	
	p1=-127.3;
	p2=45.35;
	p3=439.1;
	p4=-21.11;
	p5=-514.2;
	p6=-103.7;
	p7=200.7;
	p8=1046;
	
	h=(alt.value()-1.312e5)/9.772e4;
	speed=+p1*pow(h, 7)+p2*pow(h, 6)+p3*pow(h, 5)+p4*pow(h, 4)+p5*pow(h, 3)+p6*pow(h, 2)+p7*h+p8;
	return speed;
}


	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
					  adouble* parameters,adouble& t0, adouble& tf,
					  adouble* xad, int iphase, Workspace* workspace)
{
	return 0;
}
adouble integrand_cost(adouble* states, adouble* controls,
					  adouble* parameters, adouble& time, adouble* xad,
					  int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble cost;	//Integrand cost function
	adouble V=states[3];
	
	cost=pow(controls[0]/deg_to_rad(20), 2)+pow(controls[1]/deg_to_rad(45), 2);
	return cost;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

//Differential Equations to Define the State Dyanmics// 
void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);	//Parameters from Workspace
	adouble h, lon, lat, V, path_ang, head_ang, m, alpha, bank_ang;			//States
	adouble u_alpha, u_bank, u_mf;											//Conrol Inputs
	adouble F_T, F_N;														//Tangent and normal aerodynamic forces
	adouble C_L;															//Lift coefficients																																					
	adouble q;																//Dynamic pressure	
	adouble mach;															//Mach number
	adouble C_D;															//Drag coefficients
	adouble I_sp;															//Specific impulse
	adouble CAR, throttle;													//Capture area ratio, fuel-air equivalence ratio
	adouble term1, term2, term3, term4;
	adouble alpha_deg;
	double L, D, T;															//Lift, Drag, Thrust
	double c, r, g;															//Speed of sound, Altitude from Earth's centre, local gravity
	
	h=states[0];					//Altitude 								(ft)
	lon=states[1];					//Longitude								(rad)
	lat=states[2];					//Latitude								(rad)
	V=states[3];					//Speed									(ft/s)
	path_ang=states[4];				//Flight Path Angle						(rad)
	head_ang=states[5];				//Velocity Heading Angle				(rad)
	m=states[6];					//Mass of Vehicle						(lbs)
	alpha=states[7];				//Angle of Attack (AoA)					(rad)
	bank_ang=states[8];				//Bank Angle							(rad)
	
	u_alpha=controls[0];			//Rate of AoA							(rad/s)
	u_bank=controls[1];				//Rate of Bank Angle					(rad/s)
	//u_mf=controls[2];				//Mass Flow Rate of Propellant			(lbs/s)
	
	alpha_deg=rad_to_deg(alpha);
	r=h.value()+Param.Re;
	g=Param.mu/pow(r, 2);			
	q=amb_den(h)*V*V/2;
	
	//Calculation of the mach number
	c=speed_sound(h);
	mach=V/c;
	
	//cout << "\n" << "mach: " << mach << "\t" << "h: " << h << "\n";
	//cout << "\n" << "alpha: " << alpha_deg << "\n";
	
	//Calculation of the lift force
	MatrixXd table_mach=*Param.mach;
	MatrixXd table_C_L=*Param.C_L;
	MatrixXd table_alpha=*Param.alpha;
	smooth_bilinear_interpolation(&C_L, alpha_deg, mach, table_alpha, table_mach, table_C_L);
	L=C_L.value()*q.value()*Param.A_ref;
	
	//Calculation of the drag force
	MatrixXd table_C_D=*Param.C_D;
	smooth_bilinear_interpolation(&C_D, alpha_deg, mach, table_alpha, table_mach, table_C_D);
	D=C_D.value()*q.value()*Param.A_ref;
	
	//Calculation of thrust
	/* MatrixXd table_CAR=*Param.CAR;
	smooth_bilinear_interpolation(&CAR, alpha_deg, mach, table_alpha, table_mach, table_CAR);
	throttle=u_mf/(amb_den(h)*V*CAR*Param.A_ic*Param.phi_st);

	MatrixXd table_throttle=*Param.throttle;
	MatrixXd table_I_sp=*Param.I_sp;
	smooth_bilinear_interpolation(&I_sp, throttle, mach, table_throttle, table_mach, table_I_sp); */
	//T=Param.g0*u_mf.value()*I_sp.value();
	I_sp=300;
	T=2880000;
	
	//cout << "I_sp: " << I_sp << "\t" << "throttle: " << throttle << "\n";
		
	//Calculation of the tangent and normal acting forces
	F_T=T*cos(alpha)-D;
	F_N=T*sin(alpha)+L;
	
	//State Dynamics//
	derivatives[0]=V*sin(path_ang);
	derivatives[1]=V*cos(path_ang)*cos(head_ang)/(r*cos(lat));
	derivatives[2]=V*cos(path_ang)*sin(head_ang)/r;
	
	term1=(1/m)*F_T;
	term2=g*sin(path_ang);
	term3=r*pow(Param.omega, 2)*cos(lat)*(sin(path_ang)*cos(lat)-cos(path_ang)*sin(lat)*sin(head_ang));
	derivatives[3]=term1-term2+term3;
	
	term1=F_N*cos(bank_ang)/(m*V);
	term2=cos(path_ang)*((V/r)-(g/V));
	term3=2*Param.omega*cos(lat)*cos(head_ang);
	term4=(r/V)*pow(Param.omega, 2)*cos(lat)*(cos(path_ang)*cos(lat)+sin(path_ang)*sin(lat)*sin(head_ang));
	derivatives[4]=term1+term2+term3+term4;
	
	term1=F_N*(sin(bank_ang)/(m*V*cos(path_ang)));
	term2=(V/r)*cos(path_ang)*cos(head_ang)*tan(lat);
	term3=2*Param.omega*(tan(path_ang)*cos(lat)*sin(head_ang)-sin(lat));
	term4=(r*pow(Param.omega, 2)/(V*cos(path_ang)))*sin(lat)*cos(lat)*cos(head_ang);
	derivatives[5]=term1-term2+term3-term4;
	
	if(T>0)
	{
		derivatives[6]=-T/(I_sp*Param.g0);
	}
	else
	{
		derivatives[6]=0;
	}
	derivatives[7]=u_alpha;			
	derivatives[8]=u_bank;
	
	path[0]=throttle;
	path[0]=q;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void events(adouble* e, adouble* initial_states, adouble* final_states,
			adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
			int iphase, Workspace* workspace)
{
	e[0]=initial_states[0];
	e[1]=initial_states[1];
	e[2]=initial_states[2];
	e[3]=initial_states[3];
	e[4]=initial_states[4];
	e[5]=initial_states[5];
	e[6]=initial_states[6];
	e[7]=initial_states[7];
	e[8]=initial_states[8];

	e[9]=final_states[0];
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
	
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	

int main()
{
	int i, j;
	double TOL = 0.02;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//Declare key structures
	Alg algorithm;
	Sol solution;
	Prob problem;
	
	//Register problem name
	problem.name="Air-Breathing Vehicle Mission";
	problem.outfilename="Trajectory.txt";
	
	//Instance of Parameters structures
	Constants_ Param;
	problem.user_data=(void*) &Param;
	
	//Problem level 1 setup
	problem.nphases=1;
	problem.nlinkages=0;
	psopt_level1_setup(problem);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	

	//Level 2 setup
	problem.phases(1).nstates=9;
	problem.phases(1).ncontrols=2;
	problem.phases(1).nevents=10;
	problem.phases(1).npath=1;
	problem.phases(1).nparameters=0;
	problem.phases(1).nodes = (RowVectorXi(10) << 8, 12, 18, 24, 32, 38, 46, 58, 66, 78).finished();

	psopt_level2_setup(problem, algorithm);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//Initialize Parameters
	MatrixXd data_C_L(9, 13), data_mach(1, 13);
	MatrixXd data_C_D(9, 13);
	MatrixXd data_CAR(9, 13), data_alpha(1, 9);
	MatrixXd data_I_sp(9, 13), data_throttle(1, 9);
	
	Param.omega=0*0.00007292115856;			//(rad/s) 
	Param.mu=1.40764e16;					//(ft^3/s^2)
	Param.g0=32.174;						//(ft/s^2)						
	Param.Re=conv_dist(6378.166*1000);		//(ft)
	Param.A_ref=6000;						//(ft^2)
	Param.A_ic=300;							//(ft^2)
	Param.phi_st=0.292;
	Param.ee=00;

	Param.R[0]=1, Param.R[1]=1, Param.R[2]=1;

	
	//GHAME 3-DOF Aerodynamic Model (Lift)
	data_mach		<<	0.4,		0.6, 		0.8, 		0.9,		0.95,		1.05,		1.2,		1.5,		2,			3,			6,			12,			24;
	data_C_L		<< -0.1330,	   -0.1357,	   -0.1414,	   -0.1384,	   -0.1430,	   -0.1399,	   -0.1459,	   -0.1301,	   -0.1050,	   -0.0781,	   -0.0584,	   -0.0471,	   -0.0392,
					   -0.0418,	   -0.0424,	   -0.0454,	   -0.0412,	   -0.0371,	   -0.0313,	   -0.0409,	   -0.0401,	   -0.0327,	   -0.0283,	   -0.0257,	   -0.0240,	   -0.0200,
						0.0494,	    0.0509,	    0.0506,	    0.0560,	    0.0688,	    0.0773,	    0.0641,	    0.0499,	    0.0396,	    0.0215,	    0.0070,	   -0.0009,	   -0.0008,
						0.1406,	    0.1442,	    0.1466,	    0.1532,	    0.1747,	    0.1859,	    0.1691,	    0.1399,    	0.1119,	    0.0713,	    0.0397,	    0.0222,	    0.0184,
						0.2318,	    0.2375,	    0.2426,	    0.2504,	    0.2806,	    0.2945,	    0.2741,	    0.2299,	    0.1842,	    0.1211,	    0.0724,	    0.0453,	    0.0376,
						0.3230,	    0.3308,	    0.3386,	    0.3476,	    0.3865,	    0.4031,	    0.3791,	    0.3199,	    0.2565,	    0.1709,	    0.1051,	    0.0684,	    0.0568,
						0.4142,	    0.4241,	    0.4346,	    0.4448,	    0.4924,	    0.5117,	    0.4841,	    0.4099,	    0.3288,	    0.2207,	    0.1378,	    0.0915,    	0.0760,
						0.5054, 	0.5174,	    0.5306,	    0.5420,	    0.5983,	    0.6203,	    0.5891,	    0.4999,	    0.4011,	    0.2705,	    0.1705,	    0.1146,	    0.0952,
						0.5966,	    0.6107,	    0.6266,	    0.6392,	    0.7042,	    0.7289,	    0.6941,	    0.5899,	    0.4734,	    0.3203,	    0.2032,	    0.1377,	    0.1144;
	Param.C_L=&data_C_L;
	Param.mach=&data_mach;
	
	//GHAME 3-DOF Aerodynamic Model (Drag)
	data_C_D		<<	0.0421,	    0.0459,	    0.0492,	    0.0569,	    0.0711,	    0.0933,	    0.0903,	    0.0875,	    0.0769,	    0.0569,	    0.0400,	    0.0294,	    0.0244,
						0.0353,  	0.0363,    	0.0399,	    0.0487,	    0.0613,	    0.0838,	    0.0824,	    0.0794,	    0.0702,	    0.0516,	    0.0360,	    0.0264,	    0.0219,
						0.0343,	    0.0338,	    0.0387,	    0.0488,	    0.0610,	    0.0841,	    0.0830,	    0.0786,	    0.0692,	    0.0504,	    0.0347,	    0.0252,	    0.0209,
						0.0391,	    0.0385,	    0.0454,	    0.0573,	    0.0703,	    0.0944,	    0.0922,	    0.0852,	    0.0742,	    0.0534,	    0.0362,	    0.0259,	    0.0215,
						0.0497,	    0.0504,	    0.0602,	    0.0740,	    0.0889,	    0.1145,	    0.1101,	    0.0993,	    0.0849,	    0.0605,	    0.0403,	    0.0285,	    0.0236,
						0.0662,	    0.0696,	    0.0831,	    0.0992,	    0.1171,	    0.1445,	    0.1364,	    0.1207,	    0.1015,	    0.0718,	    0.0472,	    0.0329,	    0.0272,
						0.0884,	    0.0959,	    0.1139,	    0.1326,	    0.1548,	    0.1844,	    0.1714,	    0.1496,	    0.1239,	    0.0873,	    0.0569,	    0.0392,	    0.0323,
						0.1164,	    0.1294,	    0.1528, 	0.1744,	    0.2019,	    0.2341,	    0.2150,	    0.1859,	    0.1522,	    0.1069,	    0.0692,	    0.0474,	    0.0390,
						0.1503,	    0.1702,	    0.1997,	    0.2245,	    0.2586,	    0.2938,	    0.2672,	    0.2295,	    0.1863,	    0.1307,	    0.0843,	    0.0575,	    0.0472;
	Param.C_D=&data_C_D;
	//GHAME 3-DOF Propuslsion Model
	data_CAR		<<	1.09449,   	0.53018,	0.31459,   	0.26226,   	0.24452,	0.22157,   	0.20981,   	0.23464,   	0.34159,   	0.62377,   	1.45141,   	2.76052,   	4.68122,
						1.18766,   	0.58304,   	0.35204,   	0.29597,   	0.27698,   	0.25238,   	0.23978,   	0.26638,   	0.38098,   	0.68332,   	1.57007,   	2.97269,   	5.03057,
						1.28082,   	0.63590,   	0.38950,   	0.32969,   	0.30943,   	0.28319,   	0.26975,   	0.29813,   	0.42037,   	0.74286,   	1.68873,   	3.18486,   	5.37993,
						1.37399,   	0.68875,   	0.42696,   	0.36341,   	0.34188,   	0.31401,   	0.29973,   	0.32987,   	0.45975,   	0.80240,   	1.80739,   	3.39702,   	5.72929,
						1.46715,   	0.74161,   	0.46441,   	0.39713,   	0.37433,   	0.34482,   	0.32970,   	0.36162,   	0.49914,   	0.86194,   	1.92605,   	3.60919,   	6.07865,
						1.56032,   	0.79447,   	0.50187,   	0.43085,   	0.40679,   	0.37563,   	0.35967,   	0.39337,   	0.53852,   	0.92148,   	2.04471,   	3.82136,   	6.42801,
						1.65348,   	0.84732,   	0.53933,   	0.46457,   	0.43924,   	0.40644,   	0.38964,  	0.42511,   	0.57791,   	0.98102,   	2.16336,   	4.03352,   	6.77737,
						1.74664,   	0.90018,   	0.57678,   	0.49829,   	0.47169,   	0.43726,   	0.41962,   	0.45686,   	0.61729,   	1.04057,   	2.28202,   	4.24569,   	7.12673,	
						1.83981,   	0.95304,   	0.61424,   	0.53201,   	0.50414,   	0.46807,   	0.44959,   	0.48860,   	0.65668,   	1.10011,   	2.40068,   	4.45786,   	7.47609;
	data_alpha		<<	-3,			0,			3, 			6,			9,			12,			15,			18,			21;
	Param.CAR=&data_CAR;
	Param.alpha=&data_alpha;
	
	data_throttle	<<	0,			0.25,		0.5,		0.75,		1,			1.25,		1.5,		1.75,		2;		
	data_I_sp		<<	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,    	0.0000,
						1693.1500, 	1693.1500, 	1686.9000, 	1682.2125,	1679.0875, 	1661.9000, 	1636.9000, 	1568.1500, 	1443.1500, 	1568.1500, 	1318.1500,  755.6500,  	505.6500,
						2262.3999, 	2262.3999, 	2253.5625, 	2246.9343, 	2242.5156, 	2218.2126, 	2182.8625, 	2085.6499, 	1908.9000, 	2085.6499, 	1732.1499,  936.7750,  	583.2750,
						2699.6499, 	2699.6499, 	2688.8250, 	2680.7063, 	2675.2937, 	2645.5249, 	2602.2251, 	2483.1499, 	2266.6501, 	2483.1499, 	2050.1499, 	1075.9000,  642.9000,
						3068.1499, 	3068.1499, 	3055.6499, 	3046.2749, 	3040.0249, 	3005.6499, 	2955.6499, 	2818.1499, 	2568.1499, 	2818.1499, 	2318.1499, 	1193.1500,  693.1500,
						3392.6501, 	3392.6501, 	3378.6748, 	3368.1936, 	3361.2063, 	3322.7749, 	3266.8750, 	3113.1499, 	2833.6499, 	3113.1499, 	2554.1501, 	1296.4000,  737.4000,
						3686.8999, 	3686.8999, 	3671.5874, 	3660.1030, 	3652.4468, 	3610.3374, 	3549.0874, 	3380.6499, 	3074.3999, 	3380.6499, 	2768.1499,	1390.0250,  777.5250,
						3956.4001, 	3956.4001, 	3939.8625, 	3927.4595, 	3919.1907, 	3873.7124, 	3807.5625, 	3625.6499, 	3294.8999, 	3625.6499, 	2964.1501, 	1475.7750,  814.2750,
						4206.6499, 	4206.6499, 	4188.9751, 	4175.7188, 	4166.8813, 	4118.2749, 	4047.5750, 	3853.1499, 	3499.6499, 	3853.1499, 	3146.1499, 	1555.4000,  848.4000;
	Param.throttle=&data_throttle;						
	Param.I_sp=&data_I_sp;
										
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//State Limits
	double h_min=0, h_max=130000;											//feet
	double lon_min=deg_to_rad(-180), lon_max=deg_to_rad(180);				//rad
	double lat_min=deg_to_rad(-89), lat_max=deg_to_rad(89);					//rad
	double V_min=379, V_max=8931.6;											//ft/sec
	double path_ang_min=deg_to_rad(-89), path_ang_max=deg_to_rad(89); 		//rad
	double head_ang_min=deg_to_rad(-180), head_ang_max=deg_to_rad(180);		//rad
	double m_min=120000, m_max=300000;										//lb
	double alpha_min=deg_to_rad(-3), alpha_max=deg_to_rad(21);				//rad
	double bank_ang_min=deg_to_rad(-45), bank_ang_max=deg_to_rad(45);		//rad
	
	//Control Limits
	double u_mf_min=40, u_mf_max=360;										//lb/sec
	double u_alpha_min=deg_to_rad(-6), u_alpha_max=deg_to_rad(6);			//rad/sec
	double u_bank_min=deg_to_rad(-6), u_bank_max=deg_to_rad(6);				//rad/sec
	
	//Time Limits
	double t0=0, tf=500;
	double throttle_min=0, throttle_max=2;
	double q_min=0, q_max=2000;
	
	//Problem bounds information for Phase 1
	problem.phases(1).bounds.lower.states << h_min, lon_min, lat_min, V_min, path_ang_min, head_ang_min, m_min, alpha_min, bank_ang_min;
	problem.phases(1).bounds.upper.states << h_max, lon_max, lat_max, V_max, path_ang_max, head_ang_max, m_max, alpha_max, bank_ang_max;

	problem.phases(1).bounds.lower.controls(0) = u_alpha_min;
	problem.phases(1).bounds.upper.controls(0) = u_alpha_max;
	problem.phases(1).bounds.lower.controls(1) = u_bank_min;
	problem.phases(1).bounds.upper.controls(1) = u_bank_max;
	/* problem.phases(1).bounds.lower.controls(2) = u_mf_min;
	problem.phases(1).bounds.upper.controls(2) = u_mf_max; */
	
	problem.phases(1).bounds.lower.StartTime = t0;
	problem.phases(1).bounds.upper.StartTime = t0;
	problem.phases(1).bounds.lower.EndTime = 50;
	problem.phases(1).bounds.upper.EndTime = tf;
	
	problem.phases(1).bounds.lower.path << q_min;
	problem.phases(1).bounds.upper.path << q_max;

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
			
	//Initial Conditions
	double h0=0, hf=0;												//feet
	double lon0=deg_to_rad(0), lonf=deg_to_rad(-70);				//rad
	double lat0=deg_to_rad(0), latf=deg_to_rad(46);					//rad
	double V0=379, Vf=379;											//ft/s
	double path_ang0=deg_to_rad(0), path_angf=deg_to_rad(0);		//rad
	double head_ang0=deg_to_rad(0), head_angf=deg_to_rad(50);		//rad
	double m0=300000;												//lbs
	double alpha0=deg_to_rad(0), alphaf=deg_to_rad(0);				//rad				
	double bank_ang0=deg_to_rad(10), bank_angf=deg_to_rad(0);		//rad
			
	//Events
	problem.phases(1).bounds.lower.events << 	h0, lon0, lat0, V0, path_ang0, head_ang0, m0, alpha0, bank_ang0, 100000;	
	problem.phases(1).bounds.upper.events << 	h0, lon0, lat0, V0, path_ang0, head_ang0, m0, alpha0, bank_ang0, 100000;
												
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Initial Guess Phase 1
	problem.phases(1).guess.states=zeros(9,16);
	problem.phases(1).guess.states.row(0)=linspace(h0, h_max, 16);
	problem.phases(1).guess.states.row(1)=linspace(lon0, lon_max, 16);
	problem.phases(1).guess.states.row(2)=linspace(lat0, lat0, 16);
	problem.phases(1).guess.states.row(3)=linspace(V0, V_max, 16);
	problem.phases(1).guess.states.row(4)=linspace(path_ang0, path_ang_max, 16);
	problem.phases(1).guess.states.row(5)=linspace(head_ang0, head_ang0, 16);
	problem.phases(1).guess.states.row(6)=linspace(m0, m_min, 16);
	problem.phases(1).guess.states.row(7)=linspace(alpha0, alpha_max, 16);
	problem.phases(1).guess.states.row(8)=linspace(bank_ang0, bank_ang0, 16);

	problem.phases(1).guess.controls = zeros(3,16);
	problem.phases(1).guess.controls.row(0) = 0.3*ones(1, 16);
	problem.phases(1).guess.controls.row(1) = 0*ones(1, 16);
	//problem.phases(1).guess.controls.row(1) = u_mf_max*ones(1, 16);

	problem.phases(1).guess.time = linspace(t0, tf, 16);
		
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Register Problem Functions
	problem.integrand_cost=&integrand_cost;
	problem.endpoint_cost=&endpoint_cost;
	problem.dae=&dae;
	problem.events=&events;
	problem.linkages=&linkages;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Algorithm Options
	algorithm.nlp_method="IPOPT";
	algorithm.scaling="automatic";
	//algorithm.derivatives="numerical";
	algorithm.nlp_iter_max=1000;
	algorithm.diff_matrix = "central-differences";
	//algorithm.mesh_refinement="automatic";
	algorithm.collocation_method="Hermite-Simpson";
	algorithm.ipopt_linear_solver="ma27";
	algorithm.defect_scaling = "jacobian-based";
	//algorithm.mr_max_iterations=14;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Calling PSOPT to solve problems
	psopt(solution, problem, algorithm);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Extracting relevant variables
	MatrixXd x_ph1;
	MatrixXd t_ph1;
	MatrixXd u_ph1;
	MatrixXd x, u, t;
	
	x_ph1 = solution.get_states_in_phase(1);
	u_ph1 = solution.get_controls_in_phase(1);
	t_ph1 = solution.get_time_in_phase(1);
	
	x.resize(9, x_ph1.cols());
	x.row(0) << x_ph1.row(0);
	x.row(1) << x_ph1.row(1);
	x.row(2) << x_ph1.row(2);
	x.row(3) << x_ph1.row(3);
	x.row(4) << x_ph1.row(4);
	x.row(5) << x_ph1.row(5);
	x.row(6) << x_ph1.row(6);
	x.row(7) << x_ph1.row(7);
	x.row(8) << x_ph1.row(8);
	
	t.resize(1, t_ph1.cols());
	t << t_ph1;
	
	/* u.resize(3, u_ph1.cols());
	u.row(0)<<u_ph1.row(0);
	u.row(1)<<u_ph1.row(1);
	u.row(2)<<u_ph1.row(2); */
	
	Save(x,"x.dat");
	Save(t,"t.dat");
	Save(u,"u.dat");
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	MatrixXd altitude, longitude, latitude, speed_state, path_angle, heading_angle, attack_angle, bank_angle, mass;
	MatrixXd attack_angle_rate, bank_angle_rate, flow_rate;
    
	altitude = x.block(0,0,1,x_ph1.cols());
	longitude = x.block(1,0,1,x_ph1.cols()); 
	latitude = x.block(2,0,1,x_ph1.cols());
	speed_state = x.block(3,0,1,x_ph1.cols());
    path_angle = x.block(4,0,1,x_ph1.cols());
	heading_angle = x.block(5,0,1,x_ph1.cols());
	mass = x.block(6,0,1,x_ph1.cols());
	attack_angle = x.block(7,0,1,x_ph1.cols());
	bank_angle = x.block(8,0,1,x_ph1.cols());
	
    plot(t,altitude,problem.name, "time (s)", "Altitude (ft)");
	plot(t,longitude,problem.name, "time (s)", "Longitude (rad)");
	plot(t,latitude,problem.name, "time (s)", "Latitude (rad)");
	plot(t,speed_state,problem.name, "time (s)", "Speed (ft/s)");
	plot(t,path_angle,problem.name, "time (s)", "Path Angle (rad)");
	plot(t,heading_angle,problem.name, "time (s)", "Heading Angle (rad)");
	plot(t,attack_angle,problem.name, "time (s)", "Angle of Attack (rad)");
	plot(t,bank_angle,problem.name, "time (s)", "Bank Angle (rad)");
	plot(t,mass,problem.name, "time (s)", "Mass (lbs)");
	
	return 0;
}
