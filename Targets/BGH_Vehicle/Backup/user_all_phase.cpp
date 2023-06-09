#include<psopt.h>
using namespace PSOPT;
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

struct Constants{
	double Re, T[3]; 							//Radius of Earth, Thrust
	double I_sp[3];								//Specific impulse of boost vehicle
	double S[2];								//Vehicle reference area
	double alpha_bar[8], alpha_max[8];			//Nominal and maximum value of angle of attack
	double u_alpha_max[8];						//Maximum rate of angle of attack
	double u_bank_max[8];						//Maximum rate of bank angle
	//Lift and drag coefficients at low speed 
	MatrixXd* C_L_low;
	MatrixXd* C_D_low;
	//Angles of attack and mach numbers
	MatrixXd* alpha_low;
	MatrixXd* alpha_high;
	MatrixXd* mach_low;
	MatrixXd* mach_high;
	//Lift and Drag coefficients at high speed
	MatrixXd* C_L_high;
	MatrixXd* C_D_high;
	double omega_e, mu_e, g;					//Earth's rotation rate, Gravitational parameter, gravitational acceleration
};
typedef struct Constants Constants_;

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

//Auxiliary Functions
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
    return (180 * fmod((in), 2 * M_PI) / M_PI);
}
double deg_to_rad(double in) 
{
    return (M_PI * fmod((in), 360) / 180);
}
adouble Limit(adouble in, double upper_lim, double lower_lim)
{
	adouble out;
    if (in <= lower_lim) { out = lower_lim; }
    else if (in >= upper_lim) { out = upper_lim; }
	else { out = in; }
	return out;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

adouble amb_den(adouble alt)
{
	double scale_height=10.4e3;
	double atm_den_sea=1.225;
	adouble density=atm_den_sea*exp(-alt/scale_height);
	return density;
}
adouble map_path_ang(adouble ee_1, adouble ee_2, adouble ee_3, adouble nn)
{
	adouble path_ang;
	path_ang=atan2(0.5-ee_2*ee_2-ee_3*ee_3, sqrt((ee_1*ee_1+nn*nn)*(ee_2*ee_2+ee_3*ee_3)));
	return path_ang;
}
adouble map_head_ang(adouble ee_1, adouble ee_2, adouble ee_3, adouble nn)
{
	adouble head_ang;
	head_ang=atan2(ee_1*ee_2+ee_3*nn, ee_1*ee_3-ee_2*nn);
	return head_ang;
}
adouble map_bank_ang(adouble ee_1, adouble ee_2, adouble ee_3, adouble nn)
{
	adouble bank_ang;
	bank_ang=atan2(-ee_3*ee_1-ee_2*nn, ee_2*ee_1-ee_3*nn);
	return bank_ang;
}
adouble map_bank_rate(adouble ee_1, adouble ee_2, adouble ee_3, adouble nn,
					  adouble omega_1, adouble omega_2, adouble omega_3)
{
	adouble bank_rate;
	adouble temp1, temp2;
	temp1=(0.5-ee_2*ee_2-ee_3*ee_3)/((ee_1*ee_1+nn*nn)*(ee_2*ee_2+ee_3*ee_3));
	temp2=(omega_2*(ee_2*ee_1-ee_3*nn)+omega_3*(ee_3*ee_1+ee_2*nn));
	bank_rate=omega_1-temp1*temp2;
	return bank_rate;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
					  adouble* parameters,adouble& t0, adouble& tf,
					  adouble* xad, int iphase, Workspace* workspace)
{
	//std::cout<<"In Endpoint Cost Function"<<std::endl;
	return 0.0;
}
adouble integrand_cost(adouble* states, adouble* controls,
					  adouble* parameters, adouble& time, adouble* xad,
					  int iphase, Workspace* workspace)
{
	//std::cout<<"In Integrand Cost Function"<<std::endl;
	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble Cost;								//Integrand cost function
	adouble term1, term2, term3;				//Temporary variables
	adouble alpha;
	adouble u_alpha, omega_1, u_bank;			//Control inputs
	double u_alpha_max, omega_1_max, alpha_max, alpha_bar, u_bank_max;
	
	u_alpha=controls[0];
	u_alpha_max=deg_to_rad(Param.u_alpha_max[iphase-1]);
	term2=u_alpha/u_alpha_max;
	
	if(iphase==1 || iphase==8)
	{
		alpha=states[8];
		alpha_bar=deg_to_rad(Param.alpha_bar[iphase-1]);
		alpha_max=deg_to_rad(Param.alpha_max[iphase-1]);
		term1=(alpha-alpha_bar)/alpha_max;
		
		omega_1=controls[1];
		omega_1_max=deg_to_rad(Param.u_bank_max[iphase-1]);
		term3=omega_1/omega_1_max;
	}
	else
	{
		alpha=states[6];
		alpha_bar=deg_to_rad(Param.alpha_bar[iphase-1]);
		alpha_max=deg_to_rad(Param.alpha_max[iphase-1]);
		term1=(alpha-alpha_bar)/alpha_max;
		if(iphase==5 || iphase==6) { term1=0; }
		
		u_bank=controls[1];
		u_bank_max=deg_to_rad(Param.u_bank_max[iphase-1]);	
		term3=u_bank/u_bank_max;
	}

	Cost=pow(term1, 2)+pow(term2, 2)+pow(term3, 2);
	
	//std::cout<<"In Integrand Cost Function"<<std::endl;
	return Cost;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase, Workspace* workspace)
{
	//std::cout<<"In DAE Function"<<std::endl;
	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble C_L, C_D;						//Lift and Drag coefficients
	adouble L, D, q;						//Lift, Drag, Dynamic Pressure
	adouble L_temp, D_temp;
	double T, I_sp, S;						//Thrust, Reference Area
	adouble omega_2, omega_3;				//Angular velocity component
	adouble term1, term2, term3, term4, term5;
	adouble acc_sen;						//Sensed acceleration
	adouble Q_dot;							//Stagnation point heating rate
	adouble alpha_deg;						//Angle of Attack in degrees
	adouble temp_mass, mach;
	adouble u_alpha, u_bank, omega_1;		//Control inputs	
	/*	States of the system
		alt 					= Altitude 
		lon 					= Longitude
		lat 					= Latitude
		speed					= Magnitude of velocity
		ee_1, ee_2, ee_3, nn 	= Euler's parameters
		path_ang				= Flight Path Angle
		head_ang				= Velocity heading angle
		m						= Mass of vehicle
		alpha					= Angle of attack
		bank_ang				= Banking angle
	*/
	adouble alt, lon, lat, speed, ee_1, ee_2, ee_3, nn, path_ang, head_ang, m, alpha, bank_ang;
	
	u_alpha=controls[0];
	alt=states[0];					//Altitude
	lon=states[1];					//Longitude
	lat=states[2];					//Latitude
	speed=states[3];				//Speed
	adouble r=alt+Param.Re;			//Radial distance of vehicle from centre of the Earth
			
	//Interpolation to calculate Lift and Drag coefficients 
	alpha_deg=Limit(rad_to_deg(alpha), 25, 0);
	MatrixXd& alpha_low=*Param.alpha_low;
	MatrixXd& mach_low=*Param.mach_low;
	MatrixXd& C_L_low=*Param.C_L_low;
	MatrixXd& C_D_low=*Param.C_D_low;
	MatrixXd& alpha_high=*Param.alpha_high;
	MatrixXd& mach_high=*Param.mach_high;
	MatrixXd& C_L_high=*Param.C_L_high;
	MatrixXd& C_D_high=*Param.C_D_high;
	//cout<<"Phase: "<<iphase<<std::endl;
	if(iphase <= 4){
		mach=Limit(speed/330, 20, 0);			//Mach number
		smooth_bilinear_interpolation(&C_L, alpha_deg, mach, alpha_low, mach_low, C_L_low);
		smooth_bilinear_interpolation(&C_D, alpha_deg, mach, alpha_low, mach_low, C_D_low);
		S=Param.S[0];
	}
	else{
		mach=Limit(speed/330, 23, 3.5);
		smooth_bilinear_interpolation(&C_L, alpha_deg, mach, alpha_high, mach_high, C_L_high);
		smooth_bilinear_interpolation(&C_D, alpha_deg, mach, alpha_high, mach_high, C_D_high);
		S=Param.S[1];
	}
	q=pow(speed, 2)*amb_den(alt)/2;
	L_temp=q*S*C_L;
	D_temp=q*S*C_D;
	L=(iphase<4 || iphase>6)?L_temp:0;
	D=(iphase<4 || iphase>6)?D_temp:0;
	T=(iphase==3 || iphase==4)?Param.T[2]:((iphase==1)?Param.T[0]:((iphase==2)?Param.T[1]:0));
	I_sp=(iphase>2)?Param.I_sp[2]:((iphase==2)?Param.I_sp[1]:((iphase==1)?Param.I_sp[0]:0));
	
	if(iphase==1 || iphase==8)
	{
		omega_1=controls[1];
		ee_1=states[4];					//state[4] to state[7] --> Euler's Parameters					
		ee_2=states[5];
		ee_3=states[6];
		nn=states[7];
		alpha=states[8];				//Angle of Attack
		m=states[9];
		
		derivatives[0]=speed*(1-2*(ee_1*ee_1+ee_2*ee_2));
		derivatives[1]=(2*speed)*(ee_1*ee_2+ee_3*nn)/(r*cos(lat));
		derivatives[2]=(2*speed)*(ee_1*ee_3-ee_2*nn)/r;
		
		term1=(1/m)*(T*cos(alpha)-D);
		term2=(Param.mu_e/pow(r, 2))*(1-2*(ee_2*ee_2+ee_3*ee_3));
		term3=r*pow(Param.omega_e,2)*cos(lat)*(cos(lat)*(1-2*(ee_2*ee_2+ee_3*ee_3)));
		term4=2*sin(lat)*(ee_1*ee_3-ee_2*nn);
		derivatives[3]=term1-term2+term3-term4;
			
		//Calculation of omega_2
		term1=T*sin(alpha)/(m*speed);
		term2=2*((speed/r)-(Param.mu_e/(r*r*speed)))*(ee_1*ee_3+ee_2*nn);
		term3=4*Param.omega_e*(sin(lat)*(ee_1*ee_2-ee_3*nn)+cos(lat)*(ee_2*ee_3+ee_1*nn));
		term4=(2*r*pow(Param.omega_e, 2)*cos(lat)/speed)*(cos(lat)*(ee_1*ee_3+ee_2*nn)-sin(lat)*(0.5-ee_1*ee_1-ee_2*ee_2));
		omega_2=-term1-term2-term3-term4;
		
		//Caclulation of omega_3
		term1=(T*sin(alpha)+L)/(m*speed);
		term2=2*(speed/r-Param.mu_e/(r*r*speed))*(ee_1*ee_2-ee_3*nn);
		term3=4*Param.omega_e*(sin(lat)*(ee_1*ee_3+ee_2*nn)+cos(lat)*(0.5-ee_1*ee_1-ee_2*ee_2));
		term4=(2*r*pow(Param.omega_e, 2)/speed)*cos(lat)*(cos(lat)*(ee_1*ee_2-ee_3*nn)-sin(lat)*(ee_2*ee_3+ee_1*nn));
		omega_3=term1+term2-term3+term4;
		
		derivatives[4]=0.5*(nn*omega_1-ee_3*omega_2+ee_2*omega_3);
		derivatives[5]=0.5*(ee_3*omega_1-nn*omega_2-ee_1*omega_3);
		derivatives[6]=0.5*(-ee_2*omega_1+ee_1*omega_2+nn*omega_3);
		derivatives[7]=-0.5*(ee_1*omega_2+ee_2*omega_2+ee_3*omega_3);
		derivatives[8]=u_alpha;
		
		temp_mass=-T/(I_sp*Param.g);
		derivatives[9]=(iphase==1)?temp_mass:0;
	}
	else
	{
		u_bank=controls[1];
		path_ang=states[4];				//Flight Path Angle
		head_ang=states[5];				//Velocity Heading Angle
		alpha=states[6];				//Angle of Attack
		bank_ang=states[7];				//Banking Angle
		m=states[8];					//Mass of Vehicle
		
		derivatives[0]=speed*sin(path_ang);
		derivatives[1]=speed*cos(path_ang)*sin(head_ang)/(r*cos(lat));
		derivatives[2]=speed*cos(path_ang)*cos(head_ang)/r;
		
		term1=(1/m)*(T*cos(alpha)-D);
		term2=Param.mu_e*sin(path_ang)/pow(r, 2);
		term3=r*pow(Param.omega_e, 2)*cos(lat)*(sin(path_ang)*cos(lat)-cos(path_ang)*sin(lat)*cos(head_ang));
		derivatives[3]=term1-term2+term3;
		
		term1=cos(bank_ang)*(T*sin(alpha)+L)/(m*speed);
		term2=cos(path_ang)*((speed/r)-(Param.mu_e/(pow(r, 2)*speed)));
		term3=2*Param.omega_e*cos(lat)*sin(head_ang);
		term4=(r/speed)*pow(Param.omega_e, 2)*cos(lat)*(cos(path_ang)*cos(lat)+sin(path_ang)*sin(lat)*cos(head_ang));
		derivatives[4]=term1+term2+term3+term4;
		
		term1=(sin(bank_ang)/(m*speed*cos(path_ang)))*(T*sin(alpha)+L);
		term2=(speed/r)*cos(path_ang)*sin(head_ang)*tan(lat);
		term3=2*Param.omega_e*(tan(path_ang)*cos(lat)*cos(head_ang)-sin(lat));
		term4=(r*pow(Param.omega_e, 2)/(speed*cos(path_ang)))*sin(lat)*cos(lat)*sin(head_ang);
		derivatives[5]=term1+term2-term3+term4;
		
		derivatives[6]=u_alpha;			
		derivatives[7]=u_bank;
		
		temp_mass=-T/(I_sp*Param.g);
		derivatives[8]=(iphase>=2 && iphase<=4)?temp_mass:0;
	}
			
	//Path Constraints
	if(iphase>=7){
		acc_sen=(1/(m*Param.g))*sqrt(L*L+D*D);
		path[0]=acc_sen;
		path[1]=q;
		
		Q_dot=(199.87e6)*pow(amb_den(alt)/1.225, 0.5)*pow(speed/7.9053e3, 3.15);
		//path[2]=Q_dot;
	}
	if(iphase==1){
		path[0]=q;
	}
	
	//std::cout<<"In DAE Function"<<'\n'<<std::endl;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

adouble integrand( adouble* states, adouble* controls, adouble* parameters, adouble& time, 
				   adouble* xad, int iphase, Workspace* workspace)
{
	//std::cout<<"In Integrand Function"<<std::endl;
	
	adouble speed=states[3];
	adouble alt=states[0];
	adouble Q_dot;
	
	Q_dot=(199.87e6)*pow(amb_den(alt)/1.225, 0.5)*pow(speed/7.9053e3, 3.15);;
	//std::cout<<"In Integrand Function"<<std::endl;
	return Q_dot;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void events(adouble* e, adouble* initial_states, adouble* final_states,
			adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
			int iphase, Workspace* workspace)
{
	//std::cout<<"In Integrand Cost Function"<<'\n'<<std::endl;
	
	adouble Q;
	
	//Initial States
	adouble alt_i=initial_states[0];				//Altitude
	adouble lon_i=initial_states[1];				//Longitude
	adouble lat_i=initial_states[2];				//Latitude
	adouble	speed_i=initial_states[3];				//Speed
	adouble mass_i, alpha_i, path_ang_i, bank_ang_i, ee_1_i, ee_2_i, ee_3_i, nn_i;
	
	//Final States
	adouble alt_f=final_states[0];					//Altitude
	adouble lon_f=final_states[1];					//Longitude
	adouble lat_f=final_states[2];					//Latitude
	adouble	speed_f=final_states[3];				//Speed
	adouble alpha_f, ee_1_f, nn_f, path_ang_f, bank_ang_f;
	
	if (iphase==1){
		ee_1_i=initial_states[4];					//state[4] to state[7] --> Euler's Parameters					
		ee_2_i=initial_states[5];
		ee_3_i=initial_states[6];
		nn_i=initial_states[7];
		alpha_i=initial_states[8];					//Angle of Attack
		mass_i=initial_states[9];
		
		e[0]=alt_i;
		e[1]=lon_i;
		e[2]=lat_i;
		e[3]=speed_i;
		e[4]=mass_i;
		e[5]=alpha_i;
		e[6]=ee_1_i;
		e[7]=ee_2_i;
		e[8]=ee_3_i;
		e[9]=nn_i;
	}
	else if(iphase==2 || iphase==3)
	{
		mass_i=initial_states[8];
		e[0]=mass_i;
	}
	else if(iphase==5)
	{
		bank_ang_f=final_states[7];
		path_ang_f=final_states[4];
		alpha_f=final_states[6];
		e[0]=alpha_f;
		e[1]=bank_ang_f;
		e[2]=alt_f;
		e[3]=path_ang_f;
	}
	else if(iphase==6)
	{
		bank_ang_i=initial_states[7];
		path_ang_i=initial_states[4];
		alpha_i=initial_states[6];
		e[0]=alpha_i;
		e[1]=bank_ang_i;
		e[2]=alt_i;
		e[3]=path_ang_i;
		e[4]=alt_f;
	}
	else if(iphase>=7)
	{
		Q=integrate(integrand, xad, iphase, workspace);	//Integral Constraint
		if(iphase==7)
		{
			e[0]=Q;
			e[1]=alt_i;
		}
		else if(iphase==8)
		{
			ee_1_i=initial_states[4];					//state[4] to state[7] --> Euler's Parameters					
			ee_2_i=initial_states[5];
			ee_3_i=initial_states[6];
			nn_i=initial_states[7];
			ee_1_f=final_states[4];
			nn_f=final_states[7];
			alpha_f=final_states[8];
			
			e[0]=Q;
			e[1]=alt_f;
			e[2]=lon_f;
			e[3]=lat_f;
			e[4]=speed_f;
			e[5]=alpha_f;
			e[6]=ee_1_f;
			e[7]=nn_f;
			e[8]=pow(ee_1_i, 2)+pow(ee_2_i, 2)+pow(ee_3_i, 2)+pow(nn_i, 2);
		}
	}
	//std::cout<<"In Integrand Cost Function"<<'\n'<<std::endl;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
	//std::cout<<"In Linkages Function"<<std::endl;
	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble temp[10]={0};
	Eigen::Matrix<adouble, 8, 10> xf;
	Eigen::Matrix<adouble, 8, 10> xi;
	adouble tf[8];
	adouble ti[8];
	adouble xf_p1_path, xf_p1_head, xf_p1_bank;
	adouble xi_p8_path, xi_p8_head, xi_p8_bank;
	int phase, state;
	int offset=0;
	
	//Get initial and final states as well as times to define interior point constraints (linkages)
	for(phase=0; phase<8; phase=phase+1)
	{
		get_final_states(temp, xad, phase+1, workspace);
		xf.row(phase)<<temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],
					   temp[7],temp[8],temp[9];
		get_initial_states(temp, xad, phase+1, workspace);
		xi.row(phase)<<temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],
					   temp[7],temp[8],temp[9];
		tf[phase]=get_final_time(xad, phase+1, workspace);
		ti[phase]=get_initial_time(xad, phase+1, workspace);
	}
		
	//Time continuity constraints (3)
	for(phase=0; phase<=6; phase=phase+1)
	{
		linkages[offset]=ti[phase+1]-tf[phase];
		offset++;
	}
	
	//Conversion of Euler parameters to system states for phases 1 and 8 
	xf_p1_path=map_path_ang(xf(0,4), xf(0,5), xf(0,6), xf(0,7));
	xf_p1_head=map_head_ang(xf(0,4), xf(0,5), xf(0,6), xf(0,7));
	xf_p1_bank=map_bank_ang(xf(0,4), xf(0,5), xf(0,6), xf(0,7));
	
	xi_p8_path=map_path_ang(xi(7,4), xi(7,5), xi(7,6), xi(7,7));
	xi_p8_head=map_head_ang(xi(7,4), xi(7,5), xi(7,6), xi(7,7));
	xi_p8_bank=map_bank_ang(xi(7,4), xi(7,5), xi(7,6), xi(7,7));
	
	//State continuity constraints (56)
	for (state=0; state<8; state=state+1)
	{
		linkages[offset++]=xi(2,state)-xf(1,state);
		linkages[offset++]=xi(3,state)-xf(2,state);
		linkages[offset++]=xi(4,state)-xf(3,state);
		linkages[offset++]=xi(5,state)-xf(4,state);
		linkages[offset++]=xi(6,state)-xf(5,state);
		switch (state)
		{
			case 4:
				linkages[offset++]=xi(1,state)-xf_p1_path;
				linkages[offset++]=xi_p8_path-xf(6,state);
				break;
			case 5:
				linkages[offset++]=xi(1,state)-xf_p1_head;
				linkages[offset++]=xi_p8_head-xf(6,state);
				break;
			case 6:
				linkages[offset++]=xi(1,state)-xf(0,8);
				linkages[offset++]=xi(7,8)-xf(6,state);
				break;
			case 7:
				linkages[offset++]=xi(1,state)-xf_p1_bank;
				linkages[offset++]=xi_p8_bank-xf(6,state);
				break;
			default:
				linkages[offset++]=xi(1,state)-xf(0,state);
				linkages[offset++]=xi(7,state)-xf(6,state);
				break;
		}
	}
	
	//Mass transition constraints (1)
	linkages[offset++]=xi(3,8)-xf(2,8);
	
	//std::cout<<"In Linkages Function"<<'\n'<<std::endl;
}
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
int main()
{
	int low_offset=0, up_offset=0;
	int i, j;
	//double EQ_TOL = 0.0002;
	double inf=1e20;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//std::cout<<"Before Level1"<<std::endl;
	//Declare key structures
	Alg algorithm;
	Sol solution;
	Prob problem;
	
	//Register problem name
	problem.name="Ascent-Entry Mission";
	problem.outfilename="Trajectory.txt";
	
	//Instance of Parameters structures
	Constants_ Param;
	problem.user_data=(void*) &Param;
	
	//Problem level 1 setup
	problem.nphases=8;
	problem.nlinkages=64;
	psopt_level1_setup(problem);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//std::cout<<"Before Level2"<<std::endl;
	//Level 2 setup
	for(i=1; i<=8; i++)
	{
		problem.phases(i).nstates=(i==1 || i==8)?10:9;
		problem.phases(i).ncontrols=2;
		problem.phases(i).nevents=0;
		problem.phases(i).npath=0;
		problem.phases(i).nodes = (RowVectorXi(2) << 18, 20).finished();
	}
	
	problem.phases(1).nevents=10;
	problem.phases(2).nevents=1;
	problem.phases(3).nevents=1;
	problem.phases(5).nevents=4;
	problem.phases(6).nevents=5;
	problem.phases(7).nevents=2;
	problem.phases(8).nevents=9;

	
	problem.phases(1).npath=1;
	problem.phases(7).npath=3;
	problem.phases(8).npath=3;
	
	psopt_level2_setup(problem, algorithm);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//std::cout<<"Before Parameters Initialization"<<std::endl;
	//Initialize Parameters
	double init1[8]={0,0,0,0,0,0,11.86,1.86};					//alpha_bar
	double init2[8]={25,25,25,25,0,0,25,25};					//alpha_max
	double init3[8]={10,10,10,10,10,10,10,10};					//u_alpha_max
	double init4[8]={30,30,30,30,30,30,30,30};					//u_bank_max
	MatrixXd Mat1(6, 9), Mat2(6, 9), Mat3(6, 1), Mat4(6, 1); 
	MatrixXd Mat5(9, 1), Mat6(7, 1), Mat7(6, 7), Mat8(6, 7);
	
	double h_atm=80e3;											//Minimum altitude constant
	double acc_sen_max=12;										//Upper limit of sensed acceleration
	double q_max=126.3e3, q_min=12e3;							//Maximum and minimum dynamic pressure
	double Q_dot_max=9e6;										//Maximum stagnation point heating rate								
	double Q_max=3400e6;										//Heating load
	double t_s1=56.4, t_s2=117.1, t_s3=189.1;					//Engine burnout times					
	double t_fair=179.1;										//Fairing separation time
	double m_s2=38780, m_s3=11110;								//Mass of boost vehicle at stage 2 and stage 3 ignition
	double m_fair=400;											//Mass of fairing
	double h_peak_min=100e3; double h_peak_max=200e3;			//Window of allowable peak altitudes
	//Initial Conditions
	double t0=2.52; double alt0=167; double lat0=34.58; double lon0=-120.63; 
	double speed0=40; double ee_10=0; double ee_20=0; double ee_30=0; 
	double nn0=1; double alpha0=0; double m0=85743;
	//Final Conditions
	double altf=0.0; double latf=8.70; double lonf=-172.30; 
	double speedf=1219; double ee_1f=0; double nnf=0; double alphaf=0;
	
	for(i=0; i<8; i++)
	{
		Param.alpha_bar[i]=init1[i];
		Param.alpha_max[i]=init2[i];
		Param.u_alpha_max[i]=init3[i];
		Param.u_bank_max[i]=init4[i];
	}
	Param.omega_e=0.00007292115856; 
	Param.mu_e=398600.4405e9; 
	Param.g=9.8066498;
	Param.Re=6378.166*1000;
	Param.T[0]=2224.1*1000; Param.T[1]=1222.9*1000; Param.T[2]=289.1*1000;
	Param.I_sp[0]=282; Param.I_sp[1]=309; Param.I_sp[2]=300;
	Param.S[0]=4.307; Param.S[1]=0.48387;
	Mat1<< 	0,		0,		0,		0,		0,		0,		0,		0,		0,
			0.2,	0.25,	0.35,	0.4,	0.35,	0.25,	0.25,	0.25,	0.25,
			0.6,	0.6,	0.75,	0.9,	0.8, 	0.7,	0.7,	0.65,	0.65,
			1.1,	1.1,	1.4,	1.65,	1.5,	1.4,	1.25,	1.1,	1,
			1.6,	1.6,	2.1,	2.65,	2.5,	2.1,	1.75,	1.6,	1.5,
			2.2,	2.2,	2.9,	3.9,	3.75,	2.7,	2.3,	2.1,	2;
	Param.C_L_low=&Mat1;
	
	Mat2<< 	0.25,	0.2,	0.4,	0.5,	0.4,	0.2,	0.2,	0.2,	0.2,
			0.25,	0.2, 	0.4, 	0.6, 	0.45, 	0.25, 	0.25,	0.25,	0.25,
			0.35,	0.3,	0.5,	0.7,	0.5,	0.4,	0.4,	0.4,	0.4,
			0.5,	0.5,	0.75,	1,		0.8,	0.8,	0.8,	0.75,	0.75,
			0.8,	0.8,	1.2,	1.5,	1.5,	1.5,	1.3,	1.25,	1.25,
			1.25,	1.25,	1.75,	2.4,	2.5,	2.3,	2.15,	2,		2;
	Param.C_D_low=&Mat2;
	
	Mat3<< 	0,		5,		10,		15,		20,		25;
	Mat4<< 	0,		5,		10,		15,		20,		25;
	Param.alpha_low=&Mat3;
	Param.alpha_high=&Mat4;
	
	Mat5<< 	0,		0.5,	1.0,	1.5,	2,		3,		4,		5,		20;
	Mat6<< 	3.5,	5,		8,		10,		15,		20,		23;	 
	Param.mach_low=&Mat5;
	Param.mach_high=&Mat6;
	
	Mat7<< 	 0,		0,		0,		0,		0,		0,		0,
			 0.24,	0.24,	0.24,	0.2,	0.18,	0.18,	0.18,
			 0.44,	0.42,	0.4,	0.38,	0.38,	0.36,	0.36,
			 0.74,	0.72,	0.7,	0.64,	0.6,	0.58,	0.58,
			 1.04,	1,		0.96,	0.92,	0.88,	0.8,	0.8,
			 1.36,	1.32,	1.24,	1.16,	1.1,	1.02,	1;
	Param.C_L_high=&Mat7;

	Mat8<< 	 0.14,	0.12,	0.08,	0.06,	0.06,	0.06,	0.06,
			 0.16,	0.14,	0.1,	0.08,	0.08,	0.08,	0.08,
			 0.2,	0.18,	0.14,	0.12,	0.12,	0.12,	0.12,
			 0.3,	0.26,	0.22,	0.2,	0.2,	0.2,	0.2,
			 0.48,	0.44,	0.36,	0.32,	0.3,	0.3,	0.3,
			 0.7,	0.64,	0.54,	0.48,	0.46,	0.46,	0.46;					 
	Param.C_D_high=&Mat8;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//std::cout<<"Before Bounds Definition"<<std::endl;
	
	double alt_min=0, alt_max=inf;
	double lon_min=-inf, lon_max=inf;
	double lat_min=-inf, lat_max=inf;
	double speed_min=0, speed_max=inf;
	double ee_min=-inf, ee_max=inf;
	double nn_min=-inf, nn_max=inf;
	double bank_ang_min=-M_PI/2, bank_ang_max=M_PI/2;
	double head_ang_min=-inf, head_ang_max=inf;
	double path_ang_min=-inf, path_ang_max=inf;
	double m_min=907.186, m_max=inf;
	double temp, temp_alpha_max;
	
	//Problem bounds information
	for (i=1; i<=8; i++)
	{
		temp=-1*Param.alpha_max[i-1];
		temp_alpha_max=(i<=4)?temp:0;
		if(i>=2 && i<=7)
		{
			problem.phases(i).bounds.lower.states << alt_min, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, deg_to_rad(temp_alpha_max), 
													 bank_ang_min, m_min;
			problem.phases(i).bounds.upper.states << alt_max, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, deg_to_rad(Param.alpha_max[i-1]),
													 bank_ang_max, m_max;
			if(i<=6 && i>=4){
				problem.phases(i).bounds.lower.states(0) = h_atm;
			}
		}
		else
		{
			problem.phases(i).bounds.lower.states << alt_min, lon_min, lat_min, speed_min, ee_min, ee_min, ee_min, nn_min, deg_to_rad(temp_alpha_max), m_min;
			problem.phases(i).bounds.upper.states << alt_max, lon_max, lat_max, speed_max, ee_max, ee_max, ee_max, nn_max, deg_to_rad(Param.alpha_max[i-1]), m_max;
		}

		//Control bounds
		problem.phases(i).bounds.lower.controls(0) = -1*deg_to_rad(Param.u_alpha_max[i-1]);
		problem.phases(i).bounds.upper.controls(0) = deg_to_rad(Param.u_alpha_max[i-1]);
		problem.phases(i).bounds.lower.controls(1) = -1*deg_to_rad(Param.u_bank_max[i-1]);
		problem.phases(i).bounds.upper.controls(1) = deg_to_rad(Param.u_bank_max[i-1]);
		
	}

	//Path bounds
	problem.phases(1).bounds.upper.path(0) = q_max;
	problem.phases(1).bounds.lower.path(0) = -inf;
	
	problem.phases(7).bounds.upper.path(0) = acc_sen_max;
	problem.phases(7).bounds.lower.path(0) = -inf;
	problem.phases(7).bounds.upper.path(1) = q_min;
	problem.phases(7).bounds.lower.path(1) = -inf;
	problem.phases(7).bounds.upper.path(2) = Q_dot_max;
	problem.phases(7).bounds.lower.path(2) = -inf;
	
	problem.phases(8).bounds.upper.path(0) = acc_sen_max;
	problem.phases(8).bounds.lower.path(0) = -inf;
	problem.phases(8).bounds.upper.path(1) = inf;	
	problem.phases(8).bounds.lower.path(1) = q_min;
	problem.phases(8).bounds.upper.path(2) = Q_dot_max;
	problem.phases(8).bounds.lower.path(2) = -inf;
	
	//Time Linkage Bounds
	problem.phases(1).bounds.lower.StartTime = t0;
	problem.phases(1).bounds.upper.StartTime = t0;

	problem.phases(1).bounds.lower.EndTime = t_s1;
	problem.phases(1).bounds.upper.EndTime = t_s1;	
	//problem.phases(2).bounds.lower.StartTime = t_s1;
	//problem.phases(2).bounds.upper.StartTime = t_s1;

	problem.phases(2).bounds.lower.EndTime = t_s2;
	problem.phases(2).bounds.upper.EndTime = t_s2;
	//problem.phases(3).bounds.lower.StartTime = t_s2;
	//problem.phases(3).bounds.upper.StartTime = t_s2;
	
	problem.phases(3).bounds.lower.EndTime = t_fair;
	problem.phases(3).bounds.upper.EndTime = t_fair;
	//problem.phases(4).bounds.lower.StartTime = t_fair;
	//problem.phases(4).bounds.upper.StartTime = t_fair;
	
	problem.phases(4).bounds.lower.EndTime = t_s3;
	problem.phases(4).bounds.upper.EndTime = t_s3;
	//problem.phases(5).bounds.lower.StartTime = t_s3;
	//problem.phases(5).bounds.upper.StartTime = t_s3;
	
	/* problem.phases(5).bounds.lower.EndTime = t_s3;
	problem.phases(5).bounds.upper.EndTime = inf;
	problem.phases(6).bounds.lower.StartTime = t_s3;
	problem.phases(6).bounds.upper.StartTime = inf;
	
	problem.phases(6).bounds.lower.EndTime = t_s3;
	problem.phases(6).bounds.upper.EndTime = inf;
	problem.phases(7).bounds.lower.StartTime = t_s3;
	problem.phases(7).bounds.upper.StartTime = inf;
	
	problem.phases(7).bounds.lower.EndTime = t_s3;
	problem.phases(7).bounds.upper.EndTime = inf;
	problem.phases(8).bounds.lower.StartTime = t_s3;
	problem.phases(8).bounds.upper.StartTime = inf; */

	for(i=0; i<7; i++)
	{
		problem.bounds.lower.linkage(low_offset++)= 0;
		problem.bounds.upper.linkage(up_offset++)= 0;
	}

	//State Linkage Bounds
	for(i=0; i<56; i++)
	{
		problem.bounds.lower.linkage(low_offset++)= 0;
		problem.bounds.upper.linkage(up_offset++)= 0;
	}
				
	//Mass Transition Bounds
	problem.bounds.lower.linkage(low_offset++)= m_fair;
	problem.bounds.upper.linkage(up_offset++)= m_fair;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	//std::cout<<"Before Events"<<std::endl;
	
	//Events
	problem.phases(1).bounds.lower.events << alt0, deg_to_rad(lon0), deg_to_rad(lat0), speed0, m0, deg_to_rad(alpha0), ee_10, ee_20, ee_30, nn0;
	problem.phases(1).bounds.upper.events << alt0, deg_to_rad(lon0), deg_to_rad(lat0), speed0, m0, deg_to_rad(alpha0), ee_10, ee_20, ee_30, nn0;
	
	problem.phases(2).bounds.lower.events << m_s2;
	problem.phases(2).bounds.upper.events << m_s2;
	
	problem.phases(3).bounds.lower.events << m_s3;
	problem.phases(3).bounds.upper.events << m_s3;
	
	problem.phases(5).bounds.lower.events << 0, 0, h_peak_min, 0; 
	problem.phases(5).bounds.upper.events << 0, 0, h_peak_max, 0; 
	
	problem.phases(6).bounds.lower.events << 0, 0, h_peak_min, 0, h_atm;
	problem.phases(6).bounds.upper.events << 0, 0, h_peak_max, 0, h_atm;
	
	problem.phases(7).bounds.lower.events << 0, h_atm;
	problem.phases(7).bounds.upper.events << Q_max, h_atm;
	
	problem.phases(8).bounds.lower.events << 0, altf, deg_to_rad(lonf), deg_to_rad(latf), speedf, deg_to_rad(alphaf), ee_1f, nnf, 1;
	problem.phases(8).bounds.upper.events << Q_max, altf, deg_to_rad(lonf), deg_to_rad(latf), speedf, deg_to_rad(alphaf), ee_1f, nnf, 1;
											 		
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//std::cout<<"Before Guess"<<std::endl;
	
	//Initial Guess
	for(i=1; i<=8; i++)
	{
		if(i==1 || i==8)
		{
			problem.phases(i).guess.states=zeros(10,12);
			problem.phases(i).guess.states.row(0)=linspace(alt0, h_atm, 12);
			problem.phases(i).guess.states.row(1)=linspace(deg_to_rad(lon0), deg_to_rad(lon0), 12);
			problem.phases(i).guess.states.row(2)=linspace(deg_to_rad(lat0), deg_to_rad(lat0), 12);
			problem.phases(i).guess.states.row(3)=linspace(speed0, speedf, 12);
			problem.phases(i).guess.states.row(4)=linspace(ee_10, ee_1f, 12);
			problem.phases(i).guess.states.row(5)=linspace(ee_20, ee_20, 12);
			problem.phases(i).guess.states.row(6)=linspace(ee_30, ee_30, 12);
			problem.phases(i).guess.states.row(7)=linspace(nn0, nnf, 12);
			problem.phases(i).guess.states.row(8)=linspace(deg_to_rad(alpha0), deg_to_rad(alpha0), 12);
			problem.phases(i).guess.states.row(9)=linspace(m0, m_s2, 12);

		}
		else
		{
			problem.phases(i).guess.states=zeros(9,12);
			problem.phases(i).guess.states.row(0)=linspace(alt0, h_atm, 12);
			problem.phases(i).guess.states.row(1)=linspace(deg_to_rad(lon0), deg_to_rad(lonf), 12);
			problem.phases(i).guess.states.row(2)=linspace(deg_to_rad(lat0), deg_to_rad(latf), 12);
			problem.phases(i).guess.states.row(3)=linspace(speed0, speedf, 12);
			problem.phases(i).guess.states.row(4)=linspace(-M_PI/2, M_PI/2, 12);
			problem.phases(i).guess.states.row(5)=linspace(-M_PI, M_PI, 12);
			problem.phases(i).guess.states.row(6)=linspace(deg_to_rad(alpha0), deg_to_rad(25.0), 12);
			problem.phases(i).guess.states.row(7)=linspace(-M_PI, M_PI, 12);
			problem.phases(i).guess.states.row(8)=linspace(m_s2, m_s3, 12);
		}
		problem.phases(i).guess.controls = zeros(2,12);
		problem.phases(i).guess.controls.row(0) = zeros(1, 12);
		problem.phases(i).guess.controls.row(1) = zeros(1, 12);
	}
	
	//Phase 3
	problem.phases(3).guess.states.row(8)=linspace(m_s3, m_min, 12);
	
	//Phase 4
	problem.phases(4).guess.states.row(0)=linspace(h_atm, h_peak_min, 12);
	problem.phases(4).guess.states.row(8)=linspace(m_s3, m_min, 12);
	
	//Phase 5
	problem.phases(5).guess.states.row(0)=linspace(h_atm, h_peak_max, 12);
	problem.phases(5).guess.states.row(4)=linspace(0, 0, 12);
	problem.phases(5).guess.states.row(6)=linspace(0, 0, 12);
	problem.phases(5).guess.states.row(7)=linspace(0, 0, 12);
	problem.phases(5).guess.states.row(8)=linspace(m_s3, m_min, 12);
	
	//Phase 6
	problem.phases(6).guess.states.row(0)=linspace(h_peak_max, h_atm, 12);
	problem.phases(6).guess.states.row(4)=linspace(0, 0, 12);
	problem.phases(6).guess.states.row(6)=linspace(0, 0, 12);
	problem.phases(6).guess.states.row(7)=linspace(0, 0, 12);
	problem.phases(6).guess.states.row(8)=linspace(m_s3, m_min, 12);
	
	//Phase 7
	problem.phases(7).guess.states.row(0)=linspace(h_atm, altf, 12);
	problem.phases(7).guess.states.row(8)=linspace(m_s3, m_min, 12);

	//Phase 8
	problem.phases(8).guess.states.row(0)=linspace(h_atm, altf, 12);
	problem.phases(8).guess.states.row(1)=linspace(deg_to_rad(lonf), deg_to_rad(lonf), 12);
	problem.phases(8).guess.states.row(2)=linspace(deg_to_rad(latf), deg_to_rad(latf), 12);
	problem.phases(8).guess.states.row(8)=linspace(deg_to_rad(alphaf), deg_to_rad(alphaf), 12);
	problem.phases(8).guess.states.row(9)=linspace(m_s3-m_fair, m_min, 12);

	problem.phases(1).guess.time = linspace(t0,t_s1, 12);
	problem.phases(2).guess.time = linspace(t_s1,t_s2, 12);
	problem.phases(3).guess.time = linspace(t_s2,t_fair, 12);
	problem.phases(4).guess.time = linspace(t_fair,t_s3, 12);
	problem.phases(5).guess.time = linspace(t_s3,t_s3+100, 12);
	problem.phases(6).guess.time = linspace(t_s3+50,t_s3+200, 12);
	problem.phases(7).guess.time = linspace(t_s3+100,t_s3+300, 12);
	problem.phases(8).guess.time = linspace(t_s3+150,t_s3+400, 12);
	
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
	algorithm.derivatives="automatic";
	algorithm.nlp_iter_max=1000;
	algorithm.collocation_method="Chebyshev";
	//algorithm.mesh_refinement="automatic";
	//algorithm.ode_tolerance=1.e-5;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Calling PSOPT to solve problems
	//std::cout<<"Before PSOPT"<<std::endl;
	psopt(solution, problem, algorithm);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Extracting relevant variables
	MatrixXd x_ph1, x_ph2, x_ph3, x_ph4, x_ph5, x_ph6, x_ph7, x_ph8;
	MatrixXd t_ph1, t_ph2, t_ph3, t_ph4, t_ph5, t_ph6, t_ph7, t_ph8;
	MatrixXd u_ph1, u_ph2, u_ph3, u_ph4, u_ph5, u_ph6, u_ph7, u_ph8;
	MatrixXd x_1_8, x_2_7, u, t;
	
	x_ph1 = solution.get_states_in_phase(1);
	x_ph2 = solution.get_states_in_phase(2);
	x_ph3 = solution.get_states_in_phase(3);
	x_ph4 = solution.get_states_in_phase(4);
	x_ph5 = solution.get_states_in_phase(5);
	x_ph6 = solution.get_states_in_phase(6);
	x_ph7 = solution.get_states_in_phase(7);
	x_ph8 = solution.get_states_in_phase(8);
	
	u_ph1 = solution.get_controls_in_phase(1);
	u_ph2 = solution.get_controls_in_phase(2);
	u_ph3 = solution.get_controls_in_phase(3);
	u_ph4 = solution.get_controls_in_phase(4);
	u_ph5 = solution.get_controls_in_phase(5);
	u_ph6 = solution.get_controls_in_phase(6);
	u_ph7 = solution.get_controls_in_phase(7);
	u_ph8 = solution.get_controls_in_phase(8);
	
	t_ph1 = solution.get_time_in_phase(1);
	t_ph2 = solution.get_time_in_phase(2);
	t_ph3 = solution.get_time_in_phase(3);
	t_ph4 = solution.get_time_in_phase(4);
	t_ph5 = solution.get_time_in_phase(5);
	t_ph6 = solution.get_time_in_phase(6);
	t_ph7 = solution.get_time_in_phase(7);
	t_ph8 = solution.get_time_in_phase(8);
	
	x_1_8.resize(10, x_ph1.cols()+ x_ph8.cols() );
	x_2_7.resize(9, x_ph2.cols()+ x_ph3.cols()+ x_ph4.cols()+ x_ph5.cols()+ x_ph6.cols()+ x_ph7.cols() );
	u.resize(2, u_ph1.cols()+ u_ph2.cols()+ u_ph3.cols()+ u_ph4.cols()+ u_ph5.cols()+ u_ph6.cols()+ u_ph7.cols()+ u_ph8.cols() );
	t.resize(1, t_ph1.cols()+ t_ph2.cols()+ t_ph3.cols()+ t_ph4.cols()+ t_ph5.cols()+ t_ph6.cols()+ t_ph7.cols()+ t_ph8.cols() );
	
	x_1_8 << x_ph1, x_ph8;
	x_2_7 << x_ph2, x_ph3, x_ph4, x_ph5, x_ph6, x_ph7;
	u << u_ph1, u_ph2, u_ph3, u_ph4, u_ph5, u_ph6, u_ph7, u_ph8;
	t << t_ph1, t_ph2, t_ph3, t_ph4, t_ph5, t_ph6, t_ph7, t_ph8;
	
	Save(x_1_8,"x1.dat");
	Save(x_2_7,"x2.dat");
	Save(u,"u.dat");
	Save(t,"t.dat");
	
	return 0;
}
