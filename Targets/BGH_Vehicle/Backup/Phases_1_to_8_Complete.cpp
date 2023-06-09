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
	double lon_f, lat_f, m_s5, m_s6;
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
    return (180 * fmod(in, 2 * M_PI) / M_PI);
}
double deg_to_rad(double in) 
{
    return (M_PI * fmod(in, 360) / 180);
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
	adouble path_ang_temp, path_ang;
	path_ang=atan2(0.5-ee_2*ee_2-ee_3*ee_3, sqrt((ee_1*ee_1+nn*nn)*(ee_2*ee_2+ee_3*ee_3)));
	return path_ang;
}
double map_path_ang(double ee_1, double ee_2, double ee_3, double nn)
{
	double path_ang;
	path_ang=atan2(0.5-ee_2*ee_2-ee_3*ee_3, sqrt((ee_1*ee_1+nn*nn)*(ee_2*ee_2+ee_3*ee_3)));
	return path_ang;
}

adouble map_head_ang(adouble ee_1, adouble ee_2, adouble ee_3, adouble nn)
{
	adouble head_ang;
	head_ang=atan2(ee_1*ee_2+ee_3*nn, ee_1*ee_3-ee_2*nn);
	return head_ang;
}
double map_head_ang(double ee_1, double ee_2, double ee_3, double nn)
{
	double head_ang;
	head_ang=atan2(ee_1*ee_2+ee_3*nn, ee_1*ee_3-ee_2*nn);
	return head_ang;
}

adouble map_bank_ang(adouble ee_1, adouble ee_2, adouble ee_3, adouble nn)
{
	adouble bank_ang_temp, bank_ang;
	bank_ang=atan2(-ee_3*ee_1-ee_2*nn, ee_2*ee_1-ee_3*nn);
	return bank_ang;
}
double map_bank_ang(double ee_1, double ee_2, double ee_3, double nn)
{
	double bank_ang;
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

adouble endpoint_cost_1(adouble* initial_states, adouble* final_states,
					  adouble* parameters,adouble& t0, adouble& tf,
					  adouble* xad, int iphase, Workspace* workspace)
{
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble cost, cost1, cost2;
	double k1=0.00, k2=0.00;
	
	adouble final_lat = final_states[2];
	adouble final_lon = final_states[1];

	if(iphase==3)
	{
		cost1=(final_lat-deg_to_rad(Param.lat_f))*(final_lat-deg_to_rad(Param.lat_f));
		cost2=(final_lon-deg_to_rad(Param.lon_f))*(final_lon-deg_to_rad(Param.lon_f));
	}
	else
	{
		cost1=0;
		cost2=0;
	}
	
	cost=k1*cost1+k2*cost2;
	return cost;
}
adouble endpoint_cost_2(adouble* initial_states, adouble* final_states,
					  adouble* parameters,adouble& t0, adouble& tf,
					  adouble* xad, int iphase, Workspace* workspace)
{
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble cost, cost1, cost2;
	double k1=0.00, k2=0.00;
	
	adouble final_lat = final_states[2];
	adouble final_lon = final_states[1];

	if(iphase==3)
	{
		cost1=(final_lat-deg_to_rad(Param.lat_f))*(final_lat-deg_to_rad(Param.lat_f));
		cost2=(final_lon-deg_to_rad(Param.lon_f))*(final_lon-deg_to_rad(Param.lon_f));
	}
	else
	{
		cost1=0;
		cost2=0;
	}
	
	cost=k1*cost1+k2*cost2;
	return cost;
}
adouble endpoint_cost_3(adouble* initial_states, adouble* final_states,
					  adouble* parameters,adouble& t0, adouble& tf,
					  adouble* xad, int iphase, Workspace* workspace)
{
	return 0.0;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

adouble integrand_cost_1(adouble* states, adouble* controls,
					  adouble* parameters, adouble& time, adouble* xad,
					  int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble Cost;								//Integrand cost function
	adouble term1, term2, term3;				//Temporary variables
	adouble alpha;
	adouble u_alpha, omega_1, u_bank;			//Control inputs
	double u_alpha_max, omega_1_max, alpha_max, alpha_bar, u_bank_max;
	
	u_alpha=rad_to_deg(controls[0]);
	u_alpha_max=Param.u_alpha_max[iphase-1];
	term2=u_alpha/u_alpha_max;
	
	if(iphase==1)
	{
		alpha=rad_to_deg(states[8]);
		alpha_bar=Param.alpha_bar[0];
		alpha_max=Param.alpha_max[0];
		term1=(alpha-alpha_bar)/alpha_max;
		
		omega_1=rad_to_deg(controls[1]);
		omega_1_max=Param.u_bank_max[0];
		term3=omega_1/omega_1_max;
	}
	else
	{	
		alpha=rad_to_deg(states[6]);
		alpha_bar=Param.alpha_bar[iphase-1];
		alpha_max=Param.alpha_max[iphase-1];
		term1=(alpha-alpha_bar)/alpha_max;

		u_bank=rad_to_deg(controls[1]);
		u_bank_max=Param.u_bank_max[iphase-1];	
		term3=u_bank/u_bank_max;
	}
	
	Cost=pow(term1, 2)+pow(term2, 2)+pow(term3, 2);
	return Cost;
}
adouble integrand_cost_2(adouble* states, adouble* controls,
					  adouble* parameters, adouble& time, adouble* xad,
					  int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble Cost;								//Integrand cost function
	adouble term1, term2, term3;				//Temporary variables
	adouble alpha;
	adouble u_alpha, omega_1, u_bank;			//Control inputs
	double u_alpha_max, omega_1_max, alpha_max, alpha_bar, u_bank_max;
	
	u_alpha=rad_to_deg(controls[0]);
	u_alpha_max=Param.u_alpha_max[iphase-1];
	term2=u_alpha/u_alpha_max;
	
	if (iphase==1)
	{
		alpha=rad_to_deg(states[6]);
		alpha_bar=Param.alpha_bar[iphase-1];
		alpha_max=Param.alpha_max[iphase-1];
		term1=(alpha-alpha_bar)/alpha_max;
	}
	else
	{
		term1=0;
	}

	u_bank=rad_to_deg(controls[1]);
	u_bank_max=Param.u_bank_max[iphase-1];	
	term3=u_bank/u_bank_max;
	
	Cost=pow(term1, 2)+pow(term2, 2)+pow(term3, 2);
	return Cost;
}
adouble integrand_cost_3(adouble* states, adouble* controls,
					  adouble* parameters, adouble& time, adouble* xad,
					  int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble Cost;								//Integrand cost function
	adouble term1, term2, term3;				//Temporary variables
	adouble alpha;
	adouble u_alpha, omega_1, u_bank;			//Control inputs
	double u_alpha_max, omega_1_max, alpha_max, alpha_bar, u_bank_max;
	
	u_alpha=controls[0];
	u_alpha_max=deg_to_rad(Param.u_alpha_max[iphase-1]);
	term2=u_alpha/u_alpha_max;
	
	if(iphase==3)
	{
		alpha=states[8];
		alpha_bar=deg_to_rad(Param.alpha_bar[5+iphase-1]);
		alpha_max=deg_to_rad(Param.alpha_max[5+iphase-1]);
		term1=(alpha-alpha_bar)/alpha_max;
		
		omega_1=controls[1];
		omega_1_max=deg_to_rad(Param.u_bank_max[5+iphase-1]);
		term3=omega_1/omega_1_max;
	}
	else
	{	
		alpha=states[6];
		alpha_bar=deg_to_rad(Param.alpha_bar[5+iphase-1]);
		alpha_max=deg_to_rad(Param.alpha_max[5+iphase-1]);
		term1=(alpha-alpha_bar)/alpha_max;

		u_bank=controls[1];
		u_bank_max=deg_to_rad(Param.u_bank_max[5+iphase-1]);	
		term3=u_bank/u_bank_max;
	}
	
	Cost=pow(term1, 2)+pow(term2, 2)+pow(term3, 2);
	return Cost;
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void dae_1(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble C_L, C_D;						//Lift and Drag coefficients
	adouble L, D, q;						//Lift, Drag, Dynamic Pressure
	double T, I_sp, S;						//Thrust, Reference Area
	adouble omega_2, omega_3;				//Angular velocity component
	adouble term1, term2, term3, term4;
	adouble acc_sen;						//Sensed acceleration
	adouble Q_dot;							//Stagnation point heating rate
	adouble alpha_deg;						//Angle of Attack in degrees
	adouble mach;
	adouble u_alpha, u_bank, omega_1;		//Control inputs
	adouble u1, u2;							//Variable bounds	
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
	MatrixXd& alpha_low=*Param.alpha_low;
	MatrixXd& mach_low=*Param.mach_low;
	MatrixXd& C_L_low=*Param.C_L_low;
	MatrixXd& C_D_low=*Param.C_D_low;
	MatrixXd& alpha_high=*Param.alpha_high;
	MatrixXd& mach_high=*Param.mach_high;
	MatrixXd& C_L_high=*Param.C_L_high;
	MatrixXd& C_D_high=*Param.C_D_high;

	mach=speed/330;					//Mach number
	if(iphase==1)
	{
		alpha=states[8];				//Angle of Attack
	}
	else
	{
		alpha=states[6];
	}
	alpha_deg=rad_to_deg(alpha);

	smooth_bilinear_interpolation(&C_L, alpha_deg, mach, alpha_low, mach_low, C_L_low);
	smooth_bilinear_interpolation(&C_D, alpha_deg, mach, alpha_low, mach_low, C_D_low);
	
	S=Param.S[0];
	q=pow(speed, 2)*amb_den(alt)/2;
	L=q*S*C_L;
	D=q*S*C_D;
	T=(iphase==1)?Param.T[0]:((iphase==2)?Param.T[1]:Param.T[2]);
	I_sp=(iphase==1)?Param.I_sp[0]:((iphase==2)?Param.I_sp[1]:Param.I_sp[2]);
	
	if(iphase==1)
	{		
		omega_1=controls[1];
		ee_1=states[4];					//state[4] to state[7] --> Euler's Parameters					
		ee_2=states[5];
		ee_3=states[6];
		nn=states[7];
		m=states[9];
		u1=parameters[0];
		u2=parameters[1];
		
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
		derivatives[5]=0.5*(ee_3*omega_1+nn*omega_2-ee_1*omega_3)+u1;
		derivatives[6]=0.5*(-ee_2*omega_1+ee_1*omega_2+nn*omega_3)+u2;
		derivatives[7]=-0.5*(ee_1*omega_2+ee_2*omega_2+ee_3*omega_3);
		derivatives[8]=u_alpha;
		derivatives[9]=-T/(I_sp*Param.g);
		
		path[0]=q;
	}
	else
	{
		u_bank=controls[1];
		path_ang=states[4];				//Flight Path Angle
		head_ang=states[5];				//Velocity Heading Angle
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
		derivatives[8]=-T/(I_sp*Param.g);
	}	
}
void dae_2(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble C_L, C_D;						//Lift and Drag coefficients
	adouble L, D, q;						//Lift, Drag, Dynamic Pressure
	double T, I_sp, S;						//Thrust, Reference Area
	adouble omega_2, omega_3;				//Angular velocity component
	adouble term1, term2, term3, term4;
	adouble acc_sen;						//Sensed acceleration
	adouble Q_dot;							//Stagnation point heating rate
	adouble alpha_deg;						//Angle of Attack in degrees
	adouble mach;
	adouble u_alpha, u_bank, omega_1;		//Control inputs
	adouble u1, u2;							//Variable bounds	
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
	
	alt=states[0];					//Altitude
	lon=states[1];					//Longitude
	lat=states[2];					//Latitude
	speed=states[3];				//Speed
	path_ang=states[4];				//Flight Path Angle
	head_ang=states[5];				//Velocity Heading Angle
	alpha=states[6];
	bank_ang=states[7];				//Banking Angle
	m=(iphase==1)?states[8]:((iphase==2)?Param.m_s5:Param.m_s6);					//Mass of Vehicle
	
	u_alpha=controls[0];
	u_bank=controls[1];
			
	//Interpolation to calculate Lift and Drag coefficients 
	MatrixXd& alpha_low=*Param.alpha_low;
	MatrixXd& mach_low=*Param.mach_low;
	MatrixXd& C_L_low=*Param.C_L_low;
	MatrixXd& C_D_low=*Param.C_D_low;
	MatrixXd& alpha_high=*Param.alpha_high;
	MatrixXd& mach_high=*Param.mach_high;
	MatrixXd& C_L_high=*Param.C_L_high;
	MatrixXd& C_D_high=*Param.C_D_high;

	adouble r=alt+Param.Re;			//Radial distance of vehicle from centre of the Earth
	mach=speed/330;					//Mach number
	alpha_deg=rad_to_deg(alpha);
	
	if(iphase==1 || iphase==2)
	{
		smooth_bilinear_interpolation(&C_L, alpha_deg, mach, alpha_low, mach_low, C_L_low);
		smooth_bilinear_interpolation(&C_D, alpha_deg, mach, alpha_low, mach_low, C_D_low);
		S=Param.S[0];
	}
	else
	{
		smooth_bilinear_interpolation(&C_L, alpha_deg, mach, alpha_high, mach_high, C_L_high);
		smooth_bilinear_interpolation(&C_D, alpha_deg, mach, alpha_high, mach_high, C_D_high);
		S=Param.S[1];
	}
	
	q=pow(speed, 2)*amb_den(alt)/2;
	L=0;
	D=0;
	T=(iphase==1)?Param.T[2]:0;
	I_sp=Param.I_sp[2];
	
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
	
	if(iphase==1)
	{
		derivatives[8]=-T/(I_sp*Param.g);
	}
}
void dae_3(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase, Workspace* workspace)
{	
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble C_L, C_D;						//Lift and Drag coefficients
	adouble L, D, q;						//Lift, Drag, Dynamic Pressure
	double T, I_sp, S;						//Thrust, Reference Area
	adouble omega_2, omega_3;				//Angular velocity component
	adouble term1, term2, term3, term4;
	adouble acc_sen;						//Sensed acceleration
	adouble Q_dot;							//Stagnation point heating rate
	adouble alpha_deg;						//Angle of Attack in degrees
	adouble mach;
	adouble u_alpha, u_bank, omega_1;		//Control inputs
	adouble u1, u2;							//Variable bounds
	adouble L_temp, D_temp;
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
	MatrixXd& alpha_low=*Param.alpha_low;
	MatrixXd& mach_low=*Param.mach_low;
	MatrixXd& C_L_low=*Param.C_L_low;
	MatrixXd& C_D_low=*Param.C_D_low;
	MatrixXd& alpha_high=*Param.alpha_high;
	MatrixXd& mach_high=*Param.mach_high;
	MatrixXd& C_L_high=*Param.C_L_high;
	MatrixXd& C_D_high=*Param.C_D_high;

	mach=speed/330;					//Mach number
	if(iphase==3)
	{
		alpha=states[8];			//Angle of Attack
	}
	else
	{
		alpha=states[6];
	}
	alpha_deg=rad_to_deg(alpha);

	smooth_bilinear_interpolation(&C_L, alpha_deg, mach, alpha_high, mach_high, C_L_high);
	smooth_bilinear_interpolation(&C_D, alpha_deg, mach, alpha_high, mach_high, C_D_high);
	
	S=Param.S[1];
	q=pow(speed, 2)*amb_den(alt)/2;
	L_temp=q*S*C_L;
	D_temp=q*S*C_D;
	L=(iphase==1)?0:L_temp;
	D=(iphase==1)?0:D_temp;
	T=0;
	I_sp=Param.I_sp[2];
	m=Param.m_s5;
	
	if(iphase==3)
	{		
		omega_1=controls[1];
		ee_1=states[4];					//state[4] to state[7] --> Euler's Parameters					
		ee_2=states[5];
		ee_3=states[6];
		nn=states[7];
		u1=parameters[0];
		u2=parameters[1];
		
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
		derivatives[5]=0.5*(ee_3*omega_1+nn*omega_2-ee_1*omega_3)+u1;
		derivatives[6]=0.5*(-ee_2*omega_1+ee_1*omega_2+nn*omega_3)+u2;
		derivatives[7]=-0.5*(ee_1*omega_2+ee_2*omega_2+ee_3*omega_3);
		derivatives[8]=u_alpha;
	}
	else
	{
		u_bank=controls[1];
		path_ang=states[4];				//Flight Path Angle
		head_ang=states[5];				//Velocity Heading Angle
		bank_ang=states[7];				//Banking Angle
		
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
	}	
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void events_1(adouble* e, adouble* initial_states, adouble* final_states,
			adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
			int iphase, Workspace* workspace)
{
	//Initial States
	adouble alt_i=initial_states[0];				//Altitude
	adouble lon_i=initial_states[1];				//Longitude
	adouble lat_i=initial_states[2];				//Latitude
	adouble	speed_i=initial_states[3];				//Speed
	adouble mass_i, alpha_i, ee_1_i, ee_2_i, ee_3_i, nn_i;
	adouble alt_f;
	
	if(iphase==1)
	{
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
		e[4]=ee_1_i;
		e[5]=ee_2_i;
		e[6]=ee_3_i;
		e[7]=nn_i;
		e[8]=alpha_i;
		e[9]=mass_i;
		
	}
	else if (iphase==2 )
	{
		mass_i=initial_states[8];
		e[0]=mass_i;
	}
	else
	{
		mass_i=initial_states[8];
		alt_f=final_states[0];
		e[0]=mass_i;
		e[1]=alt_f;
	}
}
void events_2(adouble* e, adouble* initial_states, adouble* final_states,
			adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
			int iphase, Workspace* workspace)
{
	//Initial States
	adouble alt_i=initial_states[0];				//Altitude
	adouble lon_i=initial_states[1];				//Longitude
	adouble lat_i=initial_states[2];				//Latitude
	adouble	speed_i=initial_states[3];				//Speed
	adouble mass_i, alpha_i, path_ang_i, bank_ang_i, head_ang_i;
	adouble alt_f, path_ang_f, alpha_f, bank_ang_f;
	
	if(iphase==1)
	{
		path_ang_i=initial_states[4];
		head_ang_i=initial_states[5];
		alpha_i=initial_states[6];					//Angle of Attack
		bank_ang_i=initial_states[7];
		mass_i=initial_states[8];
		
		e[0]=alt_i;
		e[1]=lon_i;
		e[2]=lat_i;
		e[3]=speed_i;
		e[4]=path_ang_i;
		e[5]=head_ang_i;
		e[6]=alpha_i;
		e[7]=bank_ang_i;
		e[8]=mass_i;
	}
	else if (iphase==2) 
	{
		alt_f=final_states[0];
		path_ang_f=final_states[4];
		alpha_f=final_states[6];
		bank_ang_f=final_states[7];
	
		e[0]=alt_f;
		e[1]=path_ang_f;
		e[2]=alpha_f;
		e[3]=bank_ang_f;
	}
	else
	{
		path_ang_i=initial_states[4];
		head_ang_i=initial_states[5];
		alpha_i=initial_states[6];					//Angle of Attack
		bank_ang_i=initial_states[7];
		alt_f=final_states[0];
		
		e[0]=alt_i;
		e[1]=path_ang_i;
		e[2]=alpha_i;
		e[3]=bank_ang_i;
		e[4]=alt_f;
	}		
}
void events_3(adouble* e, adouble* initial_states, adouble* final_states,
			adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
			int iphase, Workspace* workspace)
{
	//Initial States
	adouble alt_i=initial_states[0];				//Altitude
	adouble lon_i=initial_states[1];				//Longitude
	adouble lat_i=initial_states[2];				//Latitude
	adouble	speed_i=initial_states[3];				//Speed
	adouble alpha_i, path_ang_i, bank_ang_i, head_ang_i;
	adouble ee_1_i, ee_2_i, ee_3_i, nn_i;
	adouble alt_f, lon_f, lat_f, speed_f, ee_1_f, ee_2_f, ee_3_f, nn_f, alpha_f;
	
	if(iphase==1)
	{
		path_ang_i=initial_states[4];
		head_ang_i=initial_states[5];
		alpha_i=initial_states[6];					//Angle of Attack
		bank_ang_i=initial_states[7];
		alt_f=final_states[0];
		
		e[0]=alt_i;
		e[1]=lon_i;
		e[2]=lat_i;
		e[3]=speed_i;
		e[4]=path_ang_i;
		e[5]=head_ang_i;
		e[6]=alpha_i;
		e[7]=bank_ang_i;
		e[8]=alt_f;
	}
	else if (iphase==2)
	{
		alt_i=initial_states[0];
		e[0]=alt_i;
	}
	else if(iphase==3)
	{
		ee_1_i=initial_states[4];
		ee_2_i=initial_states[5];
		ee_3_i=initial_states[6];
		nn_i=initial_states[7];
		
		alt_f=final_states[0];
		lon_f=final_states[1];
		lat_f=final_states[2];
		speed_f=final_states[3];
		ee_1_f=initial_states[4];					//state[4] to state[7] --> Euler's Parameters					
		nn_f=initial_states[7];
		alpha_f=initial_states[8];					//Angle of Attack
		
		e[0]=alt_f;
		e[1]=lon_f;
		e[2]=lat_f;
		e[3]=speed_f;
		e[4]=ee_1_f;
		e[5]=nn_f;
		e[6]=alpha_f;
		e[7]=sqrt(pow(ee_1_i, 2)+pow(ee_2_i, 2)+pow(ee_3_i, 2)+pow(nn_i, 2));
	}
}

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/

void linkages_1( adouble* linkages, adouble* xad, Workspace* workspace)
{
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble temp[10]={0};
	Eigen::Matrix<adouble, 3, 10> xf;
	Eigen::Matrix<adouble, 3, 10> xi;
	adouble tf[3];
	adouble ti[3];
	adouble xf_p1_path, xf_p1_head, xf_p1_bank;
	int phase, state;
	int offset=0;
	
	//Get initial and final states as well as times to define interior point constraints (linkages)
	for(phase=0; phase<3; phase=phase+1)
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
		
	//Time continuity constraints 
	for(phase=0; phase<2; phase=phase+1)
	{
		linkages[offset]=tf[phase]-ti[phase+1];
		offset++;
	} 
	
	//Conversion of Euler parameters to system states for phases 1 and 8 
	xf_p1_path=map_path_ang(xf(0,4), xf(0,5), xf(0,6), xf(0,7));
	xf_p1_head=map_head_ang(xf(0,4), xf(0,5), xf(0,6), xf(0,7));
	xf_p1_bank=map_bank_ang(xf(0,4), xf(0,5), xf(0,6), xf(0,7));
	
	//State continuity constraints
	for (state=0; state<8; state=state+1)
	{
		linkages[offset++]=xf(1,state)-xi(2,state);
		switch (state)
		{
			case 4:
				linkages[offset++]=xf_p1_path-xi(1,state);
				break;
			case 5:
				linkages[offset++]=xf_p1_head-xi(1,state);
				break;
			case 6:
				linkages[offset++]=xf(0,8)-xi(1,state);
				break;
			case 7:
				linkages[offset++]=xf_p1_bank-xi(1,state);
				break;
			default:
				linkages[offset++]=xf(0,state)-xi(1,state);
				break;
		}
	}
}
void linkages_2( adouble* linkages, adouble* xad, Workspace* workspace)
{
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble temp[9]={0};
	Eigen::Matrix<adouble, 3, 9> xf;
	Eigen::Matrix<adouble, 3, 9> xi;
	adouble tf[3];
	adouble ti[3];
	int phase, state;
	int offset=0;
	
	//Get initial and final states as well as times to define interior point constraints (linkages)
	for(phase=0; phase<3; phase=phase+1)
	{
		get_final_states(temp, xad, phase+1, workspace);
		xf.row(phase)<<temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],
					   temp[7],temp[8];
		get_initial_states(temp, xad, phase+1, workspace);
		xi.row(phase)<<temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],
					   temp[7],temp[8];
		tf[phase]=get_final_time(xad, phase+1, workspace);
		ti[phase]=get_initial_time(xad, phase+1, workspace);
	}
		
	//Time continuity constraints 
	for(phase=0; phase<2; phase=phase+1)
	{
		linkages[offset]=tf[phase]-ti[phase+1];
		offset++;
	} 
	
	//State continuity constraints
	for (state=0; state<8; state=state+1)
	{
		linkages[offset++]=xf(0,state)-xi(1,state);
		linkages[offset++]=xf(1,state)-xi(2,state);
	}
}
void linkages_3( adouble* linkages, adouble* xad, Workspace* workspace)
{
	Constants_& Param = *( (Constants_*) workspace->problem->user_data);
	adouble temp[10]={0};
	Eigen::Matrix<adouble, 3, 9> xf;
	Eigen::Matrix<adouble, 3, 9> xi;
	adouble tf[3];
	adouble ti[3];
	adouble xi_p8_path, xi_p8_head, xi_p8_bank;
	int phase, state;
	int offset=0;
	
	//Get initial and final states as well as times to define interior point constraints (linkages)
	for(phase=0; phase<3; phase=phase+1)
	{
		get_final_states(temp, xad, phase+1, workspace);
		xf.row(phase)<<temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],
					   temp[7],temp[8];
		get_initial_states(temp, xad, phase+1, workspace);
		xi.row(phase)<<temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6],
					   temp[7],temp[8];
		tf[phase]=get_final_time(xad, phase+1, workspace);
		ti[phase]=get_initial_time(xad, phase+1, workspace);
	}
		
	//Time continuity constraints 
	for(phase=0; phase<2; phase=phase+1)
	{
		linkages[offset]=tf[phase]-ti[phase+1];
		offset++;
	} 
	
	//Conversion of Euler parameters to system states for phases 1 and 8 
	xi_p8_path=map_path_ang(xi(2,4), xi(2,5), xi(2,6), xi(2,7));
	xi_p8_head=map_head_ang(xi(2,4), xi(2,5), xi(2,6), xi(2,7));
	xi_p8_bank=map_bank_ang(xi(2,4), xi(2,5), xi(2,6), xi(2,7));
	
	//State continuity constraints
	for (state=0; state<8; state=state+1)
	{
		linkages[offset++]=xf(0,state)-xi(1,state);
		switch (state)
		{
			case 4:
				linkages[offset++]=xf(1,state)-xi_p8_path;
				break;
			case 5:
				linkages[offset++]=xf(1,state)-xi_p8_head;
				break;
			case 6:
				linkages[offset++]=xf(1,state)-xi(2,8);
				break;
			case 7:
				linkages[offset++]=xf(1,state)-xi_p8_bank;
				break;
			default:
				linkages[offset++]=xf(1,state)-xi(2,state);
				break;
		}
	}
}
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	

int main()
{
	int low_offset=0, up_offset=0;
	int i, j;
	double TOL = 0.02;
	double mdot, temp_m;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//Declare key structures
	Alg algorithm1, algorithm2, algorithm3;
	Sol solution1, solution2, solution3;
	Prob problem1, problem2, problem3;
	
	//Register problem name
	problem1.name="Ascent-Entry Mission (Phase 1-3)";
	problem1.outfilename="Trajectory_1.txt";
	
	problem2.name="Ascent-Entry Mission (Phase 4-5)";
	problem2.outfilename="Trajectory_2.txt ";
	
	problem3.name="Ascent-Entry Mission (Phase 6-8)";
	problem3.outfilename="Trajectory_3.txt ";
	
	//Instance of Parameters structures
	Constants_ Param;
	problem1.user_data=(void*) &Param;
	problem2.user_data=(void*) &Param;
	problem3.user_data=(void*) &Param;
	
	//Problem level 1 setup
	problem1.nphases=3;
	problem1.nlinkages=2+16;
	psopt_level1_setup(problem1);
	
	problem2.nphases=3;
	problem2.nlinkages=2+16;
	psopt_level1_setup(problem2);
	
	problem3.nphases=3;
	problem3.nlinkages=2+16;
	psopt_level1_setup(problem3);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	

	// Problem1 Level 2 setup
	problem1.phases(1).nstates=10;
	problem1.phases(1).ncontrols=2;
	problem1.phases(1).nevents=10;
	problem1.phases(1).npath=1;
	problem1.phases(1).nparameters=2;	
	//problem1.phases(1).nodes = (RowVectorXi(12) << 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28).finished();
	problem1.phases(1).nodes = (RowVectorXi(12) << 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28).finished();

	problem1.phases(2).nstates=9;
	problem1.phases(2).ncontrols=2;
	problem1.phases(2).nevents=1;
	problem1.phases(2).npath=0;
	problem1.phases(2).nodes = (RowVectorXi(12) << 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26).finished();
	
	problem1.phases(3).nstates=9;
	problem1.phases(3).ncontrols=2;
	problem1.phases(3).nevents=2;
	problem1.phases(3).npath=0;
	problem1.phases(3).nodes = (RowVectorXi(12) << 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26).finished();
	
	psopt_level2_setup(problem1, algorithm1);
	
	// Problem2 Level 2 setup
	problem2.phases(1).nstates=9;
	problem2.phases(1).ncontrols=2;
	problem2.phases(1).nevents=9;
	problem2.phases(1).npath=0;
	problem2.phases(1).nodes = (RowVectorXi(4) << 4, 5, 6, 8).finished();

	problem2.phases(2).nstates=8;
	problem2.phases(2).ncontrols=2;
	problem2.phases(2).nevents=4;
	problem2.phases(2).npath=0;
	problem2.phases(2).nodes = (RowVectorXi(4) << 12, 16, 20, 24).finished();
	
	problem2.phases(3).nstates=8;
	problem2.phases(3).ncontrols=2;
	problem2.phases(3).nevents=5;
	problem2.phases(3).npath=0;
	problem2.phases(3).nodes = (RowVectorXi(4) << 12, 16, 20, 24).finished();
	
	psopt_level2_setup(problem2, algorithm2);
	
	// Problem3 Level 2 setup
	problem3.phases(1).nstates=8;
	problem3.phases(1).ncontrols=2;
	problem3.phases(1).nevents=9;
	problem3.phases(1).npath=0;
	problem3.phases(1).nodes = (RowVectorXi(4) << 18, 20, 26, 32).finished();

	problem3.phases(2).nstates=8;
	problem3.phases(2).ncontrols=2;
	problem3.phases(2).nevents=1;
	problem3.phases(2).npath=0;
	problem3.phases(2).nodes = (RowVectorXi(4) << 18, 20, 26, 32).finished();
	
	problem3.phases(3).nstates=9;
	problem3.phases(3).ncontrols=2;
	problem3.phases(3).nevents=8;
	problem3.phases(3).npath=0;
	problem3.phases(3).nodes = (RowVectorXi(4) << 18, 20, 26, 32).finished();
	
	psopt_level2_setup(problem3, algorithm3);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//Initialize Parameters
	double init1[8]={0,0,0,0,0,0,11.86,11.86};					//alpha_bar
	double init2[8]={25,25,25,25,25,25,25,25};					//alpha_max
	double init3[8]={10,10,10,10,10,10,10,10};					//u_alpha_max
	double init4[8]={30,30,30,30,30,30,30,30};					//u_bank_max
	MatrixXd Mat1(11, 9), Mat2(11, 9), Mat3(11, 1), Mat4(11, 1); 
	MatrixXd Mat5(9, 1), Mat6(7, 1), Mat7(11, 7), Mat8(11, 7);
	
	double h_atm=80e3;											//Minimum altitude constant
	double acc_sen_max=12;										//Upper limit of sensed acceleration
	double q_max=126.3e3, q_min=12e3;							//Maximum and minimum dynamic pressure
	double Q_dot_max=9e6;										//Maximum stagnation point heating rate								
	double Q_max=3400e6;										//Heating load
	double t_s1=56.4, t_s2=117.1, t_s3=189.1;					//Engine burnout times					
	double t_fair=179.1;										//Fairing separation time
	double m_s2=38780, m_s3=11110, m_s4=4600;					//Mass of boost vehicle at stage 2 and stage 3 ignition
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
	Param.lon_f=lonf;
	Param.lat_f=latf;
	Param.m_s5=4000;
	Param.m_s6=907.186;
	
	Mat1<< 	-2.2,	-2.2,	-2.9,	-3.9,	-3.75,	-2.7,	-2.3,	-2.1,	-2,
			-1.6,	-1.6,	-2.1,	-2.65,	-2.5,	-2.1,	-1.75,	-1.6,	-1.5,
			-1.1,	-1.1,	-1.4,	-1.65,	-1.5,	-1.4,	-1.25,	-1.1,	-1,
			-0.6,	-0.6,	-0.75,	-0.9,	-0.8, 	-0.7,	-0.7,	-0.65,	-0.65,
			-0.2,	-0.25,	-0.35,	-0.4,	-0.35,	-0.25,	-0.25,	-0.25,	-0.25,
			0,		0,		0,		0,		0,		0,		0,		0,		0,
			0.2,	0.25,	0.35,	0.4,	0.35,	0.25,	0.25,	0.25,	0.25,
			0.6,	0.6,	0.75,	0.9,	0.8, 	0.7,	0.7,	0.65,	0.65,
			1.1,	1.1,	1.4,	1.65,	1.5,	1.4,	1.25,	1.1,	1,
			1.6,	1.6,	2.1,	2.65,	2.5,	2.1,	1.75,	1.6,	1.5,
			2.2,	2.2,	2.9,	3.9,	3.75,	2.7,	2.3,	2.1,	2;
	Param.C_L_low=&Mat1;
	
	Mat2<< 	1.25,	1.25,	1.75,	2.4,	2.5,	2.3,	2.15,	2,		2,
			0.8,	0.8,	1.2,	1.5,	1.5,	1.5,	1.3,	1.25,	1.25,
			0.5,	0.5,	0.75,	1,		0.8,	0.8,	0.8,	0.75,	0.75,
			0.35,	0.3,	0.5,	0.7,	0.5,	0.4,	0.4,	0.4,	0.4,
			0.25,	0.2, 	0.4, 	0.6, 	0.45, 	0.25, 	0.25,	0.25,	0.25,
			0.25,	0.2,	0.4,	0.5,	0.4,	0.2,	0.2,	0.2,	0.2,
			0.25,	0.2, 	0.4, 	0.6, 	0.45, 	0.25, 	0.25,	0.25,	0.25,
			0.35,	0.3,	0.5,	0.7,	0.5,	0.4,	0.4,	0.4,	0.4,
			0.5,	0.5,	0.75,	1,		0.8,	0.8,	0.8,	0.75,	0.75,
			0.8,	0.8,	1.2,	1.5,	1.5,	1.5,	1.3,	1.25,	1.25,
			1.25,	1.25,	1.75,	2.4,	2.5,	2.3,	2.15,	2,		2;
	Param.C_D_low=&Mat2;
	
	Mat3<< 	-25, 	-20, 	-15, 	-10, 	-5, 	0,		5,		10,		15,		20,		25;
	Mat4<< 	-25, 	-20, 	-15, 	-10, 	-5, 	0,		5,		10,		15,		20,		25;
	Param.alpha_low=&Mat3;
	Param.alpha_high=&Mat4;
	
	Mat5<< 	0,		0.5,	1.0,	1.5,	2,		3,		4,		5,		20;
	Mat6<< 	3.5,	5,		8,		10,		15,		20,		23;	 
	Param.mach_low=&Mat5;
	Param.mach_high=&Mat6;
	
	Mat7<< 	-1.36,	-1.32,	-1.24,	-1.16,	-1.1,	-1.02,	-1,
			-1.04,	-1,		-0.96,	-0.92,	-0.88,	-0.8,	-0.8,
			-0.74,	-0.72,	-0.7,	-0.64,	-0.6,	-0.58,	-0.58,
			-0.44,	-0.42,	-0.4,	-0.38,	-0.38,	-0.36,	-0.36,
			-0.24,	-0.24,	-0.24,	-0.2,	-0.18,	-0.18,	-0.18,
			0,		0,		0,		0,		0,		0,		0,
			0.24,	0.24,	0.24,	0.2,	0.18,	0.18,	0.18,
			0.44,	0.42,	0.4,	0.38,	0.38,	0.36,	0.36,
			0.74,	0.72,	0.7,	0.64,	0.6,	0.58,	0.58,
			1.04,	1,		0.96,	0.92,	0.88,	0.8,	0.8,
			1.36,	1.32,	1.24,	1.16,	1.1,	1.02,	1;
	Param.C_L_high=&Mat7;

	Mat8<< 	0.7,	0.64,	0.54,	0.48,	0.46,	0.46,	0.46,
			0.48,	0.44,	0.36,	0.32,	0.3,	0.3,	0.3,
			0.3,	0.26,	0.22,	0.2,	0.2,	0.2,	0.2,
			0.2,	0.18,	0.14,	0.12,	0.12,	0.12,	0.12,
			0.16,	0.14,	0.1,	0.08,	0.08,	0.08,	0.08,
			0.14,	0.12,	0.08,	0.06,	0.06,	0.06,	0.06,
			0.16,	0.14,	0.1,	0.08,	0.08,	0.08,	0.08,
			0.2,	0.18,	0.14,	0.12,	0.12,	0.12,	0.12,
			0.3,	0.26,	0.22,	0.2,	0.2,	0.2,	0.2,
			0.48,	0.44,	0.36,	0.32,	0.3,	0.3,	0.3,
			0.7,	0.64,	0.54,	0.48,	0.46,	0.46,	0.46;					 
	Param.C_D_high=&Mat8;
	
	double alt_min=0, alt_max=h_peak_max;
	double lon_min=-2*M_PI, lon_max=2*M_PI;
	double lat_min=-(M_PI/2)+TOL, lat_max=(M_PI/2)-TOL;
	double speed_min=TOL, speed_max=6600;
	double ee_min=-1, ee_max=1;
	double nn_min=-1, nn_max=1;
	double bank_ang_min=-M_PI/3, bank_ang_max=M_PI/3;
	double head_ang_min=-M_PI, head_ang_max=M_PI;
	double path_ang_min=-(M_PI/2)+TOL, path_ang_max=(M_PI/2)-TOL;
	double m_min=907.186, m_max=m0;
	double t_max=8000;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Problem bounds information for Phase 1
	problem1.phases(1).bounds.lower.states << alt_min, lon_min, lat_min, speed_min, ee_min, ee_min, ee_min, nn_min, -deg_to_rad(Param.alpha_max[0]), m_min;
	problem1.phases(1).bounds.upper.states << h_atm, lon_max, lat_max, speed_max, ee_max, ee_max, ee_max, nn_max, deg_to_rad(Param.alpha_max[0]), m_max;
	problem1.phases(1).bounds.lower.controls << -1*deg_to_rad(Param.u_alpha_max[0]), -1*deg_to_rad(Param.u_bank_max[0]);
	problem1.phases(1).bounds.upper.controls << deg_to_rad(Param.u_alpha_max[0]), deg_to_rad(Param.u_bank_max[0]);
	problem1.phases(1).bounds.lower.parameters << -1e-5, -1e-5;
	problem1.phases(1).bounds.upper.parameters << 1e-5, 1e-5;
	problem1.phases(1).bounds.upper.path << q_max;
	problem1.phases(1).bounds.lower.StartTime = t0;
	problem1.phases(1).bounds.upper.StartTime = t0;
	problem1.phases(1).bounds.lower.EndTime = t_s1;
	problem1.phases(1).bounds.upper.EndTime = t_s1;

	//Problem bounds information for Phase 2
	problem1.phases(2).bounds.lower.states << alt_min, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, -deg_to_rad(Param.alpha_max[1]), bank_ang_min, m_min;
	problem1.phases(2).bounds.upper.states << h_atm, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, deg_to_rad(Param.alpha_max[1]), bank_ang_max, m_s2;
	problem1.phases(2).bounds.lower.controls << -1*deg_to_rad(Param.u_alpha_max[1]), -1*deg_to_rad(Param.u_bank_max[1]);
	problem1.phases(2).bounds.upper.controls << deg_to_rad(Param.u_alpha_max[1]), deg_to_rad(Param.u_bank_max[1]);
	problem1.phases(2).bounds.lower.StartTime = t_s1;
	problem1.phases(2).bounds.upper.StartTime = t_s1;
	problem1.phases(2).bounds.lower.EndTime = t_s2;
	problem1.phases(2).bounds.upper.EndTime = t_s2;

	//Problem bounds information for Phase 3
	problem1.phases(3).bounds.lower.states << alt_min, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, -deg_to_rad(Param.alpha_max[2]), bank_ang_min, m_min;
	problem1.phases(3).bounds.upper.states << h_atm, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, deg_to_rad(Param.alpha_max[2]), bank_ang_max, m_s3;
	problem1.phases(3).bounds.lower.controls << -1*deg_to_rad(Param.u_alpha_max[2]), -1*deg_to_rad(Param.u_bank_max[2]);
	problem1.phases(3).bounds.upper.controls << deg_to_rad(Param.u_alpha_max[2]), deg_to_rad(Param.u_bank_max[2]);
	problem1.phases(3).bounds.lower.StartTime = t_s2;
	problem1.phases(3).bounds.upper.StartTime = t_s2;
	problem1.phases(3).bounds.lower.EndTime = t_fair;
	problem1.phases(3).bounds.upper.EndTime = t_fair;

	for(i=0; i<18; i++)
	{
		problem1.bounds.lower.linkage(low_offset++)= 0;
		problem1.bounds.upper.linkage(up_offset++)= 0;
	}
	low_offset=0;
	up_offset=0;
	
	//Problem bounds information for Phase 4
	problem2.phases(1).bounds.lower.states << h_atm, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, -deg_to_rad(Param.alpha_max[3]), bank_ang_min, m_min;
	problem2.phases(1).bounds.upper.states << h_peak_max, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, deg_to_rad(Param.alpha_max[3]), bank_ang_max, m_s4;
	problem2.phases(1).bounds.lower.controls << -1*deg_to_rad(Param.u_alpha_max[3]), -1*deg_to_rad(Param.u_bank_max[3]);
	problem2.phases(1).bounds.upper.controls << deg_to_rad(Param.u_alpha_max[3]), deg_to_rad(Param.u_bank_max[3]);
	problem2.phases(1).bounds.lower.StartTime = t_fair;
	problem2.phases(1).bounds.upper.StartTime = t_fair;
	problem2.phases(1).bounds.lower.EndTime = t_s3;
	problem2.phases(1).bounds.upper.EndTime = t_s3;

	//Problem bounds information for Phase 5
	problem2.phases(2).bounds.lower.states << h_atm, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, 0, bank_ang_min;
	problem2.phases(2).bounds.upper.states << h_peak_max, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, deg_to_rad(Param.alpha_max[4]), bank_ang_max;
	problem2.phases(2).bounds.lower.controls << -1*deg_to_rad(Param.u_alpha_max[4]), -1*deg_to_rad(Param.u_bank_max[4]);
	problem2.phases(2).bounds.upper.controls << deg_to_rad(Param.u_alpha_max[4]), deg_to_rad(Param.u_bank_max[4]);
	problem2.phases(2).bounds.lower.StartTime = t_s3;
	problem2.phases(2).bounds.upper.StartTime = t_s3;
	problem2.phases(2).bounds.lower.EndTime = t_s3;
	problem2.phases(2).bounds.upper.EndTime = 1000;
	
	//Problem bounds information for Phase 6
	problem2.phases(3).bounds.lower.states << h_atm, lon_min, lat_min, speed_min, path_ang_min, head_ang_min, 0, bank_ang_min;
	problem2.phases(3).bounds.upper.states << h_peak_max, lon_max, lat_max, speed_max, path_ang_max, head_ang_max, deg_to_rad(Param.alpha_max[5]), bank_ang_max;
	problem2.phases(3).bounds.lower.controls << -1*deg_to_rad(Param.u_alpha_max[5]), -1*deg_to_rad(Param.u_bank_max[5]);
	problem2.phases(3).bounds.upper.controls << deg_to_rad(Param.u_alpha_max[5]), deg_to_rad(Param.u_bank_max[5]);
	problem2.phases(3).bounds.lower.StartTime = t_s3;
	problem2.phases(3).bounds.upper.StartTime = 1000;
	problem2.phases(3).bounds.lower.EndTime = t_s3;
	problem2.phases(3).bounds.upper.EndTime = 1000;

	for(i=0; i<18; i++)
	{
		problem1.bounds.lower.linkage(low_offset++)= 0;
		problem1.bounds.upper.linkage(up_offset++)= 0;
	}
	low_offset=0;
	up_offset=0;
		
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/
	
	//Events Problem1
	problem1.phases(1).bounds.lower.events << alt0, deg_to_rad(lon0), deg_to_rad(lat0), speed0, ee_10, ee_20, ee_30, nn0, deg_to_rad(alpha0), m0;
	problem1.phases(1).bounds.upper.events << alt0, deg_to_rad(lon0), deg_to_rad(lat0), speed0, ee_10, ee_20, ee_30, nn0, deg_to_rad(alpha0), m0;
	problem1.phases(1).bounds.lower.path << 0;
	problem1.phases(1).bounds.upper.path << q_max;
	
	problem1.phases(2).bounds.lower.events << m_s2;
	problem1.phases(2).bounds.upper.events << m_s2;
											 		
	problem1.phases(3).bounds.lower.events << m_s3, h_atm;
	problem1.phases(3).bounds.upper.events << m_s3, h_atm;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Initial Guess Phase 1
	mdot=m0-Param.T[0]*(t_s1-t0)/(Param.I_sp[0]*Param.g);
	problem1.phases(1).guess.states=zeros(10,5);
	problem1.phases(1).guess.states.row(0)=linspace(alt0, alt0, 5);
	problem1.phases(1).guess.states.row(1)=linspace(deg_to_rad(lon0), deg_to_rad(lon0), 5);
	problem1.phases(1).guess.states.row(2)=linspace(deg_to_rad(lat0), deg_to_rad(lat0), 5);
	problem1.phases(1).guess.states.row(3)=linspace(speed0, speed0, 5);
	problem1.phases(1).guess.states.row(4)=linspace(ee_10, ee_10, 5);
	problem1.phases(1).guess.states.row(5)=linspace(ee_20, ee_20, 5);
	problem1.phases(1).guess.states.row(6)=linspace(ee_30, ee_30, 5);
	problem1.phases(1).guess.states.row(7)=linspace(nn0, nn0, 5);
	problem1.phases(1).guess.states.row(8)=linspace(deg_to_rad(alpha0), deg_to_rad(alpha0), 5);
	problem1.phases(1).guess.states.row(9)=linspace(m0, mdot, 5);
	problem1.phases(1).guess.controls = zeros(2,5);
	problem1.phases(1).guess.controls.row(0) = -0.15*ones(1, 5);
	problem1.phases(1).guess.controls.row(1) = 0*ones(1, 5);
	problem1.phases(1).guess.parameters(0) = -1e-5;
	problem1.phases(1).guess.parameters(1) = 1e-5;
	problem1.phases(1).guess.time = linspace(t0,t_s1, 5);
	
	//Initial Guess Phase 2
	mdot=m_s2-Param.T[1]*(t_s2-t_s1)/(Param.I_sp[1]*Param.g);
	problem1.phases(2).guess.states=zeros(9,5);
	problem1.phases(2).guess.states.row(0)=linspace(alt0, alt0, 5);
	problem1.phases(2).guess.states.row(1)=linspace(deg_to_rad(lon0), deg_to_rad(lon0), 5);
	problem1.phases(2).guess.states.row(2)=linspace(deg_to_rad(lat0), deg_to_rad(lat0), 5);
	problem1.phases(2).guess.states.row(3)=linspace(speed0, speed0, 5);
	problem1.phases(2).guess.states.row(4)=linspace(1, 1, 5);
	problem1.phases(2).guess.states.row(5)=linspace(0, 0, 5);
	problem1.phases(2).guess.states.row(6)=linspace(deg_to_rad(alpha0), deg_to_rad(alpha0), 5);
	problem1.phases(2).guess.states.row(7)=linspace(0, 0, 5);
	problem1.phases(2).guess.states.row(8)=linspace(m_s2, mdot, 5);
	problem1.phases(2).guess.controls = zeros(2,5);
	problem1.phases(2).guess.controls.row(0) = -0.15*ones(1, 5);
	problem1.phases(2).guess.controls.row(1) = 0*ones(1, 5);
	problem1.phases(2).guess.time = linspace(t_s1,t_s2, 5);
	
	//Initial Guess Phase 3
	mdot=m_s3-Param.T[2]*(t_fair-t_s2)/(Param.I_sp[2]*Param.g);
	problem1.phases(3).guess.states=zeros(9,5);
	problem1.phases(3).guess.states.row(0)=linspace(h_atm, h_atm, 5);
	problem1.phases(3).guess.states.row(1)=linspace(deg_to_rad(lon0), deg_to_rad(lon0), 5);
	problem1.phases(3).guess.states.row(2)=linspace(deg_to_rad(lat0), deg_to_rad(lat0), 5);
	problem1.phases(3).guess.states.row(3)=linspace(speed0, speed0, 5);
	problem1.phases(3).guess.states.row(4)=linspace(0.5, 0.5, 5);
	problem1.phases(3).guess.states.row(5)=linspace(0, 0, 5);
	problem1.phases(3).guess.states.row(6)=linspace(deg_to_rad(alpha0), deg_to_rad(alpha0), 5);
	problem1.phases(3).guess.states.row(7)=linspace(0, 0, 5);
	problem1.phases(3).guess.states.row(8)=linspace(m_s3, mdot, 5);
	problem1.phases(3).guess.controls = zeros(2,5);
	problem1.phases(3).guess.controls.row(0) = -0.15*ones(1, 5);
	problem1.phases(3).guess.controls.row(1) = 0*ones(1, 5);
	problem1.phases(3).guess.time = linspace(t_s2,t_fair, 5);
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Register Problem Functions
	problem1.integrand_cost=&integrand_cost_1;
	problem1.endpoint_cost=&endpoint_cost_1;
	problem1.dae=&dae_1;
	problem1.events=&events_1;
	problem1.linkages=&linkages_1;
	
	problem2.integrand_cost=&integrand_cost_2;
	problem2.endpoint_cost=&endpoint_cost_2;
	problem2.dae=&dae_2;
	problem2.events=&events_2;
	problem2.linkages=&linkages_2;
	
	problem3.integrand_cost=&integrand_cost_3;
	problem3.endpoint_cost=&endpoint_cost_3;
	problem3.dae=&dae_3;
	problem3.events=&events_3;
	problem3.linkages=&linkages_3;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Algorithm Options
	algorithm1.nlp_method="IPOPT";
	algorithm1.scaling="automatic";
	algorithm1.derivatives="automatic";
	algorithm1.nlp_iter_max=800;
	algorithm1.diff_matrix = "central-differences";
	//algorithm1.mesh_refinement="automatic";
	//algorithm1.collocation_method="Hermite-Simpson";
	//algorithm1.collocation_method="Chebyshev";
	algorithm1.ipopt_linear_solver="ma27";
	//algorithm1.mr_max_iterations=10;
	
	algorithm2.nlp_method="IPOPT";
	algorithm2.scaling="automatic";
	algorithm2.derivatives="automatic";
	algorithm2.nlp_iter_max=800;
	algorithm2.diff_matrix = "central-differences";
	//algorithm2.mesh_refinement="automatic";
	algorithm2.collocation_method="Hermite-Simpson";
	//algorithm2.collocation_method="Chebyshev";
	algorithm2.ipopt_linear_solver="ma27";
	//algorithm2.mr_max_iterations=10;
	
	algorithm3.nlp_method="IPOPT";
	algorithm3.scaling="automatic";
	algorithm3.derivatives="automatic";
	algorithm3.nlp_iter_max=800;
	algorithm3.diff_matrix = "central-differences";
	//algorithm3.mesh_refinement="automatic";
	//algorithm3.collocation_method="Hermite-Simpson";
	//algorithm3.collocation_method="Chebyshev";
	algorithm3.ipopt_linear_solver="ma27";
	//algorithm3.mr_max_iterations=10;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Calling PSOPT to solve Problem1
	while (1)
	{
		psopt(solution1, problem1, algorithm1);
		if (!solution1.error_flag)
		{
			break;
		}
	}
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	//Extracting relevant variables from Solution1
	MatrixXd x_ph1, x_ph2, x_ph3, x_ph4, x_ph5, x_ph6, x_ph7, x_ph8, x_mass;
	MatrixXd t_ph1, t_ph2, t_ph3, t_ph4, t_ph5, t_ph6, t_ph7, t_ph8;
	MatrixXd u_ph1, u_ph2, u_ph3, u_ph4, u_ph5, u_ph6, u_ph7, u_ph8;
	MatrixXd x, u, t, t_mass;
	MatrixXd conv1, conv2, conv3, conv4;
	int t3f, t6f;
	
	x_ph1 = solution1.get_states_in_phase(1);
	u_ph1 = solution1.get_controls_in_phase(1);
	t_ph1 = solution1.get_time_in_phase(1);
	
	conv1.resize(1, x_ph1.cols());
	conv2.resize(1, x_ph1.cols());
	conv3.resize(1, x_ph1.cols());
	conv4.resize(1, x_ph1.cols());
	for (i=0; i<x_ph1.cols(); i++)
	{
		conv1(0, i)=map_path_ang(x_ph1(4, i), x_ph1(5, i), x_ph1(6, i), x_ph1(7, i));
		conv2(0, i)=map_head_ang(x_ph1(4, i), x_ph1(5, i), x_ph1(6, i), x_ph1(7, i));
		conv3(0, i)=map_bank_ang(x_ph1(4, i), x_ph1(5, i), x_ph1(6, i), x_ph1(7, i));
	}
	
	x_ph2 = solution1.get_states_in_phase(2);
	u_ph2 = solution1.get_controls_in_phase(2);
	t_ph2 = solution1.get_time_in_phase(2);
	
	x_ph3 = solution1.get_states_in_phase(3);
	u_ph3 = solution1.get_controls_in_phase(3);
	t_ph3 = solution1.get_time_in_phase(3);
	t3f	= x_ph3.cols()-1;
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/		
	
	//Events for Problem2
	problem2.phases(1).bounds.lower.events << x_ph3(0, t3f), x_ph3(1, t3f), x_ph3(2, t3f), x_ph3(3, t3f), x_ph3(4, t3f), x_ph3(5, t3f), x_ph3(6, t3f), x_ph3(7, t3f), m_s4;
	problem2.phases(1).bounds.upper.events << x_ph3(0, t3f), x_ph3(1, t3f), x_ph3(2, t3f), x_ph3(3, t3f), x_ph3(4, t3f), x_ph3(5, t3f), x_ph3(6, t3f), x_ph3(7, t3f), m_s4;
	//cout << "\n" << "Initial Conditions for Phase 4: " << problem2.phases(1).bounds.lower.events << "\n";
	
	problem2.phases(2).bounds.lower.events << h_peak_min, 0, 0, 0;
	problem2.phases(2).bounds.upper.events << h_peak_max, 0, 0, 0;
											 		
	problem2.phases(3).bounds.lower.events << h_peak_min, 0, 0, 0, h_atm;
	problem2.phases(3).bounds.upper.events << h_peak_max, 0, 0, 0, h_atm;
	
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/		

	//Initial Guess Phase 4
	mdot=m_s4-Param.T[2]*(t_s3-t_fair)/(Param.I_sp[2]*Param.g);
	problem2.phases(1).guess.states=zeros(10,5);
	problem2.phases(1).guess.states.row(0)=linspace(x_ph3(0, t3f), x_ph3(0, t3f), 5);
	problem2.phases(1).guess.states.row(1)=linspace(x_ph3(1, t3f), x_ph3(1, t3f), 5);
	problem2.phases(1).guess.states.row(2)=linspace(x_ph3(2, t3f), x_ph3(2, t3f), 5);
	problem2.phases(1).guess.states.row(3)=linspace(x_ph3(3, t3f), x_ph3(3, t3f), 5);
	problem2.phases(1).guess.states.row(4)=linspace(x_ph3(4, t3f), x_ph3(4, t3f), 5);
	problem2.phases(1).guess.states.row(5)=linspace(x_ph3(5, t3f), x_ph3(5, t3f), 5);
	problem2.phases(1).guess.states.row(6)=linspace(x_ph3(6, t3f), x_ph3(6, t3f), 5);
	problem2.phases(1).guess.states.row(7)=linspace(x_ph3(7, t3f), x_ph3(7, t3f), 5);
	problem2.phases(1).guess.states.row(8)=linspace(m_s4, mdot, 5);
	problem2.phases(1).guess.controls = zeros(2,5);
	problem2.phases(1).guess.controls.row(0) = -0.15*ones(1, 5);
	problem2.phases(1).guess.controls.row(1) = -0.3*ones(1, 5);
	problem2.phases(1).guess.time = linspace(t_fair,t_s3, 5);
	
	//Initial Guess Phase 5
	problem2.phases(2).guess.states=zeros(10,5);
	problem2.phases(2).guess.states.row(0)=linspace(h_peak_min, h_peak_max, 5);
	problem2.phases(2).guess.states.row(1)=linspace(x_ph3(0, t3f), x_ph3(0, t3f), 5);
	problem2.phases(2).guess.states.row(2)=linspace(x_ph3(1, t3f), x_ph3(1, t3f), 5);
	problem2.phases(2).guess.states.row(3)=linspace(0, 0, 5);
	problem2.phases(2).guess.states.row(4)=linspace(x_ph3(3, t3f), x_ph3(3, t3f), 5);
	problem2.phases(2).guess.states.row(5)=linspace(0, 0, 5);
	problem2.phases(2).guess.states.row(6)=linspace(0, 0, 5);
	problem2.phases(2).guess.states.row(7)=linspace(0, 0, 5);
	problem2.phases(2).guess.controls = zeros(2,5);
	problem2.phases(2).guess.controls.row(0) = 0*ones(1, 5);
	problem2.phases(2).guess.controls.row(1) = -0.3*ones(1, 5);
	problem2.phases(2).guess.time = linspace(t_s3,800, 5);
	
	//Initial Guess Phase 6
	problem2.phases(3).guess.states=zeros(10,5);
	problem2.phases(3).guess.states.row(0)=linspace(h_peak_max, h_atm, 5);
	problem2.phases(3).guess.states.row(1)=linspace(x_ph3(0, t3f), x_ph3(0, t3f), 5);
	problem2.phases(3).guess.states.row(2)=linspace(x_ph3(1, t3f), x_ph3(1, t3f), 5);
	problem2.phases(3).guess.states.row(3)=linspace(0, 0, 5);
	problem2.phases(3).guess.states.row(4)=linspace(x_ph3(3, t3f), x_ph3(3, t3f), 5);
	problem2.phases(3).guess.states.row(5)=linspace(0, 0, 5);
	problem2.phases(3).guess.states.row(6)=linspace(0, 0, 5);
	problem2.phases(3).guess.states.row(7)=linspace(0, 0, 5);
	problem2.phases(3).guess.controls = zeros(2,5);
	problem2.phases(3).guess.controls.row(0) = 0*ones(1, 5);
	problem2.phases(3).guess.controls.row(1) = -0.3*ones(1, 5);
	problem2.phases(3).guess.time = linspace(400,1000, 5);

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/		
		
	//Calling PSOPT to solve Problem2
	while (2)
	{
		psopt(solution2, problem2, algorithm2);
		if (!solution2.error_flag)
		{
			break;
		}
	}
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/		

	//Extracting relevant variables from Solution2
	x_ph4 = solution2.get_states_in_phase(1);
	u_ph4 = solution2.get_controls_in_phase(1);
	t_ph4 = solution2.get_time_in_phase(1);
		
	x_ph5 = solution2.get_states_in_phase(2);
	u_ph5 = solution2.get_controls_in_phase(2);
	t_ph5 = solution2.get_time_in_phase(2);
	
	x_ph6 = solution2.get_states_in_phase(3);
	u_ph6 = solution2.get_controls_in_phase(3);
	t_ph6 = solution2.get_time_in_phase(3);
	t6f	= x_ph6.cols()-1;

/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/		
	
	x.resize(10, x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	x_mass.resize(1, x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols());
	x.row(0) << x_ph1.row(0), x_ph2.row(0), x_ph3.row(0), x_ph4.row(0), x_ph5.row(0), x_ph6.row(0);
	x.row(1) << x_ph1.row(1), x_ph2.row(1), x_ph3.row(1), x_ph4.row(1), x_ph5.row(1), x_ph6.row(1);
	x.row(2) << x_ph1.row(2), x_ph2.row(2), x_ph3.row(2), x_ph4.row(2), x_ph5.row(2), x_ph6.row(2);
	x.row(3) << x_ph1.row(3), x_ph2.row(3), x_ph3.row(3), x_ph4.row(3), x_ph5.row(3), x_ph6.row(3);	
	x.row(4) << conv1, x_ph2.row(4), x_ph3.row(4), x_ph4.row(4), x_ph5.row(4), x_ph6.row(4);
	x.row(5) << conv2, x_ph2.row(5), x_ph3.row(5), x_ph4.row(5), x_ph5.row(5), x_ph6.row(5);
	x.row(6) << x_ph1.row(8), x_ph2.row(6), x_ph3.row(6), x_ph4.row(6), x_ph5.row(6), x_ph6.row(6);
	x.row(7) << conv3, x_ph2.row(7), x_ph3.row(7), x_ph4.row(7), x_ph5.row(7), x_ph6.row(7);
	x_mass << x_ph1.row(9), x_ph2.row(8), x_ph3.row(8), x_ph4.row(8);
	
	t.resize(1, t_ph1.cols()+t_ph2.cols()+t_ph3.cols()+t_ph4.cols()+t_ph5.cols()+t_ph6.cols());
	t << t_ph1, t_ph2, t_ph3, t_ph4, t_ph5, t_ph6;
	t_mass.resize(1, t_ph1.cols()+t_ph2.cols()+t_ph3.cols()+t_ph4.cols());
	t_mass << t_ph1, t_ph2, t_ph3, t_ph4;
	
	u.resize(2, u_ph1.cols()+u_ph2.cols()+u_ph3.cols());
	u.row(0)<<u_ph1.row(0), u_ph2.row(0), u_ph3.row(0);
	u.row(1)<<u_ph1.row(1), u_ph2.row(1), u_ph3.row(1);
	
	Save(x,"x.dat");
	Save(t,"t.dat");
	Save(u,"u.dat");
	
/*************************************************************************************************************************************************/
/*************************************************************************************************************************************************/	
	
	MatrixXd altitude, longitude, latitude, speed_state, path_angle, heading_angle, attack_angle, bank_angle;
	MatrixXd attack_angle_rate, bank_angle_rate;
    
	altitude = x.block(0,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	longitude = x.block(1,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols()); 
	latitude = x.block(2,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	speed_state = x.block(3,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
    path_angle = x.block(4,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	heading_angle = x.block(5,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	attack_angle = x.block(6,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	bank_angle = x.block(7,0,1,x_ph1.cols()+x_ph2.cols()+x_ph3.cols()+x_ph4.cols()+x_ph5.cols()+x_ph6.cols());
	
    plot(t,altitude,problem1.name, "time (s)", "Altitude (m)");
	plot(t,longitude,problem1.name, "time (s)", "Longitude (rad)");
	plot(t,latitude,problem1.name, "time (s)", "Latitude (rad)");
	plot(t,speed_state,problem1.name, "time (s)", "Speed (m/s)");
	plot(t,path_angle,problem1.name, "time (s)", "Path Angle (rad)");
	plot(t,heading_angle,problem1.name, "time (s)", "Heading Angle (rad)");
	plot(t,attack_angle,problem1.name, "time (s)", "Angle of Attack (rad)");
	plot(t,bank_angle,problem1.name, "time (s)", "Bank Angle (rad)");
	plot(t_mass,x_mass,problem1.name, "time (s)", "Mass (kg)");
	
	return 0;
}
