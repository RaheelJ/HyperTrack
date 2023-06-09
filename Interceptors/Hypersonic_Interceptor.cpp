#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "Hypersonic_Interceptor.h"
#include <iostream>
#include "xml/pugixml.hpp"

void Hypersonic::Interceptor::Reset() {
    /* Variables to calculate metrics */
    dist_LOS = -1;
    speed = 0;
    altitude = 0;
    dist_covered = 0;
    time_taken = 0;
    accuracy = 0;
    num = 1;

    head_ang = 0;
    path_ang = 0;
    x_dist = 0;
    y_dist = 0;
    energy = 0;

    /* Command related data */
    bank_ang = 0;
    L = 0;
    lift_coeff = 0;
    attack_ang = 0;
    n_long = 0;
    n_lat = 0;
    Vm = 0;
    D = 0;

    /* Measurement related data */
    delta_V = 0;
    s = 0;
    delta_h = 0;
    q_dot_long = 0;
    q_long = 0;
    q_dot_lat = 0;
    q_lat = 0;
    phase = 0;

	q_lat_last = 0;
	q_long_last = 0;
	t_last = 0;

    is_first_run = true;
    Last_Target_State.id = -1;
    my_terminated = false;
	my_killed = false;
    self_destruct = false;
    valid_csv = false;
}
bool Hypersonic::Interceptor::Initialize(const std::string& xi_init_file_name, std::string& xo_message) {
    pugi::xml_document Doc;
    int error = 0;
    try
    {
        pugi::xml_parse_result Result = Doc.load_file(xi_init_file_name.c_str());	/* Read "Parameters.xml" */
        if (Result) {}
        else
        {
			throw "File not found or formatting error!";	                        /* Throw exception if file not read */
        }
    }
    catch (const char* error_2)
    {
        cerr << error_2 << "\n";
		xo_message = "File not found or formatting error!";	                    /* Throw exception if file not read */
		return false;                                                           /* Print the error message */
    }
    pugi::xml_node Root = Doc.child("Parameters");                              /* Accessing Root Node */
    pugi::xml_node Common = Root.child("Common");
    pugi::xml_node Algo = Root.child("Hypersonic");
    if (Root && Common && Algo) {}
    else {
        xo_message = "Root Node (Parameters) not found!";
        return false;
    }

	pugi::xml_node Position = Common.child("Position");
    pugi::xml_node Param_Boost = Common.child("Boost_Data");
    pugi::xml_node Exec_Param = Root.child("Selection_Execution");
    pugi::xml_node Aiming_Data = Algo.child("Aiming_Data");
    pugi::xml_node Param_Lift = Algo.child("Param_Lift");
    pugi::xml_node Corridor = Algo.child("Corridor");
    pugi::xml_node Param_Jet = Algo.child("Param_Jet");
    if (Position && Param_Boost && Aiming_Data && Param_Lift && Corridor && Param_Jet && Exec_Param) {}
    else {
        xo_message = "One of the child nodes is missing!";
        return false;
    }

    /* Initialize parameters from "Position" node */
    Initial_State.y_or_lon = Read_Param(Position.attribute("Longitude").as_double(), -180, 180, error);
    Initial_State.x_or_lat = Read_Param(Position.attribute("Latitude").as_double(), -90, 90, error);
    Initial_State.z_or_alt = Read_Param(Position.attribute("Altitude").as_double(), 0, 6378100, error);

    /* Initialize parameters from "Boost_Data" node */
    boost_vel = Read_Param(stod(Param_Boost.child_value("Boost_Speed")), 0.001, 100000, error);
    boost_alt = Read_Param(stod(Param_Boost.child_value("Boost_Altitude")), 0, 6378100, error);

    /* Initialize parameters from "Aiming_Data" node */
    K1 = Read_Param(stod(Aiming_Data.child_value("K1")), 0, 100, error);
    K2 = Read_Param(stod(Aiming_Data.child_value("K2")), 0, 100, error);

    /* Initialize parameters from "Param_Lift" node */
    k1 = Read_Param(stod(Param_Lift.child_value("k1")), 0, 1000, error);
    k2 = Read_Param(stod(Param_Lift.child_value("k2")), 0, 100, error);
    m = Read_Param(stod(Param_Lift.child_value("Mass")), 0, stod(Param_Lift.child_value("Mass")), error);
    Sa = Read_Param(stod(Param_Lift.child_value("Reference_Area")), 0.001, stod(Param_Lift.child_value("Mass")), error);
    Re = 6371000;
    atm_constant = 1.225;
    Hs = 8500;
    g = 9.8;

    /* Initialize parameters from "Corridor" node */
    s0 = Read_Param(stod(Exec_Param.child_value("Start_Distance")), 1000000, 10000000, error);
    s1 = Read_Param(stod(Corridor.child_value("Glide_Distance")), 800000, s0, error);
    s2 = Read_Param(stod(Corridor.child_value("Aim_Distance")), 200000, s1, error);
    s3 = Read_Param(stod(Corridor.child_value("Kill_Distance")), 20000, s2, error);
    bias = Read_Param(stod(Corridor.child_value("Bias")), 0, 1000, error);
    slope = Read_Param(stod(Corridor.child_value("Slope")), 0, 2000, error);
    dead_band = Read_Param(stod(Corridor.child_value("Dead_Band")), s2, s1, error);
    gain = Read_Param(stod(Corridor.child_value("Gain")), 0, 60, error);

    /* Initialize parameters from "Param_Jet" node */
    N_long = Read_Param(stod(Param_Jet.child_value("Load_Long")), 0, 6, error);
    N_lat = Read_Param(stod(Param_Jet.child_value("Load_Lat")), 0, 6, error);
    delta_Q = Read_Param(stod(Param_Jet.child_value("ON_Threshold")), 0, 0.6, error);
    delta_q = Read_Param(stod(Param_Jet.child_value("OFF_Threshold")), 0, 0.2, error);

    /* Initialize parameters from "Execution" node */
    hit_radius = Read_Param(stod(Exec_Param.child_value("Stop_Distance")), 0, 5000, error);
    energy_lim = Read_Param(stod(Exec_Param.child_value("Energy_Limit")), 0, 10000000, error);
    //save_csv = Exec_Param.attribute("Save_CSV").as_bool();

    /* Check for out of range parameters */
    if (error == 1) {
        xo_message = to_string(error) + " parameter is out of range";
        return false;
    }
    else if (error > 1) {
        xo_message = to_string(error) + " parameters are out of range";
        return false;
    }
    
    /* Read the relatioship between Lift Coefficient and Angle of Attack */
    //vector<double> LC(18, 0);                                                               /* Array for Lift Coefficients Data of CAV-H Model */
    //vector<double> MN = { 3.5,3.5,3.5, 5,5,5, 8,8,8, 10,10,10, 15,15,15, 20,20,20 };        /* Array for Speed Data of CAV-H Model */
    //vector<double> AoA = { 10,15,20, 10,15,20, 10,15,20, 10,15,20, 10,15,20, 10,15,20 };    /* Array for Angle of Attack Data of CAV-H Model */
    //pugi::xml_attribute_iterator ait;
    //pugi::xml_node_iterator it;
    //it = pugi::xml_node_iterator(Lift_Coefficients);
    //int i = 0;
    //for (ait = it->attributes_begin(); ait != it->attributes_end(); ++ait)
    //{
    //    LC[i] = ait->as_double();
    //    i = i + 1;
    //}
    //Interpolate.setData(LC, MN, AoA);

	Initial_State.id = -1;
	Initial_State.time = 0;
	Initial_State.state_in_geo = true;
    State = Initial_State;

    Reset();
	return true;
};
void Hypersonic::Interceptor::Reinitialize() {
    /* Reinitialize parameters from "Position" node */
    State = Initial_State;
    Reset();
}

void Hypersonic::Interceptor::Boost() {                     
    State.z_or_alt = boost_alt;
    Vm = boost_vel;
    head_ang = rad_to_deg(q_lat);
    path_ang = rad_to_deg(q_long);
};

double Hypersonic::Interceptor::V_Lim() {				 /* Returns  Velocity Bound for Banking Commands based on s */
    double bound;
    bound = bias + s / s1 * slope;
    if (s <= dead_band) { bound = bias; }
    return bound;
};
double Hypersonic::Interceptor::Bank_Ref() {			     /* Returns reference for Banking Command based on s */
    double ref;
    ref = 80 - gain * s / s1;
    return ref;
};
void Hypersonic::Interceptor::Calc_Bank() {			     /* Calulates Banking Command for Gliding Phase*/
    if (delta_V >= V_Lim()) {
        bank_ang = Bank_Ref();
    }
    else if (delta_V > -V_Lim() && delta_V < V_Lim()) {   /* Bank_Ref(s) * signum(bank_ang) */
        bank_ang = Bank_Ref() * ((bank_ang > 0) - (bank_ang < 0));
    }
    else {
        bank_ang = -Bank_Ref();
    }
}
double Hypersonic::Interceptor::Atm_Den() {			     /* Returns Atmospheric Density */
    return atm_constant * exp(-State.z_or_alt / Hs);
};
void Hypersonic::Interceptor::Calc_Lift() {			     /* Calculate Lift using Modified QEGC and then Calulate Lift Coefficient */
    double r = State.z_or_alt + Re;
    L = (1 / cosd(bank_ang)) * ((-k1 * delta_h / s) - (k2 * deg_to_rad(path_ang)) - (pow(Vm, 2) * cosd(path_ang) / r) + (g * cosd(path_ang)));
    lift_coeff = 2 * m * L / (Sa * Atm_Den() * pow(Vm, 2));
}
void Hypersonic::Interceptor::Calc_AoA() {			     /* Calculate Angle of Attack from Lift Coefficient */
    double Mach = Vm / 343;
    //attack_ang = Interpolate(lift_coeff, Mach);
}
void Hypersonic::Interceptor::Glide()
{
    Calc_Bank();
    Calc_Lift();
}

void Hypersonic::Interceptor::Aim() {				     /* Calculate Banking and Lift commands using PNG in the Aiming mode */
    double Lcos_bank, Lsin_bank, tan_bank;
    double r = State.z_or_alt + Re;

    Lcos_bank = Vm * K1 * q_dot_long + g - pow(Vm, 2) / r;
    Lsin_bank = Vm * K2 * q_dot_lat - (pow(Vm, 2) * sind(head_ang) * tand(State.x_or_lat) / r);
    tan_bank = Lsin_bank / Lcos_bank;
    bank_ang = rad_to_deg(atan(tan_bank));
    L = Lcos_bank / cosd(bank_ang);
    lift_coeff = 2 * m * L / (Sa * Atm_Den() * pow(Vm, 2));
};

void Hypersonic::Interceptor::Jet_Law(double q_dot, double N, double& n) {	/* Switching Control Law of the Reaction-Jet */
    //n = N * ((q_dot > 0) - (q_dot < 0));

    double N_signum;
    N_signum = N * ((q_dot > 0) - (q_dot < 0));
    if (n == 0 && abs(q_dot) <= delta_Q) { n = 0; }
    else if (n == 0 && abs(q_dot) > delta_Q) { n = N_signum; }
    else if (n != 0 && abs(q_dot) < delta_q) { n = 0; }
    else if (n != 0 && abs(q_dot) >= delta_q) { n = N_signum; }
}
void Hypersonic::Interceptor::Kill() {				    /* Change Velocity direction using Reaction Jet in Hit to Kill Phase */
    Jet_Law(q_dot_long, N_long, n_long);
    Jet_Law(q_dot_lat, N_lat, n_lat);
};

void Hypersonic::Interceptor::Course::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] -- Velocity heading angle (head_ang)
        x[1] -- Flight-path angle (path_ang)
        x[2] -- Radial distance from centre of the Earth (Pos.altitude + Re)
        x[3] -- Longitude (Pos.longitude)
        x[4] -- Latitude (Pos.latitude) 
        x[5] -- Speed 
        x[6] -- x distance
        x[7] -- y distance */

    dx[0] = (1 / x[5]) * (n_lat * g + L * sind(bank_ang) / cos(x[1]) + pow(x[5], 2) * cos(x[1]) * sin(x[0]) * tan(x[4]) / x[2]);
    dx[1] = (1 / x[5]) * (n_long * g + L * cosd(bank_ang) - g * cos(x[1]) + pow(x[5], 2) * cos(x[1]) / x[2]);
    dx[2] = x[5] * sin(x[1]);
    dx[3] = x[5] * cos(x[1]) * sin(x[0]) / (x[2] * cos(x[4]));
    dx[4] = x[5] * cos(x[1]) * cos(x[0]) / x[2];
    dx[5] = 0;
    dx[6] = x[5] * cos(x[1]) * sin(x[0]);
    dx[7] = x[5] * cos(x[1]) * cos(x[0]);
}
void Hypersonic::Interceptor::Check_Param() {
    bool condition_1 = (Vm < 0.000001 && phase > 1) || (State.z_or_alt < 0) || (abs(State.x_or_lat) >= 89.999999);
    bool condition_2 = (89.999999 <= abs(path_ang) && abs(path_ang) <= 90.000001);
    bool condition_3 = (269.999999 <= abs(path_ang) && abs(path_ang) <= 270.000001);
    self_destruct = condition_1 || condition_2 || condition_3;
}
void Hypersonic::Interceptor::Track(double t_meas) {
    state_t x = { 0,0,0,0,0,0,0,0 };
    x[0] = deg_to_rad(head_ang);    /* Initial state values for the current iteration */
    x[1] = deg_to_rad(path_ang);
    x[2] = State.z_or_alt + Re;
    x[3] = deg_to_rad(State.y_or_lon);
    x[4] = deg_to_rad(State.x_or_lat);
    x[5] = Vm;
    x[6] = x_dist;
    x[7] = y_dist;

    double t = State.time;
    double step_size = t_meas - t;
    RK4 Integrator;

    int last_phase = phase;

    phase = (int)(s <= s0) + (int)(s <= s1) + (int)(s <= s2) + (int)(s <= s3);
    if ((last_phase == 0) && (phase > 1)) { phase = 1; }
    
    Check_Param();
    if (self_destruct) {
        return;
    }
       
    while (t < t_meas - 0.000001)                /* While flight_time <= t_meas */
    {
        if (phase == 0) {                        /* Do nothing if Target is out of range */
            t = t_meas;
            break;
        }
        else if (phase == 1) {                   /* Enter the Boost phase */
            if (save_csv) {
                Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[2], rad_to_deg(x[3]), rad_to_deg(x[4]), L, bank_ang,
                       s, delta_h, delta_V, rad_to_deg(q_lat), rad_to_deg(q_long), n_lat, n_long });
                valid_csv = true;
            }
            t = t_meas;
            Boost();
            x[0] = deg_to_rad(head_ang);
            x[1] = deg_to_rad(path_ang);
            x[2] = State.z_or_alt + Re;
            x[5] = Vm;
            //cout << "Boost Mode \n";
            break;
        }
        else if (abs(State.time - t) < 0.000001) {    /* First calculate the commands based on the mode */
            switch (phase)
            {
            case 2:
                Glide();
                //cout << "Glide Mode \n";
                break;
            case 3:
                Aim();
                //cout << "Aim Mode \n";
                break;
            case 4:
                Kill();
                //cout << "Kill Mode \n";
                break;
            case 5:
                //cout << "Test Flight with Initial Conditions";
                break;
            default:
                break;
            }
            Calc_AoA();                             /* Find Angle of Attack from Lift Coefficient using the table */
        }
        if (save_csv) {
            Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[2], rad_to_deg(x[3]), rad_to_deg(x[4]), L, bank_ang,
                   s, delta_h, delta_V, rad_to_deg(q_lat), rad_to_deg(q_long), n_lat, n_long });
            valid_csv = true;
        }
        Integrator(Course(), x, t, step_size);      /* Solve the ODEs to find system states */
        energy = sqrt(pow(x[6], 2) + pow(x[7], 2) + pow(boost_alt - (x[2] - Re), 2));
        //Limit(x[0], M_PI_2, -M_PI_2); Limit(x[1], M_PI_2, -M_PI_2); Limit(x[3], M_PI, -M_PI); Limit(x[4], M_PI_2, -M_PI_2);
    }

    head_ang = rad_to_deg(x[0]);                    /* Save the system states in appropriate variables for next iteration */
    path_ang = rad_to_deg(x[1]);
    Vm = x[5];
    x_dist = x[6];
    y_dist = x[7];

    State.id = 0;
    State.time = t;
    State.z_or_alt = x[2] - Re;
    State.y_or_lon = rad_to_deg(x[3]);
    State.x_or_lat = rad_to_deg(x[4]);

    if (phase > 1) {
        State.x_vel = Vm * cosd(path_ang) * sind(head_ang);
        State.y_vel = Vm * cosd(path_ang) * cosd(head_ang);
        State.z_vel = Vm * sind(path_ang);
    }
    else {
        State.x_vel = 0;
        State.y_vel = 0;
        State.z_vel = 0;
    }
}
void Hypersonic::Interceptor::Fine_Calculations(objectstate Last_Interceptor_State, objectstate Target_State) {
    double fine_time_steps = 2000, i = 0;
    double t1 = Last_Interceptor_State.time, t2 = Target_State.time;

    objectstate Fine_Target = Last_Target_State;
    double fine_target_lat = (Target_State.x_or_lat - Fine_Target.x_or_lat) / fine_time_steps;
    double fine_target_long = (Target_State.y_or_lon - Fine_Target.y_or_lon) / fine_time_steps;
    double fine_target_alt = (Target_State.z_or_alt - Fine_Target.z_or_alt) / fine_time_steps;

    objectstate Fine_Interceptor = Last_Interceptor_State;
    double fine_interceptor_lat = (State.x_or_lat - Fine_Interceptor.x_or_lat) / fine_time_steps;
    double fine_interceptor_long = (State.y_or_lon - Fine_Interceptor.y_or_lon) / fine_time_steps;
    double fine_interceptor_alt = (State.z_or_alt - Fine_Interceptor.z_or_alt) / fine_time_steps;

    for (i = t1 + (t2 - t1) / fine_time_steps; i < (t2 - 0.000001); i = i + (t2 - t1) / fine_time_steps) {
        Fine_Target.time = i;
        Fine_Target.x_or_lat = Fine_Target.x_or_lat + fine_target_lat;
        Fine_Target.y_or_lon = Fine_Target.y_or_lon + fine_target_long;
        Fine_Target.z_or_alt = Fine_Target.z_or_alt + fine_target_alt;

        Fine_Interceptor.time = i;
        Fine_Interceptor.x_or_lat = Fine_Interceptor.x_or_lat + fine_interceptor_lat;
        Fine_Interceptor.y_or_lon = Fine_Interceptor.y_or_lon + fine_interceptor_long;
        Fine_Interceptor.z_or_alt = Fine_Interceptor.z_or_alt + fine_interceptor_alt;

        Calc_Distance(Fine_Interceptor, Fine_Target, s, q_lat);
        q_lat = deg_to_rad(q_lat + 360 * (q_lat < 0));
        delta_h = Fine_Interceptor.z_or_alt - Fine_Target.z_or_alt;
        if (save_csv) {
            Rec1({ i, head_ang, path_ang, Fine_Interceptor.z_or_alt + Re, Fine_Interceptor.y_or_lon, Fine_Interceptor.x_or_lat, L, bank_ang,
                s, delta_h, delta_V, rad_to_deg(q_lat), rad_to_deg(q_long), n_lat, n_long });
            valid_csv = true;
        }

        my_terminated = Is_Hit() || self_destruct;
        my_killed = Is_Hit();
        if (my_terminated)
        {
            /* Extract variables required to calculate metrics */
            if (s < dist_LOS || dist_LOS < 0)
            {
                dist_LOS = s;
                dist_covered = energy;
                time_taken = State.time;
            }
            speed = Vm;
            altitude = State.z_or_alt;
            accuracy = 1;
            num = 1;

            break;
        }
    }
}
void Hypersonic::Interceptor::First_Run(objectstate Target_State) {
    State.time = Target_State.time;
    Calc_Distance(State, Target_State, s, q_lat);
    if (q_lat < 0) { q_lat = q_lat + 360; }
    q_lat = deg_to_rad(q_lat);
    delta_h = State.z_or_alt - Target_State.z_or_alt;
    q_long = atan(-delta_h / s);

    t_last = State.time;
    q_lat_last = q_lat;
    q_long_last = q_long;
    Last_Target_State = Target_State;

    is_first_run = false;
}
void Hypersonic::Interceptor::Update_Target_State(objectstate Target_State)
{
    int last_phase;
    double t_meas = Target_State.time;
    double delta_t = t_meas - t_last;
    objectstate Last_Interceptor_State;
    double delta_q_lat;

	if (delta_t <= 0)
	{
		return;
	}

    if (my_terminated || (t_meas < 0))
    {
        self_destruct = 1;
        my_terminated = 1;
        return;
    }

    if (is_first_run)
    {
        First_Run(Target_State);
        return;
    }

    Last_Interceptor_State = State;
    last_phase = phase;
    Track(t_meas);

    if (s < 20000) {
        Fine_Calculations(Last_Interceptor_State, Target_State);
        if (my_terminated)
        {
            return;
        }
    }

    Calc_Distance(State, Target_State, s, q_lat);
    delta_h = State.z_or_alt - Target_State.z_or_alt;
    my_terminated = Is_Hit() || Is_Missed(last_phase) || self_destruct || (energy > energy_lim);
    my_killed = Is_Hit();
    //cout << "Radius: " << sqrt(pow(s, 2) + pow(delta_h, 2)) << "\n";
    //cout << "Energy: " << energy << "\n";

    q_lat = deg_to_rad(q_lat + 360 * (q_lat < 0));
    q_long = atan(-delta_h / s);

    if (delta_t > 0.000001) {
        delta_q_lat = (q_lat - q_lat_last) + ((q_lat - q_lat_last) > 6) * (-2 * M_PI) + ((q_lat - q_lat_last) < -6) * (2 * M_PI);
        q_dot_long = (q_long - q_long_last) / delta_t;
        q_dot_lat = (delta_q_lat) / delta_t;
        delta_V = q_dot_lat * s;
    }
    
    t_last = t_meas;
    q_lat_last = q_lat;
    q_long_last = q_long;
    Last_Target_State = Target_State;

    /* Extract variables required to calculate metrics */
    if (s < dist_LOS || dist_LOS < 0)
    {
        dist_LOS = s;
        dist_covered = energy;
        time_taken = State.time;
    }
    speed = Vm;
    altitude = State.z_or_alt;
    accuracy = 1;
    num = 1;
};

bool Hypersonic::Interceptor::Get_State(objectstate& xo_state) {
	bool valid_state = Is_Started() && !my_killed;
    objectstate State_Out;
    State_Out = State;
    State_Out.y_or_lon = deg_to_long(State.y_or_lon);
    State_Out.x_or_lat = deg_to_lat(State.x_or_lat);
	if (valid_state)
	{
		xo_state = State_Out;
	}
	return valid_state;
}
void Hypersonic::Interceptor::Get_Record() {
    if (save_csv && valid_csv)
    {   
        Rec1.csv("Output/Hypersonic_States", { "Time", "Heading_Angle", "Flight_Path_Angle", "Radial_Distance", "Longitude", "Latitude", "Lift", "Bank_Angle", "LOS_Distance_Lateral",
                    "Altitude_Error", "Delta_V_Lateral", "LOS_Angle_Lateral", "LOS_Angle_Longitudnal", "Reaction_Jet_Lateral", "Reaction_Jet_Longitudnal" });
    }
    
}

bool Hypersonic::Interceptor::Is_Started() {
    return (phase > 0);
}
bool Hypersonic::Interceptor::Is_Hit() {
    double radius = sqrt(pow(s, 2) + pow(delta_h, 2));
    return (radius < hit_radius);
}
bool Hypersonic::Interceptor::Is_Missed(int last_phase) {
    //return (phase < last_phase);
    if (phase == 3 && last_phase == 4)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

bool Hypersonic::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
{
	xo_interceptor = my_terminated;
	xo_target = my_killed;
	if (xo_target)
	{
		xo_target_id = Last_Target_State.id;
		
	}
	else
	{
		xo_target_id = -1;
	}

	return true;
}
