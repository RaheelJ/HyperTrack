#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "Coordinated_Interceptor.h"
#include <iostream>
#include "xml/pugixml.hpp"

void Coordinated::Interceptor::Reset() {
    /* Variables to calculate metrics */
    dist_LOS = -1;
    speed = 0;
    altitude = 0;
    dist_covered = 0;
    time_taken = 0;
    accuracy = 1;
    num = 1;

    command_L = 0;
    head_ang_L = 0;
    R_L = 0;
    q_L = 0;
    phase = 0;
    P_L[0] = 0;
    P_L[1] = 0;
    energy = 0;

    u_F = MatrixXd::Constant(2, n, 0);
    R_F.assign(n, 0);
    P_F_Last = MatrixXd::Constant(2, n, 0);

    q_dot = 0;
    lead_ang_L = 0;
    lead_ang_T = 0;
    V_L = 0;
    a_T = 0;
    head_ang_T = 0;
    head_ang_F = 0;
    V_F = 0;
    u_x = 0;
    u_y = 0;
    t_last = 0;
    q_last = 0;
        
    Last_Target_State.id = -1;
    my_terminated = false;
    my_killed = false;
    self_destruct = false;
    valid_csv = false;
    is_first_run = true;
}
void Coordinated::Interceptor::Initialize_Followers() {
    State_F.resize(n);
    R_F.resize(n);
    P_F.resize(2, n);
    P_F_Last.resize(2, n);
    Initial_P_F.resize(2, n);
    u_F.resize(2, n);
    W.resize(n + 1, n + 1);
    W = MatrixXd::Constant(n + 1, n + 1, 0);
    W(0, n) = 1;
    int i;
    for (i = 1; i < n; i++) {
       W(i, i - 1) = 1;
    }
    alpha.resize(n);
    alpha = VectorXd::Constant(n, 0);
    alpha(0) = 1;
    beta.resize(n);
    gamma.resize(n);
}

bool Coordinated::Interceptor::Initialize(const string& XML_Parameters, string& message) {
    pugi::xml_document Doc;
    int error = 0;
    try
    {
        pugi::xml_parse_result Result = Doc.load_file(XML_Parameters.c_str());	/* Read "Parameters.xml" */
        if (Result) {}
        else {
            throw "File not found or formatting error!";	                    /* Throw exception if file not read */
        }
    }
    catch (const char* error_2)
    {
        message = "File not found or formatting error!";
        cerr << error_2 << "\n";                        /* Print the error message */
        return false;
    }
    pugi::xml_node Root = Doc.child("Parameters");      /* Accessing Root Node */
    pugi::xml_node Common = Root.child("Common");
    pugi::xml_node Algo = Root.child("Coordinated");
    if (Root && Common && Algo) {}
    else {
        message = "Root Node (Parameters) not found!";
        return false;
    }

    pugi::xml_node Position = Common.child("Position");
    pugi::xml_node Boost_Data = Common.child("Boost_Data");
    pugi::xml_node Exec_Param = Root.child("Selection_Execution");
    pugi::xml_node Leader = Algo.child("Leader");
    if (Position && Boost_Data && Exec_Param && Leader) {}
    else {
        message = "One of the child nodes is missing!";
        return false;
    }

    /* Initialize parameters from "Position" and "Boost_Data" node */
    Initial_State_L.x_or_lat = Read_Param(Position.attribute("Latitude").as_double(), -90, 90, error);
    Initial_State_L.y_or_lon = Read_Param(Position.attribute("Longitude").as_double(), -180, 180, error);
    Initial_State_L.z_or_alt = Read_Param(Position.attribute("Altitude").as_double(), 0, 6378100, error);
    Initial_Velocity_L = Read_Param(stod(Boost_Data.child_value("Boost_Speed")), 0, 10000, error);

    /* Initialize parameters from "Leader" node */
    N = Read_Param(stod(Leader.child_value("Proportional_Gain")), 0, 100, error);
    K = Read_Param(stod(Leader.child_value("Scale_Coefficient")), 0, 100, error);
    //n = (int)Read_Param(Leader.attribute("Total_Followers").as_double(), 1, 10, error);
    Re = 6378100;

    /* Initialize parameters from "Execution" node */
    R_start = Read_Param(stod(Exec_Param.child_value("Start_Distance")), 0, 10000000, error);
    R_stop = Read_Param(stod(Exec_Param.child_value("Stop_Distance")), 0, 5000, error);
    energy_lim = Read_Param(stod(Exec_Param.child_value("Energy_Limit")), 0, 10000000, error);
    //save_csv = Exec_Param.attribute("Save_CSV").as_bool();

	n=10;
    Initialize_Followers();
	
	n=0;
	int i=0;
	for(pugi::xml_node Follower = Algo.first_child(); Follower; Follower = Follower.next_sibling())
	{
		if (n>0)
		{				
			Initial_P_F(0, i) = stod(Follower.child_value("x"));
			Initial_P_F(1, i) = stod(Follower.child_value("y"));
            //cout << Initial_P_F(0, i) << "\t" << Initial_P_F(1, i) << "\n";
			error = error + (int)(abs(Initial_P_F(0, i)) < 1000 || abs(Initial_P_F(0, i)) > 1000000) + (int)(abs(Initial_P_F(1, i)) < 1000 || abs(Initial_P_F(1, i)) > 1000000);
			beta(i) = Read_Param(stod(Follower.child_value("Beta_Gain")), 0, 10, error);
			gamma(i) = Read_Param(stod(Follower.child_value("Gamma_Gain")), 0, 10, error);
			i=i+1;
		}
		n=n+1;
	}
    n = n-1;

	// int i;
    // const string temp1 = "Follower_";
    // string temp2;
    // int temp3;
    // for (i = 0; i < n; i = i + 1) {
        // temp2 = temp1 + to_string(i+1);
        // Follower[i] = Algo.child(temp2.c_str());
        // Initial_P_F(0, i) = Follower[i].attribute("x").as_double();
        // Initial_P_F(1, i) = Follower[i].attribute("y").as_double();
        // error = error + (int)(abs(Initial_P_F(0, i)) < 1000 || abs(Initial_P_F(0, i)) > 1000000) + (int)(abs(Initial_P_F(1, i)) < 1000 || abs(Initial_P_F(1, i)) > 1000000);
        // beta(i) = Read_Param(Follower[i].attribute("Beta_Gain").as_double(), 0, 10, error);
        // gamma(i) = Read_Param(Follower[i].attribute("Gamma_Gain").as_double(), 0, 10, error);
    // }
    
    /* Check for out of range parameters */
    if (error == 1) {
        message = to_string(error) + " parameter is out of range";
        return false;
    }
    else if (error > 1) {
        message = to_string(error) + " parameters are out of range";
        return false;
    }

    Initial_State_L.id = -1;
    Initial_State_L.time = 0;
    Initial_State_L.state_in_geo = true;
    State_L = Initial_State_L;
    radial_distance = State_L.z_or_alt + Re;

    for (i = 0; i < n; i++) {
        P_F(all, i) = Initial_P_F(all, i);
        Calc_Geo(State_F[i], {P_F(0, i), P_F(1, i)}, Initial_State_L);
        State_F[i].z_or_alt = Initial_State_L.z_or_alt;
        State_F[i].id = -1;
        State_F[i].time = 0;
        State_F[i].state_in_geo = true;
    }

    Reset();
    return true;
};
void Coordinated::Interceptor::Reinitialize() {
    State_L = Initial_State_L;
    int i;
    for (i = 0; i < n; i++) {
        P_F(all, i) = Initial_P_F(all, i);
        Calc_Geo(State_F[i], { P_F(0, i), P_F(1, i) }, Initial_State_L);
        State_F[i].z_or_alt = Initial_State_L.z_or_alt;
        State_F[i].id = -1;
        State_F[i].time = 0;
        State_F[i].state_in_geo = true;
    }
    Reset();
}

void Coordinated::Interceptor::MPNG() {
    double cos_L = cos(lead_ang_L);
    if (abs(cos_L) < 0.000001) {
        self_destruct = true;
        return;
    }
    double cos_T = cos(lead_ang_T);
    command_L = V_L * q_dot * (1 + N / cos_L) + a_T * cos_T / cos_L;
    Limit(command_L, 10 * 9.8, -10 * 9.8);
}
void Coordinated::Interceptor::Course_L::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] -- Latitude (radians)
        x[1] -- Longitude (radians)
        x[2] -- Velocity Heading Angle of Interceptor (radians)
        x[3] -- Speed 
        x[4] -- x distance
        x[5] -- y distance */

    dx[0] = x[3] * cos(x[2]) / radial_distance;
    dx[1] = x[3] * sin(x[2]) / (radial_distance * cos(x[0]));
    dx[2] = command_L / x[3];
    dx[3] = 0;
    dx[4] = x[3] * sin(x[2]);
    dx[5] = x[3] * cos(x[2]);
}
void Coordinated::Interceptor::Course_F::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] -- Latitude (radians)
        x[1] -- Longitude (radians)
        x[2] -- x Coordinate
        x[3] -- y Coordinate  */

    dx[0] = V_F * cos(head_ang_F) / radial_distance;
    dx[1] = V_F * sin(head_ang_F) / (radial_distance * cos(x[0]));
    dx[2] = u_x;
    dx[3] = u_y;
}
void Coordinated::Interceptor::CGS(int i) {
    int j = 0;
    double temp1, temp2;
    Vector2d P_j_dot, P_L_dot, temp3, temp4, temp5;
    j = 0;
    temp1 = alpha(i) * beta(i);
    temp2 = temp1 + W(i, all).sum();
    if (temp2 < 0.001) {
        self_destruct = true;
        return;
    }
    P_L_dot << State_L.x_vel, State_L.y_vel;
    while (!W(i, j)) {
        j = j + 1;
        if (j >= n) {
            j = 0;
            break;
        }
    }
    P_j_dot << State_F[j].x_vel, State_F[j].y_vel;
    temp5(0) = P_L[0];
    temp5(1) = P_L[1];
    temp3 = P_j_dot - K * gamma(i) * (P_F(all, i) - P_F_Last(all, j));
    temp4 = P_L_dot - K * gamma(i) * (P_F(all, i) - temp5);
    u_F(all, i) = (1 / temp2) * W(i, j) * temp3 + (1 / temp2) * temp1 * temp4;
}
void Coordinated::Interceptor::Check_Param() {
    bool condition_1 = (V_L < 0.000001 && phase > 1) || (State_L.z_or_alt < 0) || (abs(State_L.x_or_lat) >= 89.999999);
    self_destruct = condition_1;
}
void Coordinated::Interceptor::Track_L(double t_meas) {
    state_t x = { 0,0,0,0,0,0 };
    x[0] = deg_to_rad(State_L.x_or_lat);           /* Initial State_L values for the current iteration */
    x[1] = deg_to_rad(State_L.y_or_lon);
    x[2] = deg_to_rad(head_ang_L);
    x[3] = V_L;
    x[4] = P_L[0];
    x[5] = P_L[1];

    double t = State_L.time;
    double step_size = t_meas - t;
    RK4 Integrator;

    int last_phase = phase;
    phase =  (int)(R_L < (R_start + 10000)) + (int)(R_L < R_start) + (int)(R_L < 50000);
    if (last_phase == 0 && phase > 1) { phase = 1; }

    Check_Param();
    if (self_destruct) {
        return;
    }

    while (t < (t_meas - 0.000001))             /* While flight_time <= t_meas */
    {
        if (phase == 0) {
            t = t_meas;
            break;
        }
        else if (phase == 1) {
            if (save_csv) {
                Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), rad_to_deg(x[2]), 0, R_L, rad_to_deg(q_L), rad_to_deg(q_dot), rad_to_deg(lead_ang_L), rad_to_deg(lead_ang_T), 
                       rad_to_deg(a_T), rad_to_deg(head_ang_T) });
                valid_csv = true;
            }

            V_L = Initial_Velocity_L;
            head_ang_L = rad_to_deg(q_L);
            x[2] = deg_to_rad(head_ang_L);
            x[3] = V_L;
            t = t_meas;
            break;
        }
        else if (abs(State_L.time - t) < 0.000001) {
            MPNG();
        }
        if (save_csv)
        {
            Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), rad_to_deg(x[2]), command_L, R_L, rad_to_deg(q_L), rad_to_deg(q_dot), rad_to_deg(lead_ang_L), rad_to_deg(lead_ang_T),
                       rad_to_deg(a_T), rad_to_deg(head_ang_T) });
            valid_csv = true;
        }
        Integrator(Course_L(), x, t, step_size);        /* Solve the ODEs to find system states */   
        energy = sqrt(pow(x[4], 2) + pow(x[5], 2));
    }

    head_ang_L = rad_to_deg(x[2]);
    P_L[0] = x[4];
    P_L[1] = x[5];
    V_L = x[3];

    State_L.id = 0;
    State_L.time = t;
    State_L.y_or_lon = rad_to_deg(x[1]);
    State_L.x_or_lat = rad_to_deg(x[0]);
    if (phase > 1)
    {
        State_L.x_vel = V_L * sind(head_ang_L);
        State_L.y_vel = V_L * cosd(head_ang_L);
    }
    else {
        State_L.x_vel = 0;
        State_L.y_vel = 0;
    }
}
void Coordinated::Interceptor::Track_F(double t_meas) {
    int i;
    state_t x = { 0,0,0,0 };
    double t, step_size;
    bool condition = 0;
    RK4 Integrator;

    for (i = 0; i < n; i++)
    {
        x[0] = deg_to_rad(State_F[i].x_or_lat);           /* Initial State_L values for the current iteration */
        x[1] = deg_to_rad(State_F[i].y_or_lon);
        x[2] = P_F(0, i);
        x[3] = P_F(1, i);

        t = State_F[i].time;
        step_size = t_meas - t;

        bool condition = (State_F[i].z_or_alt < 0) || (abs(State_F[i].x_or_lat) >= 89.999999);
        self_destruct = condition;
        if (self_destruct) {
            return;
        }
        else if (phase == 0 || phase == 1) {
            t = t_meas;
            State_F[i].time = t;
            State_F[i].id = i + 1;
            continue;
        }

        if (save_csv && i == 0)
        {
            Rec2({ t, State_F[0].x_or_lat, State_F[0].y_or_lon, State_F[1].x_or_lat, State_F[1].y_or_lon, State_F[2].x_or_lat, State_F[2].y_or_lon, State_F[3].x_or_lat, State_F[3].y_or_lon });
            valid_csv = true;
        }
        CGS(i);
        if (self_destruct) {
            return;
        }
        head_ang_F = atan2(u_F(0, i), u_F(1, i));
        u_x = u_F(0, i);
        u_y = u_F(1, i);
        V_F = sqrt(pow(u_x, 2) + pow(u_y, 2));
        Integrator(Course_F(), x, t, step_size);        /* Solve the ODEs to find system states */

        State_F[i].id = i + 1;
        State_F[i].time = t;
        State_F[i].y_or_lon = rad_to_deg(x[1]);
        State_F[i].x_or_lat = rad_to_deg(x[0]);
        State_F[i].x_vel = u_F(0, i);
        State_F[i].y_vel = u_F(1, i);
        P_F(0, i) = x[2];
        P_F(1, i) = x[3];
    }
}

void Coordinated::Interceptor::First_Run(objectstate Target_State) {
    State_L.time = Target_State.time;
    int i;
    double q_F;
    for (i = 0; i < n; i++) {
        State_F[i].time = Target_State.time;
        Calc_Distance(State_F[i], Target_State, R_F[i], q_F);
    }
    Calc_Distance(State_L, Target_State, R_L, q_L);
    if (q_L < 0) { q_L = q_L + 360; }
    q_L = deg_to_rad(q_L);

    t_last = State_L.time;
    q_last = q_L;
    Last_Target_State = Target_State;

    is_first_run = false;
}
void Coordinated::Interceptor::Fine_Calculations_L(objectstate Last_Leader_State, objectstate Target_State) {
    double fine_time_steps = 2000, i = 0;
    double t1 = Last_Leader_State.time, t2 = Target_State.time;

    objectstate Fine_Target = Last_Target_State;
    double fine_target_lat = (Target_State.x_or_lat - Fine_Target.x_or_lat) / fine_time_steps;
    double fine_target_long = (Target_State.y_or_lon - Fine_Target.y_or_lon) / fine_time_steps;
    double fine_target_alt = (Target_State.z_or_alt - Fine_Target.z_or_alt) / fine_time_steps;

    objectstate Fine_Interceptor = Last_Leader_State;
    double fine_interceptor_lat = (State_L.x_or_lat - Fine_Interceptor.x_or_lat) / fine_time_steps;
    double fine_interceptor_long = (State_L.y_or_lon - Fine_Interceptor.y_or_lon) / fine_time_steps;
    double fine_interceptor_alt = (State_L.z_or_alt - Fine_Interceptor.z_or_alt) / fine_time_steps;


    for (i = t1 + (t2 - t1) / fine_time_steps; i < (t2 - 0.000001); i = i + (t2 - t1) / fine_time_steps) {
        Fine_Target.time = i;
        Fine_Target.x_or_lat = Fine_Target.x_or_lat + fine_target_lat;
        Fine_Target.y_or_lon = Fine_Target.y_or_lon + fine_target_long;
        Fine_Target.z_or_alt = Fine_Target.z_or_alt + fine_target_alt;

        Fine_Interceptor.time = i;
        Fine_Interceptor.x_or_lat = Fine_Interceptor.x_or_lat + fine_interceptor_lat;
        Fine_Interceptor.y_or_lon = Fine_Interceptor.y_or_lon + fine_interceptor_long;
        Fine_Interceptor.z_or_alt = Fine_Interceptor.z_or_alt + fine_interceptor_alt;

        Calc_Distance(Fine_Interceptor, Fine_Target, R_L, q_L);
        q_L = deg_to_rad(q_L + 360 * (q_L < 0));
        if (save_csv) {
            Rec1({ i, Fine_Interceptor.x_or_lat, Fine_Interceptor.y_or_lon, head_ang_L, command_L, R_L, rad_to_deg(q_L), rad_to_deg(q_dot), rad_to_deg(lead_ang_L), rad_to_deg(lead_ang_T),
                       rad_to_deg(a_T), rad_to_deg(head_ang_T) });
            valid_csv = true;
        }

        my_terminated = Is_Hit() || self_destruct;
        my_killed = Is_Hit();
        if (my_terminated)
        {
            /* Extract variables required to calculate metrics */
            if (R_L < dist_LOS || dist_LOS < 0)
            {
                dist_LOS = R_L;
                if (R_F[Find_Min(R_F)] < R_L)
                {
                    dist_LOS = R_F[Find_Min(R_F)];
                }
                dist_covered = energy;
                time_taken = State_L.time;
            }
            speed = V_L;
            altitude = State_L.z_or_alt;
            accuracy = 1;
            num = 1;

            break;
        }
    }
}
void Coordinated::Interceptor::Fine_Calculations_F(objectstate Last_State, objectstate Target_State, int i) {
    double fine_time_steps = 2000, j = 0;
    double t1 = Last_State.time, t2 = Target_State.time;
    double q_F;

    objectstate Fine_Target = Last_Target_State;
    double fine_target_lat = (Target_State.x_or_lat - Fine_Target.x_or_lat) / fine_time_steps;
    double fine_target_long = (Target_State.y_or_lon - Fine_Target.y_or_lon) / fine_time_steps;
    double fine_target_alt = (Target_State.z_or_alt - Fine_Target.z_or_alt) / fine_time_steps;

    objectstate Fine_Interceptor = Last_State;
    double fine_interceptor_lat = (State_F[i].x_or_lat - Fine_Interceptor.x_or_lat) / fine_time_steps;
    double fine_interceptor_long = (State_F[i].y_or_lon - Fine_Interceptor.y_or_lon) / fine_time_steps;
    double fine_interceptor_alt = (State_F[i].z_or_alt - Fine_Interceptor.z_or_alt) / fine_time_steps;


    for (j = t1 + (t2 - t1) / fine_time_steps; j < (t2 - 0.000001); j = j + (t2 - t1) / fine_time_steps) {
        Fine_Target.time = j;
        Fine_Target.x_or_lat = Fine_Target.x_or_lat + fine_target_lat;
        Fine_Target.y_or_lon = Fine_Target.y_or_lon + fine_target_long;
        Fine_Target.z_or_alt = Fine_Target.z_or_alt + fine_target_alt;

        Fine_Interceptor.time = j;
        Fine_Interceptor.x_or_lat = Fine_Interceptor.x_or_lat + fine_interceptor_lat;
        Fine_Interceptor.x_or_lat = Fine_Interceptor.y_or_lon + fine_interceptor_long;
        Fine_Interceptor.z_or_alt = Fine_Interceptor.z_or_alt + fine_interceptor_alt;

        Calc_Distance(Fine_Interceptor, Fine_Target, R_F[i], q_F);
        /*if (save_csv)
        {
            if (i == 0) { Rec2({ j, Fine_Interceptor.x_or_lat, Fine_Interceptor.x_or_lat, State_F[1].x_or_lat, State_F[1].y_or_lon }); }
            if (i == 1) { Rec2({ j, State_F[0].x_or_lat, State_F[0].y_or_lon, Fine_Interceptor.x_or_lat, Fine_Interceptor.x_or_lat }); }
            valid_csv = true;
        }*/

        my_terminated = Is_Hit() || self_destruct;
        my_killed = Is_Hit();
        if (my_terminated)
        {
            /* Extract variables required to calculate metrics */
            if (R_L < dist_LOS || dist_LOS < 0)
            {
                dist_LOS = R_L;
                if (R_F[Find_Min(R_F)] < dist_LOS)
                {
                    dist_LOS = R_F[Find_Min(R_F)];
                }
                dist_covered = energy;
                time_taken = State_L.time;
            }
            speed = V_L;
            altitude = State_L.z_or_alt;
            accuracy = 1;
            num = 1;

            break;
        }
    }
}
void Coordinated::Interceptor::Calc_Target_Overload(objectstate Target_State, double delta_t) {
    double head_ang_T_last = atan2(Last_Target_State.x_vel, Last_Target_State.y_vel);
    head_ang_T_last = head_ang_T_last + 2 * M_PI * (head_ang_T_last < 0);
    
    double rate_head_ang_T = (head_ang_T - head_ang_T_last) / delta_t;
    a_T = rate_head_ang_T * sqrt(pow(Last_Target_State.x_vel, 2) + pow(Last_Target_State.y_vel, 2));
    /* double A_x = (Target_State.x_vel - Last_Target_State.x_vel) / delta_T;
    double A_y = (Target_State.y_vel - Last_Target_State.y_vel) / delta_T;
    double A_z = (Target_State.z_vel - Last_Target_State.z_vel) / delta_T;
    double theta = atan2(A_x, A_y) - atan2(Target_State.x_vel, Target_State.y_vel);
    double A = sqrt(pow(A_x, 2) + pow(A_y, 2));
    a_T = A * sin(theta); */
}

void Coordinated::Interceptor::Update_Target_State(objectstate Target_State)
{
    int i;
    int last_phase;
    double t_meas = Target_State.time;
    double delta_t = t_meas - t_last;
    objectstate Last_Leader_State;
    vector<objectstate> Last_Follower_State;
    Last_Follower_State.resize(n);
    double delta_q;

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
    
    last_phase = phase;
    Last_Leader_State = State_L;
    Last_Follower_State = State_F;
    P_F_Last = P_F;
    Track_F(t_meas);
    Track_L(t_meas);

    if (R_L < 20000) {
        Fine_Calculations_L(Last_Leader_State, Target_State);
        for (i = 0; i < n; i++) {
            if (R_F[i] < 20000) {
                Fine_Calculations_F(Last_Follower_State[i], Target_State, i);
            }
        }
        if (my_terminated)
        {
            return;
        }
    }

    /* Calculation of Quantities required for Coordinated Guidance Command */
    Calc_Distance(State_L, Target_State, R_L, q_L);
    q_L = deg_to_rad(q_L + 360 * (q_L < 0));
    for (i = 0; i < n; i++) {
        double q_F;
        Calc_Distance(State_F[i], Target_State, R_F[i], q_F);
    }
    my_terminated = Is_Hit() || Is_Missed(last_phase) || self_destruct || (energy > energy_lim);
    my_killed = Is_Hit();

    lead_ang_L = deg_to_rad(head_ang_L) - q_L;
    head_ang_T = atan2(Target_State.x_vel, Target_State.y_vel);
    head_ang_T = head_ang_T + 2 * M_PI * (head_ang_T < 0);
    lead_ang_T = head_ang_T - q_L;

    if (delta_t > 0.000001) {
        //delta_q = (q_L - q_last) + ((q_L - q_last) > 6) * (-2 * M_PI) + ((q_L - q_last) < -6) * (2 * M_PI);
        delta_q = (q_L - q_last);
        q_dot = (delta_q) / delta_t;
        Calc_Target_Overload(Target_State, delta_t);
    }

    t_last = t_meas;
    q_last = q_L;
    Last_Target_State = Target_State;

    /* Extract variables required to calculate metrics */
    if (R_L < dist_LOS || dist_LOS < 0)
    {
        dist_LOS = R_L;
        if (R_F[Find_Min(R_F)] < R_L)
        {
            dist_LOS = R_F[Find_Min(R_F)];
        }
        dist_covered = energy;
        time_taken = State_L.time;
    }
    speed = V_L;
    altitude = State_L.z_or_alt;
    accuracy = 1;
    num = n+1;
};

bool Coordinated::Interceptor::Get_State(objectstate& xo_state) {
    bool valid_state = Is_Started() && !my_killed;
    objectstate State_Out;
    State_Out = State_L;
    State_Out.y_or_lon = deg_to_long(State_L.y_or_lon);
    State_Out.x_or_lat = deg_to_lat(State_L.x_or_lat);
    if (valid_state)
    {
        xo_state = State_Out;
    }

        /*if (xo_state.id == 0 || n == 0) {
            xo_state = State_L;
            return valid_state;
        }
        else if (xo_state.id > 0) {
            int temp;
            temp = xo_state.id;
            xo_state = State_F[temp-1];
        }*/

    return valid_state;
}
void Coordinated::Interceptor::Get_Record() {
    if (save_csv && valid_csv)
    {
        Rec1.csv("Output/Coordinated_Leader_States", { "Time", "Latitude", "Longitude", "Heading_Angle", "Normal_Acceleration_Command", "LOS_Distance_Rate", "LOS_Angle", "LOS_Angle_Rate", "Leading_Angle", "Leading_Angle_Target",
                  "Target_Overload", "Heading_Angle_Target"});
        Rec2.csv("Output/Coordinated_Follower_States", { "Time", "F_1_Latitude", "F_1_Longitude", "F_2_Latitude", "F_2_Longitude"});
    }

}

bool Coordinated::Interceptor::Is_Started() {
    return (phase > 0);
}
bool Coordinated::Interceptor::Is_Hit() {
    int i, count;
    count = (int)(R_L < R_stop);
    for (i = 0; i < n; i++) {
        count = count + (int)(R_F[i] < R_stop);
    }
    return (count >= 1);
}
bool Coordinated::Interceptor::Is_Missed(int last_phase) {
    //return (phase < last_phase);
    return 0;
}

bool Coordinated::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
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
