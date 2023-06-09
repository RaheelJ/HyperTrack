#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "Fractional_Interceptor.h"
#include <iostream>
#include "xml/pugixml.hpp"

void Fractional::Interceptor::Reset() {
    /* Variables to calculate metrics */
    dist_LOS = -1;
    speed = 0;
    altitude = 0;
    dist_covered = 0;
    time_taken = 0;
    accuracy = 1;
    num = 1;

    command = 0;
    head_ang = 0;
    R = 0;
    q = 0;
    phase = 0;
    x_dist = 0;
    y_dist = 0;
    energy = 0;

    q_dot = 0;
    V_R = 0;
    V_M = 0;
    q_dot_t.assign(size, 0);
    interval_t.assign(size, 0);
    t_last = 0;
    q_last = 0;
    R_last = 0;
    GL_q_dot = 0;
        
    Last_Target_State.id = -1;
    my_terminated = false;
    my_killed = false;
    self_destruct = false;
    valid_csv = false;
    is_first_run = true;
}
bool Fractional::Interceptor::Initialize(const string& XML_Parameters, string& message) {
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
    pugi::xml_node Algo = Root.child("Fractional");
    if (Root && Common && Algo) {}
    else {
        message = "Root Node (Parameters) not found!";
        return false;
    }

    pugi::xml_node Position = Common.child("Position");
    pugi::xml_node Boost_Data = Common.child("Boost_Data");
    pugi::xml_node Exec_Param = Root.child("Selection_Execution");
    if (Position && Boost_Data && Exec_Param) {}
    else {
        message = "One of the child nodes is missing!";
        return false;
    }

    /* Initialize parameters from "Position" node */
    Initial_State.x_or_lat = Read_Param(Position.attribute("Latitude").as_double(), -90, 90, error);
    Initial_State.y_or_lon = Read_Param(Position.attribute("Longitude").as_double(), -180, 180, error);
    Initial_State.z_or_alt = Read_Param(Position.attribute("Altitude").as_double(), 0, 6378100, error);

    /* Initialize parameters from "FCG_Data" node */
    order = Read_Param(stod(Algo.child_value("Fractional_Order")), 0, 10, error);
    size = (int)Read_Param(stod(Algo.child_value("Size")), 1, 20, error);
    K_P = Read_Param(stod(Algo.child_value("Proportional_Gain")), 0, 100, error);
    K_D = Read_Param(stod(Algo.child_value("Differential_Gain")), 0, 100, error);
    Initial_Velocity = Read_Param(stod(Boost_Data.child_value("Boost_Speed")), 0, 10000, error);
    Re = 6378100;

    /* Initialize parameters from "Execution" node */
    R_start = Read_Param(stod(Exec_Param.child_value("Start_Distance")), 0, 10000000, error);
    R_stop = Read_Param(stod(Exec_Param.child_value("Stop_Distance")), 0, 5000, error);
    energy_lim = Read_Param(stod(Exec_Param.child_value("Energy_Limit")), 0, 10000000, error);
    //save_csv = Exec_Param.attribute("Save_CSV").as_bool();

    /* Check for out of range parameters */
    if (error == 1) {
        message = to_string(error) + " parameter is out of range";
        return false;
    }
    else if (error > 1) {
        message = to_string(error) + " parameters are out of range";
        return false;
    }

    q_dot_t.resize(size, 0);
    interval_t.resize(size, 0);

    Initial_State.id = -1;
    Initial_State.time = 0;
    Initial_State.state_in_geo = true;
    State = Initial_State;
    radial_distance = State.z_or_alt + Re;
    
    Reset();
    return true;
};
void Fractional::Interceptor::Reinitialize() {
    State = Initial_State;
    Reset();
}

void Fractional::Interceptor::FCG() {
    if (q_dot_t.size() < size) {
        self_destruct = true;
        return;
    }
    int n;
    GL_q_dot = 0;
    for (n = 0; n < size; n++)
    {
        if (interval_t[n] > 0.000001) {
            GL_q_dot = GL_q_dot + tgamma(-order + 1) / (tgamma(n) * tgamma(-order + n + 1)) * q_dot_t[n] / interval_t[n];
        }
    }
    command = K_P * abs(V_R) * q_dot + K_D * abs(V_R) * GL_q_dot;
    Limit(command, 10 * 9.8, -10 * 9.8);
}
void Fractional::Interceptor::Course::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] -- Latitude (radians)
        x[1] -- Longitude (radians)
        x[2] -- Velocity Heading Angle of Interceptor (radians)
        x[3] -- Speed 
        x[4] -- x distance
        x[5] -- y distance */

    dx[0] = x[3] * cos(x[2]) / radial_distance;
    dx[1] = x[3] * sin(x[2]) / (radial_distance * cos(x[0]));
    dx[2] = command / x[3];
    dx[3] = 0;
    dx[4] = x[3] * sin(x[2]);
    dx[5] = x[3] * cos(x[2]);
}
void Fractional::Interceptor::Check_Param() {
    bool condition_1 = (V_M < 0.000001 && phase > 1) || (State.z_or_alt < 0) || (abs(State.x_or_lat) >= 89.999999);
    self_destruct = condition_1;
}
void Fractional::Interceptor::Track(double t_meas) {
    state_t x = { 0,0,0,0,0,0 };
    x[0] = deg_to_rad(State.x_or_lat);           /* Initial state values for the current iteration */
    x[1] = deg_to_rad(State.y_or_lon);
    x[2] = deg_to_rad(head_ang);
    x[3] = V_M;
    x[4] = x_dist;
    x[5] = y_dist;

    double t = State.time;
    double step_size = t_meas - t;
    RK4 Integrator;

    int last_phase = phase;
    phase =  (int)(R < (R_start + 10000)) + (int)(R < R_start) + (int)(R < 50000);
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
                Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), rad_to_deg(x[2]), 0, R, V_R, rad_to_deg(q), rad_to_deg(q_dot), rad_to_deg(GL_q_dot) });
                valid_csv = true;
            }

            V_M = Initial_Velocity;
            head_ang = rad_to_deg(q);
            x[2] = deg_to_rad(head_ang);
            x[3] = V_M;
            t = t_meas;
            break;
        }
        else if (abs(State.time - t) < 0.000001) {
            FCG();
        }
        if (save_csv)
        {
            Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), rad_to_deg(x[2]), rad_to_deg(command / V_M), R, V_R, rad_to_deg(q), rad_to_deg(q_dot), rad_to_deg(GL_q_dot) });
            valid_csv = true;
        }
        Integrator(Course(), x, t, step_size);        /* Solve the ODEs to find system states */   
        energy = sqrt(pow(x[4], 2) + pow(x[5], 2));
    }

    head_ang = rad_to_deg(x[2]);
    x_dist = x[4];
    y_dist = x[5];
    V_M = x[3];

    State.id = 0;
    State.time = t;
    State.y_or_lon = rad_to_deg(x[1]);
    State.x_or_lat = rad_to_deg(x[0]);
    if (phase > 1)
    {
        State.y_vel = V_M * cosd(head_ang);
        State.x_vel = V_M * sind(head_ang);
    }
    else {
        State.x_vel = 0;
        State.y_vel = 0;
    }
}

void Fractional::Interceptor::First_Run(objectstate Target_State) {
    State.time = Target_State.time;
    Calc_Distance(State, Target_State, R, q);
    if (q < 0) { q = q + 360; }
    q = deg_to_rad(q);

    t_last = State.time;
    q_last = q;
    R_last = R;
    Last_Target_State = Target_State;

    is_first_run = false;
}
void Fractional::Interceptor::Fine_Calculations(objectstate Last_Interceptor_State, objectstate Target_State) {
    double fine_time_steps = 100, i = 0;
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

        Calc_Distance(Fine_Interceptor, Fine_Target, R, q);
        q = deg_to_rad(q + 360 * (q < 0));
        if (save_csv) {
            Rec1({ i, Fine_Interceptor.x_or_lat, Fine_Interceptor.y_or_lon, head_ang, rad_to_deg(command / V_M), R, V_R, rad_to_deg(q), rad_to_deg(q_dot), rad_to_deg(GL_q_dot) });
            valid_csv = true;
        }

        my_terminated = Is_Hit() || self_destruct;
        my_killed = Is_Hit();
        if (my_terminated)
        {
            /* Extract variables required to calculate metrics */
            if (R < dist_LOS || dist_LOS < 0)
            {
                dist_LOS = R;
                dist_covered = energy;
                time_taken = State.time;
            }
            speed = V_M;
            altitude = State.z_or_alt;
            accuracy = 1;
            num = 1;

            break;
        }
    }
}

void Fractional::Interceptor::Update_Target_State(objectstate Target_State)
{
    int last_phase;
    double t_meas = Target_State.time;
    double delta_t = t_meas - t_last;
    objectstate Last_Interceptor_State;
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
    Last_Interceptor_State = State;
    Track(t_meas);

    if (R < 20000) {
        Fine_Calculations(Last_Interceptor_State, Target_State);
        if (my_terminated)
        {
            return;
        }
    }

    /* Calculation of Quantities required for Fractional Guidance Command */
    Calc_Distance(State, Target_State, R, q);
    q = deg_to_rad(q + 360 * (q < 0));
    my_terminated = Is_Hit() || Is_Missed(last_phase) || self_destruct || (energy > energy_lim);
    my_killed = Is_Hit();

    if (delta_t > 0.000001) {
        delta_q = (q - q_last) + ((q - q_last) > 6) * (-2 * M_PI) + ((q - q_last) < -6) * (2 * M_PI);
        q_dot = (delta_q) / delta_t;
        q_dot_t.insert(q_dot_t.begin(), q_dot);
        q_dot_t.resize(size);
        V_R = (R - R_last) / delta_t;
        interval_t.insert(interval_t.begin(), delta_t);
        interval_t.resize(size);
    }

    t_last = t_meas;
    q_last = q;
    R_last = R;
    Last_Target_State = Target_State;

    /* Extract variables required to calculate metrics */
    if (R < dist_LOS || dist_LOS < 0)
    {
        dist_LOS = R;
        dist_covered = energy;
        time_taken = State.time;
    }
    speed = V_M;
    altitude = State.z_or_alt;
    accuracy = 1;
    num = 1;
};

bool Fractional::Interceptor::Get_State(objectstate& xo_state) {
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
void Fractional::Interceptor::Get_Record() {
    if (save_csv && valid_csv)
    {
        Rec1.csv("Output/Fractional_States", { "Time", "Latitude", "Longitude", "Heading_Angle", "Normal_Acceleration_Command", "LOS_Distance", "LOS_Distance_Rate", "LOS_Angle", "LOS_Angle_Rate", "LOS_Angle_Acceleration"});
    }

}

bool Fractional::Interceptor::Is_Started() {
    return (phase > 0);
}
bool Fractional::Interceptor::Is_Hit() {
    return (R < R_stop);
}
bool Fractional::Interceptor::Is_Missed(int last_phase) {
    //return (phase < last_phase);
    return 0;
}

bool Fractional::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
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
