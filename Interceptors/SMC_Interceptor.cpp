#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "SMC_Interceptor.h"
#include <iostream>
#include "xml/pugixml.hpp"
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>

using namespace Eigen;

void SMC::Interceptor::Reset() {
    /* Variables to calculate metrics */
    dist_LOS = -1;
    speed = 0;
    altitude = 0;
    dist_covered = 0;
    time_taken = 0;
    accuracy = 1;
    num = 1;

    /* Static Variable */
    command_lat = 0;
    command_lon = 0;
    head_ang = 0;
    path_ang = 0;

    // Interceptor States
    R = 0;
    q_lon = 0;
    q_lat = 0;
    lead_lat_m = 0;
    lead_lon_m = 0;
    phase = 0;
    x_dist = 0;
    y_dist = 0;
    energy = 0;

    /* Target States */
    Last_Target_State.id = -1;
    lead_lat_t = 0;
    lead_lon_t = 0;

    // Measurement Data
    q_lat_dot = 0; 
    q_lon_dot = 0;
        
    // Execution Flags
    my_terminated = false;
    my_killed = false;
    self_destruct = false;
    valid_csv = false;
    is_first_run = true;
}
bool SMC::Interceptor::Initialize(const string& XML_Parameters, string& message) {
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
    pugi::xml_node Algo = Root.child("SMC");
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
    Initial_Velocity = Read_Param(stod(Boost_Data.child_value("Boost_Speed")), 0, 10000, error);
    Initial_Altitude = Read_Param(stod(Boost_Data.child_value("Boost_Altitude")), 0, 1000000, error);

    /* Initialize Parameters for Guidance Law */
    alpha = Read_Param(stod(Algo.child_value("alpha")), -100, 100, error);
    p = Read_Param(stod(Algo.child_value("p")), -100, 100, error);
    r0 = Read_Param(stod(Algo.child_value("r0")), -100, 100, error);
    r1 = Read_Param(stod(Algo.child_value("r1")), -100, 100, error);
    r2 = Read_Param(stod(Algo.child_value("r2")), -100, 100, error);
    n1 = Read_Param(stod(Algo.child_value("n11")), -100, 100, error);
    n2 = Read_Param(stod(Algo.child_value("n22")), -100, 100, error);
    k1 = Read_Param(stod(Algo.child_value("k11")), -100, 100, error);
    k2 = Read_Param(stod(Algo.child_value("k22")), -100, 100, error);

    /* Initialize parameters from "Execution" node */
    R_start = Read_Param(stod(Exec_Param.child_value("Start_Distance")), 0, 10000000, error);
    R_stop = Read_Param(stod(Exec_Param.child_value("Stop_Distance")), 0, 5000, error);
    energy_lim = Read_Param(stod(Exec_Param.child_value("Energy_Limit")), 0, 10000000, error);
    //save_csv = Exec_Param.attribute("Save_CSV").as_bool();
    //save_csv = true;

    /* Check for out of range parameters */
    if (error == 1) {
        message = to_string(error) + " parameter is out of range";
        return false;
    }
    else if (error > 1) {
        message = to_string(error) + " parameters are out of range";
        return false;
    }

    Initial_State.id = -1;
    Initial_State.time = 0;
    Initial_State.state_in_geo = true;
    State = Initial_State;
    
    Reset();
    return true;
};
void SMC::Interceptor::Reinitialize() {
    State = Initial_State;
    Reset();
}

void SMC::Interceptor::SMCGuidance() {
    Matrix2d B, N, K, Q;
    Vector2d E, F, S, Y, NS, U;
    double Ts = 4e-6;           //Sampling Time for the Guidance Law

    B <<    1 / Vm,     0,
            0,          1 / (Vm * cos(lead_lon_m));
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    N <<    n1,         0,
            0,          n2;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    K <<    k1,         0,
            0,          k2;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    F(0) = -q_lat_dot * sin(q_lon) * sin(lead_lat_m) - q_lon_dot * cos(lead_lat_m);
    F(1) = q_lat_dot * sin(q_lon) * cos(lead_lat_m) * tan(lead_lon_m) - q_lon_dot * sin(lead_lat_m) * tan(lead_lon_m) - q_lat_dot * cos(q_lon);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    E(0) = -q_lat_dot * sin(q_lon) * sin(lead_lat_t) - q_lon_dot * cos(lead_lat_t);
    E(1) = q_lat_dot * sin(q_lon) * cos(lead_lat_t) * tan(lead_lon_t) - q_lon_dot * sin(lead_lat_t) * tan(lead_lon_t) - q_lat_dot * cos(q_lon);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    S(0) = lead_lon_m - n1 * lead_lon_t;
    S(1) = lead_lat_m - n2 * lead_lat_t;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Q <<    signum(S(0)),   0,
            0,              signum(S(1));
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Y(0) = lead_lon_t;
    Y(1) = lead_lat_t;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    NS(0) = r0 * (exp(r1 * pow(abs(S(0)), p) - r2));
    NS(1) = r0 * (exp(r1 * pow(abs(S(1)), p) - r2));
    
    /* Guidance Law */
    U = -B.inverse() * Ts * (F - N * E + K * S + Q * (alpha * Y + NS));

    /* Formulating a Guidance Law for Long as well as Close Range */
    double g1 = 4;
    double g2 = 6;
    command_lon = (R <= 12000) * U(0) + (R > 12000) * -g1 * q_lon_dot * Vm;
    command_lat = (R <= 12000) * U(1) + (R > 12000) * g2 * q_lat_dot * Vm;
    Limit(command_lon, 10 * 9.8, -10 * 9.8);
    Limit(command_lat, 10 * 9.8, -10 * 9.8);
}
void SMC::Interceptor::Course::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] -- Latitude (radians)
        x[1] -- Longitude (radians)
        x[2] -- Velocity Heading Angle of Interceptor (radians)
        x[3] -- Speed 
        x[4] -- x distance
        x[5] -- y distance
        x[6] -- Altitude (m)
        x[7] -- Flight Path Angle (radians) */

    dx[0] = x[3] * cos(x[7]) * cos(x[2]) / (Re + x[6]);
    dx[1] = x[3] * cos(x[7]) * sin(x[2]) / ((Re + x[6]) * cos(x[0]));
    dx[2] = command_lat / (x[3] * cos(x[7]));
    dx[3] = 0;
    dx[4] = x[3] * cos(x[7]) * sin(x[2]);
    dx[5] = x[3] * cos(x[7]) * cos(x[2]);
    dx[6] = x[3] * sin(x[7]);
    dx[7] = command_lon / x[3];
}
void SMC::Interceptor::Check_Param() {
    bool condition1 = (Vm < 0.000001 && phase > 1) || (State.z_or_alt < 0) || (abs(State.x_or_lat) >= 89.999999);
    bool condition2 = (abs(lead_lon_m) >= 89.999999) || (abs(lead_lon_t) >= 89.999999);
    self_destruct = condition1 || condition2;
}
void SMC::Interceptor::Track(double t_meas) {
    state_t x = { 0,0,0,0,0,0,0,0 };
    x[0] = deg_to_rad(State.x_or_lat);           /* Initial state values for the current iteration */
    x[1] = deg_to_rad(State.y_or_lon);
    x[2] = deg_to_rad(head_ang);
    x[3] = Vm;
    x[4] = x_dist;
    x[5] = y_dist;
    x[6] = State.z_or_alt;
    x[7] = deg_to_rad(path_ang);

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
                Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[6], head_ang, path_ang, rad_to_deg(q_lat), rad_to_deg(q_lon), command_lat, command_lon });
                valid_csv = true;
            }
            State.z_or_alt = Initial_Altitude;
            Vm = Initial_Velocity;
            head_ang = rad_to_deg((q_lat >= M_PI) * (q_lat - M_PI) + (q_lat < M_PI) * (q_lat + M_PI));
            path_ang = rad_to_deg(-q_lon);

            x[2] = deg_to_rad(head_ang);
            x[7] = deg_to_rad(path_ang);
            x[3] = Vm;
            x[6] = State.z_or_alt;
            t = t_meas;
            break;
        }
        else if (abs(State.time - t) < 0.000001) {
            SMCGuidance();
        }
        if (save_csv)
        {
            Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[6], head_ang, path_ang, rad_to_deg(q_lat), rad_to_deg(q_lon), command_lat, command_lon });
            valid_csv = true;
        }
        Integrator(Course(), x, t, step_size);        /* Solve the ODEs to find system states */   
        energy = sqrt(pow(x[4], 2) + pow(x[5], 2));
    }

    head_ang = rad_to_deg(x[2]);
    path_ang = rad_to_deg(x[7]);
    x_dist = x[4];
    y_dist = x[5];
    Vm = x[3];

    State.id = 0;
    State.time = t;
    State.z_or_alt = x[6];
    State.y_or_lon = rad_to_deg(x[1]);
    State.x_or_lat = rad_to_deg(x[0]);
    if (phase > 1)
    {
        State.y_vel = Vm * cosd(path_ang) * cosd(head_ang);
        State.x_vel = Vm * cosd(path_ang) * sind(head_ang);
        State.z_vel = Vm * sind(path_ang);
    }
    else {
        State.x_vel = 0;
        State.y_vel = 0;
        State.z_vel = 0;
    }
}

void SMC::Interceptor::First_Run(objectstate Target_State) {
    State.time = Target_State.time;
    
    /* Calculate LOS Angle in Lateral and Longitudnal Plane */
    Calc_LOS(Target_State, State, R, q_lat, q_lon);
    lead_lat_m = deg_to_rad(head_ang) - q_lat;
    lead_lon_m = deg_to_rad(path_ang) - q_lon;

    double head_ang_t = Calc_Head_Ang(Target_State);
    double path_ang_t = Calc_Path_Ang(Target_State);
    lead_lat_t = head_ang_t - q_lat;
    lead_lon_t = path_ang_t - q_lon;

    Last_Target_State = Target_State;
    is_first_run = false;
}
void SMC::Interceptor::Fine_Calculations(objectstate Last_Interceptor_State, objectstate Target_State) {
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

        Calc_LOS(Fine_Target, Fine_Interceptor, R, q_lat, q_lon);
        if (save_csv) {
            Rec1({ i, Fine_Interceptor.x_or_lat, Fine_Interceptor.y_or_lon, Fine_Interceptor.z_or_alt, head_ang, path_ang, rad_to_deg(q_lat), rad_to_deg(q_lon), command_lat, command_lon });
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
            speed = Vm;
            altitude = State.z_or_alt;
            accuracy = 1;
            num = 1;

            break;
        }
    }
}

void SMC::Interceptor::Update_Target_State(objectstate Target_State)
{
    int last_phase;
    objectstate Last_Interceptor_State;
    double delta_q_lat, delta_q_lon;
    double t_last = State.time;
    double q_lat_last = q_lat;
    double q_lon_last = q_lon;
    double R_last = R;
        
    double t_meas = Target_State.time;
    double delta_t = t_meas - t_last;

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

    /* Calculation of Quantities required for SMC Guidance Command */
    Calc_LOS(Target_State, State, R, q_lat, q_lon);
    lead_lat_m = deg_to_rad(head_ang) - q_lat;
    lead_lon_m = deg_to_rad(path_ang) - q_lon;

    double head_ang_t = Calc_Head_Ang(Target_State);
    double path_ang_t = Calc_Path_Ang(Target_State);
    lead_lat_t = head_ang_t - q_lat;
    lead_lon_t = path_ang_t - q_lon;

    /* Check the Termination Conditions */
    my_terminated = Is_Hit() || Is_Missed(last_phase) || self_destruct || (energy > energy_lim);
    my_killed = Is_Hit();

    if (delta_t > 0.000001) {
        delta_q_lat = (q_lat - q_lat_last) + ((q_lat - q_lat_last) > 6) * (-2 * M_PI) + ((q_lat - q_lat_last) < -6) * (2 * M_PI);
        q_lat_dot = (delta_q_lat) / delta_t;
        delta_q_lon = (q_lon - q_lon_last);
        q_lon_dot = (delta_q_lon) / delta_t;
    }

    Last_Target_State = Target_State;

    /* Extract variables required to calculate metrics */
    if (R < dist_LOS || dist_LOS < 0)
    {
        dist_LOS = R;
        dist_covered = energy;
        time_taken = State.time;
    }
    speed = Vm;
    altitude = State.z_or_alt;
    accuracy = 1;
    num = 1;
};

bool SMC::Interceptor::Get_State(objectstate& xo_state) {
    bool valid_state = Is_Started() && !my_killed;
    objectstate State_Out;
    State_Out = State;
    State_Out.y_or_lon = deg_to_long(State.y_or_lon);
    State_Out.x_or_lat = deg_to_lat(State.x_or_lat);
    State_Out.z_or_alt = State.z_or_alt;
    if (valid_state)
    {
        xo_state = State_Out;
    }
    return valid_state;
}
void SMC::Interceptor::Get_Record() {
    if (save_csv && valid_csv)
    {
        Rec1.csv("Output/SMC_States", { "Time", "Latitude", "Longitude", "Altitude", "Heading_Angle", "Path_Angle", "LOS_Angle_Lat", "LOS_Angle_Lon", "Command_Lat", "Command_Lon"});
    }

}

bool SMC::Interceptor::Is_Started() {
    return (phase > 0);
}
bool SMC::Interceptor::Is_Hit() {
    return (R < R_stop);
}
bool SMC::Interceptor::Is_Missed(int last_phase) {
    //return (phase < last_phase);
    return 0;
}

bool SMC::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
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
