#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "Coordinated_TPNG_Interceptor.h"
#include <iostream>
#include "xml/pugixml.hpp"
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>

using namespace Eigen;

void Coordinated_TPNG::Interceptor::Reset() {
    /* Variables to calculate metrics */
    dist_LOS = -1;
    speed = 0;
    altitude = 0;
    dist_covered = 0;
    time_taken = 0;
    accuracy = 1;
    num = 1;

    /* Static Variable */
    aM_lon = 0;
    aM_lat = 0;

    // Interceptor States
    LOS_dist.resize(num_M);
    LOS_dist.resize(num_M);
    LOS_ang_lat.resize(num_M);
    LOS_ang_lon.resize(num_M);
    phase.resize(num_M);
    x_dist.resize(num_M);
    y_dist.resize(num_M);
    energy.resize(num_M);
    head_ang.resize(num_M);
    path_ang.resize(num_M);
    Rec.resize(num_M);

    /* Target States */
    Last_Target_State.id = -1;

    // Measurement Data
    LOS_ang_lat_dot.resize(num_M);
    LOS_ang_lon_dot.resize(num_M);
    LOS_dist_dot.resize(num_M);
    ZEM_0.resize(num_M);
    ZEM.resize(num_M);
    delta_aZEM.resize(num_M);

    // Parameters of Cooperative Coverage Strategy
    a_M_max = 0;
    aM_max = 0;

    // Parameters of TPNG Law 
    command_lat.resize(num_M);
    command_lon.resize(num_M);
    t_go.resize(num_M);

    // Execution Flags
    my_terminated = false;
    my_killed = false;
    self_destruct = false;
    valid_csv = false;
    is_first_run = true;
}
bool Coordinated_TPNG::Interceptor::Initialize(const string& XML_Parameters, string& message) {
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
    pugi::xml_node Algo = Root.child("TPNG");
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

    /* Initialize parameters from "Boost_Data" node */
    Initial_Velocity = Read_Param(stod(Boost_Data.child_value("Boost_Speed")), 0, 10000, error);
    Initial_Altitude = Read_Param(stod(Boost_Data.child_value("Boost_Altitude")), 0, 1000000, error);

    /* Initialize Parameters for Guidance Law */
    N = Read_Param(stod(Algo.child_value("N")), 0.2, 100, error);
    num_M_max = Read_Param(stod(Algo.child_value("num_M_max")), 0, 1, error);
    aT_max = Read_Param(stod(Algo.child_value("aT_max")), 2, 20, error);
    cov_ratio = Read_Param(stod(Algo.child_value("cov_ratio")), 0.34, 100, error);
    
    /* Initialize parameters from "Execution" node */
    dist_start = Read_Param(stod(Exec_Param.child_value("Start_Distance")), 0, 10000000, error);
    dist_stop = Read_Param(stod(Exec_Param.child_value("Stop_Distance")), 0, 5000, error);
    energy_lim = Read_Param(stod(Exec_Param.child_value("Energy_Limit")), 0, 10000000, error);
    //save_csv = Exec_Param.attribute("Save_CSV").as_bool();
    //save_csv = true;

    /* Initialize parameters from "Position" node */
    int i = 0;
    Initial_State.x_or_lat = Read_Param(Position.attribute("Latitude").as_double(), -90, 90, error);
    Initial_State.y_or_lon = Read_Param(Position.attribute("Longitude").as_double(), -180, 180, error);
    Initial_State.z_or_alt = Read_Param(Position.attribute("Altitude").as_double(), 0, 6378100, error);
    Initial_State.id = -1;
    Initial_State.time = 0;
    Initial_State.state_in_geo = true;

    Coop_Coverage(num_M_max);
    for (i = 0; i < num_M; i++)
    {
        State.push_back(Initial_State);
    }
    
    /* Check for out of range parameters */
    if (error == 1) {
        message = to_string(error) + " parameter is out of range";
        return false;
    }
    else if (error > 1) {
        message = to_string(error) + " parameters are out of range";
        return false;
    }
    
    Reset();
    return true;
};
void Coordinated_TPNG::Interceptor::Reinitialize() {
    int i = 0;
    for (i = 0; i < num_M; i++)
    {
        State[i] = Initial_State;
    }
    Reset();
}

double Coordinated_TPNG::Interceptor::sin_DOmF(int n)
{
    double A, B, C, Out;
    A = sin(n);
    B = (pow(cov_ratio, 2) - 1) / (4 * pow(cov_ratio, 2) * cos(M_PI / n));
    C = sin(acos(B + 1 / (2 * cov_ratio)));

    Out = A / C;
    return Out;
}
void Coordinated_TPNG::Interceptor::Coop_Coverage(int n_max)
{
    double ACB;
    double scale_factor, lon, lat;
    int i, test_num;
    int step;           /* Step Number of the Algorithm */
    step = 2 * (cov_ratio >= (sqrt(2) / 2)) + 3 * ((cov_ratio < (sqrt(2) / 2)) && (cov_ratio >= 0.577)) + 4 * (cov_ratio < 0.577);
    switch (step)
    {
    case 2:
        num_M = round(M_PI / asin(cov_ratio));
        CoC.resize(num_M);
        for (i = 1; i <= num_M; i++)
        {
            scale_factor = cos(M_PI / num_M) - sqrt(pow(cov_ratio, 2) - pow(sin(M_PI / num_M), 2));
            lat = cos(2 * M_PI * i / num_M);
            lon = sin(2 * M_PI * i / num_M);

            CoC[i - 1].lat = scale_factor * lat;
            CoC[i - 1].lon = scale_factor * lon;
        }
        break;
    case 3:
        ACB = 2 * acos(1 / (2 * cov_ratio));
        num_M = round(2 * M_PI / (ACB));
        CoC.resize(num_M);
        for (i = 1; i <= num_M; i++)
        {
            scale_factor = cos(M_PI / num_M) - sqrt(pow(cov_ratio, 2) - pow(sin(M_PI / num_M), 2));
            lat = cos(2 * M_PI * i / num_M);
            lon = sin(2 * M_PI * i / num_M);

            CoC[i - 1].lat = scale_factor * lat;
            CoC[i - 1].lon = scale_factor * lon;
        }
        break;
    case 4:
        for (test_num = 5; test_num <= n_max; test_num++)
        {
            if ((cov_ratio - sin_DOmF(test_num) < 0.04) && (cov_ratio - sin_DOmF(test_num) > -0.04))
            {
                num_M = test_num;
                break;
            }
        }
        if (num_M < 5)
        {
            num_M = n_max;
        }
        CoC.resize(num_M);
        CoC[0].lat = 0;
        CoC[0].lon = 0;
        for (i = 1; i <= num_M - 1; i++)
        {
            scale_factor = 2 * cov_ratio * cos(M_PI / (num_M - 1));
            lat = cos(2 * M_PI * i / (num_M - 1));
            lon = sin(2 * M_PI * i / (num_M - 1));

            CoC[i].lat = scale_factor * lat;
            CoC[i].lon = scale_factor * lon;
        }
        break;
    default:
        break;
    }
    /* Number of Interceptors Needed */
    // num_M = (N-2)*cov_ratio/N;
}

Point_3D Coordinated_TPNG::Interceptor::Calc_ZEM(int index)
{
    int i = index;
    Point_3D Wi, Ri, Wi_Ri;
    Point_3D ZEMi;

    t_go[i] = LOS_dist[i] / abs(LOS_dist_dot[i]);

    Wi.x = LOS_ang_lat_dot[i] * sin(LOS_ang_lon[i]);
    Wi.y = LOS_ang_lat_dot[i] * cos(LOS_ang_lon[i]);
    Wi.z = LOS_ang_lon_dot[i];
    Ri.x = LOS_dist[i] * cos(LOS_ang_lon[i]) * cos(LOS_ang_lat[i]);
    Ri.y = LOS_dist[i] * sin(LOS_ang_lon[i]);
    Ri.z = LOS_dist[i] * cos(LOS_ang_lon[i]) * -sin(LOS_ang_lat[i]);
    Wi_Ri = CrossProduct(Wi, Ri);

    ZEMi.x = Wi_Ri.x * t_go[i];
    ZEMi.y = Wi_Ri.y * t_go[i];
    ZEMi.z = Wi_Ri.z * t_go[i];

    return ZEMi;
}
void Coordinated_TPNG::Interceptor::Calc_Delta_aZEM(int index)
{
    int i = index;
    double norm_aZEM, norm_aZEM_0;
    Point_3D aZEM, aZEM_0;

    aZEM.x = -2 * ZEM[i].x / pow(t_go[i], 2);
    aZEM.y = -2 * ZEM[i].y / pow(t_go[i], 2);
    aZEM.z = -2 * ZEM[i].z / pow(t_go[i], 2);
    aZEM_0.x = -2 * ZEM_0[i].x / pow(t_go[i], 2);
    aZEM_0.y = -2 * ZEM_0[i].y / pow(t_go[i], 2);
    aZEM_0.z = -2 * ZEM_0[i].z / pow(t_go[i], 2);
    norm_aZEM = sqrt(pow(aZEM.x, 2) + pow(aZEM.y, 2) + pow(aZEM.y, 2));
    norm_aZEM_0 = sqrt(pow(aZEM_0.x, 2) + pow(aZEM_0.y, 2) + pow(aZEM_0.y, 2));
    delta_aZEM[i] = norm_aZEM - norm_aZEM_0;
}
void Coordinated_TPNG::Interceptor::Calc_AIM(int index)
{
    int i=index;
    Point_2D a_bias;

    a_bias = CoC[i];
    
    AIM[i].x = 0;
    AIM[i].y = a_bias.lon * pow(t_go[i], 2) / N;
    AIM[i].z = a_bias.lat * pow(t_go[i], 2) / N;

    //AIM[i].mag = sqrt(pow(CoC[i].lat, 2) + pow(CoC[i].lon, 2));
    //AIM[i].ang = atan2(CoC[i].lon, CoC[i].lat);
}
void Coordinated_TPNG::Interceptor::Adjust_AIM(int index)
{
    int i = index;

    Calc_Delta_aZEM(i);
    int i_min = Find_Min(delta_aZEM);
    double K = exp(delta_aZEM[i_min] - delta_aZEM[i]);

    AIM[i].y = AIM[i].y * K + AIM[i_min].y * (1 - K);
    AIM[i].z = AIM[i].z * K + AIM[i_min].z * (1 - K);
}
void Coordinated_TPNG::Interceptor::TPNG(int index) 
{
    int i=index;

    ZEM[i] = Calc_ZEM(i);
    command_lon[i] = N * (ZEM[i].y - AIM[i].y) / pow(t_go[i], 2);
    command_lat[i] = N * (ZEM[i].z - AIM[i].z) / pow(t_go[i], 2);

    Adjust_AIM(i);    
}
void Coordinated_TPNG::Interceptor::Course::operator()(const state_t& x, state_t& dx, const double) {
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
    dx[2] = aM_lat / (x[3] * cos(x[7]));
    dx[3] = 0;
    dx[4] = x[3] * cos(x[7]) * sin(x[2]);
    dx[5] = x[3] * cos(x[7]) * cos(x[2]);
    dx[6] = x[3] * sin(x[7]);
    dx[7] = aM_lon / x[3];
}
void Coordinated_TPNG::Interceptor::Check_Param(int index) {
    bool condition1 = false;
    condition1 = condition1 || (speed_M < 0.000001 && phase[index] > 1) || (State[index].z_or_alt < 0) || (abs(State[index].x_or_lat) >= 89.999999);
    self_destruct = condition1;
}
void Coordinated_TPNG::Interceptor::Track(double t_meas, int index) {
    int i, last_phase;
    state_t x = { 0,0,0,0,0,0,0,0 };
    double t, step_size;
    RK4 Integrator;

    i = index;
    x[0] = deg_to_rad(State[i].x_or_lat);           /* Initial state values for the current iteration */
    x[1] = deg_to_rad(State[i].y_or_lon);
    x[2] = deg_to_rad(head_ang[i]);
    x[3] = speed_M;
    x[4] = x_dist[i];
    x[5] = y_dist[i];
    x[6] = State[i].z_or_alt;
    x[7] = deg_to_rad(path_ang[i]);

    t = State[0].time;
    step_size = t_meas - t;
        
    last_phase = phase[i];
    phase[i] = (int)((LOS_dist[i] < (dist_start + 10000)) || If_Found(phase, 1)) + (int)((LOS_dist[i] < dist_start) || If_Found(phase, 2)) + (int)(LOS_dist[i] < 50000);
    if (last_phase == 0 && phase[i] > 1) { phase[i] = 1; }

    Check_Param(i);
    if (self_destruct) {
        return;
    }

    while (t < (t_meas - 0.000001))             /* While flight_time <= t_meas */
    {
        if (phase[i] == 0) {
            t = t_meas;
            break;
        }
        else if (phase[i] == 1) {
            if (save_csv) {
                Rec[i]({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[6], head_ang[i], path_ang[i], rad_to_deg(LOS_ang_lat[i]), rad_to_deg(LOS_ang_lon[i]), command_lat[i], command_lon[i]});
                valid_csv = true;
            }
            ZEM_0[i] = Calc_ZEM(i);
            Calc_AIM(i);

            State[i].z_or_alt = Initial_Altitude;
            speed_M = Initial_Velocity;
            head_ang[i] = rad_to_deg(LOS_ang_lat[i]);
            path_ang[i] = rad_to_deg(LOS_ang_lon[i]);

            x[2] = deg_to_rad(head_ang[i]);
            x[7] = deg_to_rad(path_ang[i]);
            x[3] = speed_M;
            x[6] = State[i].z_or_alt;
            t = t_meas;
            break;
        }
        else if (abs(State[i].time - t) < 0.000001) {
            TPNG(i);
            aM_lon = command_lon[i];
            aM_lat = command_lat[i];
        }
        if (save_csv)
        {
            Rec[i]({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[6], head_ang[i], path_ang[i], rad_to_deg(LOS_ang_lat[i]), rad_to_deg(LOS_ang_lon[i]), command_lat[i], command_lon[i] });
            valid_csv = true;
        }
        Integrator(Course(), x, t, step_size);        /* Solve the ODEs to find system states */
        energy[i] = sqrt(pow(x[4], 2) + pow(x[5], 2));
    }
    head_ang[i] = rad_to_deg(x[2]);
    path_ang[i] = rad_to_deg(x[7]);
    x_dist[i] = x[4];
    y_dist[i] = x[5];
    speed_M = x[3];

    State[i].id = 0;
    State[i].time = t;
    State[i].z_or_alt = x[6];
    State[i].y_or_lon = rad_to_deg(x[1]);
    State[i].x_or_lat = rad_to_deg(x[0]);
    if (phase[i] > 1)
    {
        State[i].y_vel = speed_M * cosd(path_ang[i]) * cosd(head_ang[i]);
        State[i].x_vel = speed_M * cosd(path_ang[i]) * sind(head_ang[i]);
        State[i].z_vel = speed_M * sind(path_ang[i]);
    }
    else {
        State[i].x_vel = 0;
        State[i].y_vel = 0;
        State[i].z_vel = 0;
    }
}   


void Coordinated_TPNG::Interceptor::First_Run(objectstate Target_State) {
    int i;

    for (i = 0; i < num_M; i++)
    {
        /* Calculate LOS Angle in Lateral and Longitudnal Plane */
        Calc_LOS(State[i], Target_State, LOS_dist[i], LOS_ang_lat[i], LOS_ang_lon[i]);
        LOS_dist_lat[i] = LOS_dist[i] * cos(LOS_ang_lon[i]);

        State[i].time = Target_State.time;
    }

    Last_Target_State = Target_State;
    is_first_run = false;
}
void Coordinated_TPNG::Interceptor::Fine_Calculations(objectstate Last_Interceptor_State, objectstate Target_State, int index) {
    double fine_time_steps = 1000, i = 0;
    int j, best_j;
    double t1, t2;
    objectstate Fine_Target, Fine_Interceptor;
    double fine_target_lat, fine_target_long, fine_target_alt;
    double fine_interceptor_lat, fine_interceptor_long, fine_interceptor_alt;

    j = index;
    t1 = Last_Interceptor_State.time;
    t2 = Target_State.time;

    Fine_Target = Last_Target_State;
    fine_target_lat = (Target_State.x_or_lat - Fine_Target.x_or_lat) / fine_time_steps;
    fine_target_long = (Target_State.y_or_lon - Fine_Target.y_or_lon) / fine_time_steps;
    fine_target_alt = (Target_State.z_or_alt - Fine_Target.z_or_alt) / fine_time_steps;


    Fine_Interceptor = Last_Interceptor_State;
    fine_interceptor_lat = (State[j].x_or_lat - Fine_Interceptor.x_or_lat) / fine_time_steps;
    fine_interceptor_long = (State[j].y_or_lon - Fine_Interceptor.y_or_lon) / fine_time_steps;
    fine_interceptor_alt = (State[j].z_or_alt - Fine_Interceptor.z_or_alt) / fine_time_steps;

    for (i = t1 + (t2 - t1) / fine_time_steps; i < (t2 - 0.000001); i = i + (t2 - t1) / fine_time_steps) {
        Fine_Target.time = i;
        Fine_Target.x_or_lat = Fine_Target.x_or_lat + fine_target_lat;
        Fine_Target.y_or_lon = Fine_Target.y_or_lon + fine_target_long;
        Fine_Target.z_or_alt = Fine_Target.z_or_alt + fine_target_alt;

        Fine_Interceptor.time = i;
        Fine_Interceptor.x_or_lat = Fine_Interceptor.x_or_lat + fine_interceptor_lat;
        Fine_Interceptor.y_or_lon = Fine_Interceptor.y_or_lon + fine_interceptor_long;
        Fine_Interceptor.z_or_alt = Fine_Interceptor.z_or_alt + fine_interceptor_alt;

        Calc_LOS(Fine_Interceptor, Fine_Target, LOS_dist[j], LOS_ang_lat[j], LOS_ang_lon[j]);
        LOS_dist_lat[j] = LOS_dist[j] * cos(LOS_ang_lon[j]);
        if (save_csv) {
            Rec[j]({ i, Fine_Interceptor.x_or_lat, Fine_Interceptor.y_or_lon, Fine_Interceptor.z_or_alt, head_ang[j], path_ang[j], rad_to_deg(LOS_ang_lat[j]), rad_to_deg(LOS_ang_lon[j]), command_lat[j], command_lon[j]});
            valid_csv = true;
        }

        my_terminated = Is_Hit(j) || self_destruct;
        my_killed = Is_Hit(j);
        if (my_terminated)
        {
            best_j = Find_Min(LOS_dist);
            /* Extract variables required to calculate metrics */
            dist_LOS = LOS_dist[best_j];
            speed = speed_M;
            altitude = State[best_j].z_or_alt;
            dist_covered = energy[best_j];
            time_taken = State[best_j].time;
            accuracy = 1;
            num = num_M;

            break;
        }
    }

}

void Coordinated_TPNG::Interceptor::Update_Target_State(objectstate Target_State)
{
    int i, last_phase, i_min;
    objectstate Last_Interceptor_State;
    double delta_q_lat, delta_q_lon, delta_R, delta_R_lat;
    double t_last, q_lat_last, q_lon_last, R_last, R_lat_last;
    double t_meas, delta_t;

    t_meas = Target_State.time;
	if (t_meas - Last_Target_State.time <= 0)
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

    for (i = 0; i < num_M; i++)
    {
        t_last = State[i].time;
        q_lat_last = LOS_ang_lat[i];
        q_lon_last = LOS_ang_lon[i];
        R_last = LOS_dist[i];
        R_lat_last = LOS_dist_lat[i];

        delta_t = t_meas - t_last;

        last_phase = phase[i];
        Last_Interceptor_State = State[i];
        Track(t_meas, i);

        if (LOS_dist[i] < 20000) {
            Fine_Calculations(Last_Interceptor_State, Target_State, i);
            if (my_terminated)
            {
                return;
            }
        }

        /* Calculation of Quantities required for CPNG Command */
        Calc_LOS(State[i], Target_State, LOS_dist[i], LOS_ang_lat[i], LOS_ang_lon[i]);
        LOS_dist_lat[i] = LOS_dist[i] * cos(LOS_ang_lon[i]);

        /* Check the Termination Conditions */
        my_terminated = Is_Hit(i) || Is_Missed(last_phase, i) || self_destruct || (energy[i] > energy_lim);
        my_killed = Is_Hit(i);

        /* Calculate the Measurement Data */
        if (delta_t > 0.000001) {
            delta_q_lat = (LOS_ang_lat[i] - q_lat_last);
            delta_q_lon = (LOS_ang_lon[i] - q_lon_last);
            LOS_ang_lat_dot[i] = (delta_q_lat) / delta_t;
            LOS_ang_lon_dot[i] = (delta_q_lon) / delta_t;

            delta_R = LOS_dist[i] - R_last;
            delta_R_lat = LOS_dist_lat[i] - R_lat_last;
            LOS_dist_dot[i] = delta_R / delta_t;
            LOS_dist_lat_dot[i] = delta_R_lat / delta_t;
        }
    }
    Last_Target_State = Target_State;    

    i_min = Find_Min(delta_aZEM);
    /* Extract variables required to calculate metrics */
    dist_LOS = LOS_dist[i_min];
    speed = speed_M;
    altitude = State[i_min].z_or_alt;
    dist_covered = energy[i_min];
    time_taken = State[i_min].time;
    accuracy = 1;
    num = num_M;
};

bool Coordinated_TPNG::Interceptor::Get_State(objectstate& xo_state) {
    bool valid_state = Is_Started() && !my_killed;
    objectstate State_Out;
    int i_min = Find_Min(delta_aZEM);

    State_Out = State[i_min];
    State_Out.y_or_lon = deg_to_long(State[i_min].y_or_lon);
    State_Out.x_or_lat = deg_to_lat(State[i_min].x_or_lat);
    if (valid_state)
    {
        xo_state = State_Out;
    }
    return valid_state;
}
void Coordinated_TPNG::Interceptor::Get_Record() {
    int i;
    if (save_csv && valid_csv)
    {
        for (i = 0; i < num_M; i++)
        {
            Rec[i].csv("Output/TPNG_States_"+std::to_string(i), {"Time", "Latitude", "Longitude", "Altitude", "Heading_Angle", "Path_Angle", "LOS_Angle_Lat", "LOS_Angle_Lon", "Command_Lat", "Command_Lon"});
        }  
    }
}

bool Coordinated_TPNG::Interceptor::Is_Started() {
    int i_min = Find_Min(delta_aZEM);
    return (phase[i_min] > 0);
}
bool Coordinated_TPNG::Interceptor::Is_Hit(int index) {
    return (LOS_dist[index] < dist_stop);
}
bool Coordinated_TPNG::Interceptor::Is_Missed(int last_phase, int index) {
    //return (phase[index] < last_phase);
    return 0;
}

bool Coordinated_TPNG::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
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
