#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Auxiliary.h"

using namespace GeographicLib;

void Read_CSV(double t_meas, string file_name, objectstate& Target_State)
{
    fstream fin;                        // File pointer
    int time = 0, count = 0;

    fin.open(file_name);                // Open an existing file

    vector<string> row;                 // Read the data from the file as String Vector
    string line, word;
    getline(fin, line, '\n');

    while (fin) {
        row.clear();
        getline(fin, line, '\n');       // read an entire row and store it in a string variable 'line'
        stringstream s(line);           // used for breaking words

        // read every column data of a row and store it in a string variable, 'word'
        while (getline(s, word, ',')) {
            row.push_back(word);        // add all the column data of a row to a vector
        }

        time = stoi(row[0]);            // convert string to integer for comparision
        if (time == (int)t_meas) {      // Compare the time
            count = 1;
            Target_State.x_or_lat = stod(row[2]);
            Target_State.y_or_lon = stod(row[3]);
            Target_State.z_or_alt = stod(row[4]);
            Target_State.x_vel = stod(row[5]);
            Target_State.y_vel = stod(row[6]);
            Target_State.z_vel = stod(row[7]);
            Target_State.time = t_meas;
            break;
        }
    }
    if (count == 0)
        cout << "Record not found\n";
}
void Read_CSV(double t_meas, string file_name, TargetState& Target_State)
{
    fstream fin;                        // File pointer
    int time = 0, count = 0;

    fin.open(file_name);                // Open an existing file

    vector<string> row;                 // Read the data from the file as String Vector
    string line, word;
    getline(fin, line, '\n');

    while (fin) {
        row.clear();
        getline(fin, line, '\n');       // read an entire row and store it in a string variable 'line'
        stringstream s(line);           // used for breaking words

        // read every column data of a row and store it in a string variable, 'word'
        while (getline(s, word, ',')) {
            row.push_back(word);        // add all the column data of a row to a vector
        }

        time = stoi(row[0]);            // convert string to integer for comparision
        if (time == (int)t_meas) {      // Compare the time
            count = 1;
            Target_State.x_or_lat = stod(row[2]);
            Target_State.y_or_lon = stod(row[3]);
            Target_State.z_or_alt = stod(row[4]);
            Target_State.x_vel = stod(row[5]);
            Target_State.y_vel = stod(row[6]);
            Target_State.z_vel = stod(row[7]);
            Target_State.time = t_meas;
            break;
        }
    }
    if (count == 0)
        cout << "Record not found\n";
}

void Calc_Local(objectstate Geo_Position, array<double, 2>& Local_Position) {    /* Calculate Local Cartesian Coordinates from Geographic Coordinates */
    double temp_alt;
    Geocentric Earth(Constants::WGS84_a(), Constants::WGS84_f());
    LocalCartesian Local(0, 0, 0, Earth);                                       /* Reference Position for Local Cartesian Coordinates */
    Local.Forward(Geo_Position.x_or_lat, Geo_Position.y_or_lon, 0, Local_Position[0], Local_Position[1], temp_alt);
}
void Calc_Geo(objectstate& Geo_Position, array<double, 2> Local_Position, objectstate Reference) {
    double temp_alt;
    Geocentric Earth(Constants::WGS84_a(), Constants::WGS84_f());
    LocalCartesian Local(Reference.x_or_lat, Reference.y_or_lon, 0, Earth);     /* Reference Position for Local Cartesian Coordinates */
    Local.Reverse(Local_Position[0], Local_Position[1], 0, Geo_Position.x_or_lat, Geo_Position.y_or_lon, temp_alt);
}
bool Calc_Distance(objectstate Interceptor_State, objectstate Target_State, double& xo_dist, double& xo_angle) {
	Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
    double angle2;
	geod.Inverse(Interceptor_State.x_or_lat, Interceptor_State.y_or_lon, Target_State.x_or_lat, Target_State.y_or_lon, xo_dist, xo_angle, angle2);
	return true;
}
bool Calc_Distance(TargetState Interceptor_State, TargetState Target_State, double& xo_dist, double& xo_angle)
{
    Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
    double angle2;
    geod.Inverse(Interceptor_State.x_or_lat, Interceptor_State.y_or_lon, Target_State.x_or_lat, Target_State.y_or_lon, xo_dist, xo_angle, angle2);
    return true;
}
double Calc_Path_Ang(objectstate State)
{   //Angle is in radians
    double Vxy = sqrt(pow(State.x_vel, 2) + pow(State.y_vel, 2));
    double path_ang = atan(State.z_vel / Vxy);
    return path_ang;
}
double Calc_Path_Ang(TargetState State)
{   //Angle is in radians
    double Vxy = sqrt(pow(State.x_vel, 2) + pow(State.y_vel, 2));
    double path_ang = atan(State.z_vel / Vxy);
    return path_ang;
}
double Calc_Head_Ang(objectstate State)
{   //Angle is in radians
    double head_ang = atan2(State.x_vel, State.y_vel);
    return head_ang;
}
double Calc_Head_Ang(TargetState State)
{   //Angle is in radians
    double head_ang = atan2(State.x_vel, State.y_vel);
    return head_ang;
}
bool Calc_LOS(objectstate State1, objectstate State2, double& dist_LOS, double& angle_lat_LOS, double& angle_lon_LOS)
{   // Calculates the angles in radians

    double dist_lat_LOS;
    Calc_Distance(State1, State2, dist_lat_LOS, angle_lat_LOS);
    if (angle_lat_LOS < 0) { angle_lat_LOS = angle_lat_LOS + 360; }
    angle_lat_LOS = deg_to_rad(angle_lat_LOS);

    double delta_h = State2.z_or_alt - State1.z_or_alt;
    angle_lon_LOS = atan(delta_h / dist_lat_LOS);

    dist_LOS = sqrt(pow(dist_lat_LOS, 2) + pow(delta_h, 2));
    return true;
}
bool Calc_LOS(TargetState State1, TargetState State2, double& dist_LOS, double& angle_lat_LOS, double& angle_lon_LOS)
{
    double dist_lat_LOS;
    Calc_Distance(State1, State2, dist_lat_LOS, angle_lat_LOS);
    if (angle_lat_LOS < 0) { angle_lat_LOS = angle_lat_LOS + 360; }
    angle_lat_LOS = deg_to_rad(angle_lat_LOS);

    double delta_h = State2.z_or_alt - State1.z_or_alt;
    angle_lon_LOS = atan(delta_h / dist_lat_LOS);

    dist_LOS = sqrt(pow(dist_lat_LOS, 2) + pow(delta_h, 2));
    return true;
}
double Calc_Speed(objectstate State)
{
    double speed;
    speed = sqrt(pow(State.x_vel, 2) + pow(State.y_vel, 2) + pow(State.z_vel, 2));

    return speed;
}

void Limit(double& in, double upper_lim, double lower_lim)
{
    if (in <= lower_lim) { in = lower_lim; }
    else if (in >= upper_lim) { in = upper_lim; }
}
double signum(double in) {
    double out = ((in > 0) - (in < 0));
    return out;
}

double sind(double in) {                /* Sine Function with Argument in Degrees */
    return sin(fmod((in), 360) * M_PI / 180);
}
double cosd(double in) {                /* Cosine Function with Argument in Degrees */
    return cos(fmod((in), 360) * M_PI / 180);
}
double tand(double in) {                /* Cosine Function with Argument in Degrees */
    return tan(fmod((in), 360) * M_PI / 180);
}

double rad_to_deg(double in) {
    return 180 * fmod((in), 2 * M_PI) / M_PI;
}
double deg_to_rad(double in) {
    return M_PI * fmod((in), 360) / 180;
}

double Read_Param(double param, double lower_lim, double upper_lim, int& error)
{
    double out = 0;
    if (param <= upper_lim && param >= lower_lim) {
        out = param;
    }
    else {
        error = error + 1;
    }
    return out;
}
double deg_to_long(double in) {
    double out = 0;
    if (in <= -180) {
        out = in + 360;
    }
    else if (in >= 180) {
        out = out - 360;
    }
    else {
        out = in;
    }
    return out;
}
double deg_to_lat(double in) {
    double out;
    if (-270 <= in && in <= -90) {
        out = -180 - in;
    }
    else if (in <= -270) {
        out = 360 + in;
    }
    else if (90 <= in && in <= 270) {
        out = 180 - in;
    }
    else if (in >= 270) {
        out = -360 + in;
    }
    else {
        out = in;
    }
    return out;
}

double Derivative(double y1, double y2, double t1, double t2)
{
    double out;
    if ((t2 - t1) > 0)
    {
        out = (y2 - y1) / (t2 - t1);
    }
    else
        out = 0;
    return out;
}
double Scale(double in, double upper_lim, double lower_lim)
{
    double out = 1;
    if (in <= upper_lim && in >= lower_lim)
    {
        out = 1;
    }
    else if (in > upper_lim)
    {
        out = upper_lim / in;
    }
    else
    {
        out = lower_lim / in;
    }
    return out;
}
double Max(double in1, double in2)
{
    double out=0;
    if (in1 > in2)
    {
        out = in1;
    }
    else
    {
        out = in2;
    }
    return out;
}
double Min(double in1, double in2)
{
    double out = 0;
    if (in1 < in2)
    {
        out = in1;
    }
    else
    {
        out = in2;
    }
    return out;
}
int Find_Min(std::vector<double> X)
{
    int i, min_index = 0;
    double min_value = X[0];
    for (i = 1; i < X.size(); i++)
    {
        min_index = (X[i] < min_value) ? i : min_index;
        min_value = X[min_index];
    }
    return min_index;
}
std::vector<int> findNonZeroMinInd(const Eigen::MatrixXd& input_M)
{
    std::vector<int> indeces = { 0, 0 }; // index for the position of the element 
    double non_zero_min;
    bool get_min = false;
    for (int i = 0; i < input_M.rows(); ++i)
    {
        for (int j = 0; j < input_M.cols(); ++j)
        {
            if (input_M(i, j) != 0)
            {
                non_zero_min = input_M(i, j);
                indeces[0] = i;
                indeces[1] = j;
                get_min = true;
                break;
            }
        }
        if (get_min == true)
            break;
    }
    for (int i = 0; i < input_M.rows(); ++i)
    {
        for (int j = 0; j < input_M.cols(); ++j)
        {
            if (input_M(i, j) != 0 && input_M(i, j) < non_zero_min)
            {
                non_zero_min = input_M(i, j);
                indeces[0] = i;
                indeces[1] = j;
            }
        }
    }
    return indeces;
}

Point_3D CrossProduct(Point_3D V1, Point_3D V2)
{
    Point_3D Vout;
    Vout.x = V1.y * V2.z - V1.z * V2.y;
    Vout.y = V1.z * V2.x - V1.x * V2.z;
    Vout.z = V1.x * V2.y - V1.y * V2.x;

    return Vout;
}
bool If_Found(std::vector<int> collection, int key)
{
    int i;
    bool found = false;
    for (i = 0; i < collection.size(); i++)
    {
        if (collection[i] >= key)
        {
            found = true;
            break;
        }
    }
    return found;
}