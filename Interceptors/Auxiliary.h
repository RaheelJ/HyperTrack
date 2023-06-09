#pragma once

#include <Interceptor/i_Interceptor.hpp>
#include <Interceptor/i_target.hpp>
#include <string>
#include <array>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>

using namespace std;

struct Point_2D
{
    double lon;
    double lat;
};
struct Point_2D_Polar
{
    double mag;
    double ang;
};
struct Point_3D
{
    double x;
    double y;
    double z;
};
Point_3D CrossProduct(Point_3D, Point_3D);

void Read_CSV(double, string, objectstate&);
void Read_CSV(double, string, TargetState&);
void Calc_Local(objectstate, array<double, 2>&);
void Calc_Geo(objectstate&, array<double, 2>, objectstate);
bool Calc_Distance(objectstate Interceptor_State, objectstate Target_State, double& xo_dist, double& xo_angle);
bool Calc_Distance(TargetState Interceptor_State, TargetState Target_State, double& xo_dist, double& xo_angle);
double Calc_Path_Ang(objectstate State);
double Calc_Path_Ang(TargetState State);
double Calc_Head_Ang(objectstate State);
double Calc_Head_Ang(TargetState State);
bool Calc_LOS(objectstate State1, objectstate State2, double& dist_LOS, double& angle_lat_LOS, double& angle_lon_LOS);
bool Calc_LOS(TargetState State1, TargetState State2, double& dist_LOS, double& angle_lat_LOS, double& angle_lon_LOS);
double Calc_Speed(objectstate);

void Limit(double&, double, double);
double signum(double);

double sind(double);
double cosd(double);
double tand(double);

double rad_to_deg(double);
double deg_to_rad(double);

double Read_Param(double, double, double, int&);
double deg_to_long(double);
double deg_to_lat(double);

double Derivative(double, double, double, double);
double Scale(double, double, double);
double Max(double, double);
double Min(double, double);
int Find_Min(std::vector<double>);
bool If_Found(std::vector<int>, int);
std::vector<int> findNonZeroMinInd(const Eigen::MatrixXd&);

class Interceptor_Template {
public:
    double dist_LOS;
    double speed;
    double altitude;
    double dist_covered;
    double time_taken;
    double accuracy;
    double num;
    double energy_i;

    virtual ~Interceptor_Template() {};
    virtual bool Initialize(const std::string&, std::string&) = 0;
    virtual void Reinitialize() = 0;
    virtual void Update_Target_State(objectstate) = 0;
    virtual bool Get_State(objectstate&) = 0;
    virtual void Get_Record() = 0;
    virtual bool Is_Started() = 0;
    virtual bool Is_Terminated(bool&, bool&, int&) = 0;
};

class Target_Template {
public:
    std::vector<double> dist_destination;
    std::vector<double> time_destination;
    
    virtual ~Target_Template() {};
    virtual bool 	initialize(const std::string&, std::string&) = 0;
    virtual bool 	reinitialize() = 0;
    virtual bool 	get_state(double, std::vector<TargetState>&) = 0;
    virtual bool    get_target_property(int xi_target_id, TargetProperty&) = 0;
    virtual bool 	update_target(std::vector<TargetChange>&) = 0;
    //virtual void 	get_performance_measures(TargetPerfMeasure&) = 0;
    virtual int     get_expected_number_of_targets() = 0;
    virtual double  get_start_time() = 0;
    virtual bool    has_more_data() = 0;
};

