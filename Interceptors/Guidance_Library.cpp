#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "Guidance_Library.h"
#include <iostream>
#include "xml/pugixml.hpp"
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>

using namespace Eigen;

void GuidanceLib::Interceptor::Reset() {
    /* Variables to calculate metrics */
    dist_LOS = -1;
    speed = 0;
    altitude = 0;
    dist_covered = 0;
    time_taken = 0;
    accuracy = 1;
    num = 1;

    /* Static Variable */
    command_x = 0;
    command_y = 0;
    command_z = 0;

    // Interceptor States
    R = 0;
    q_lon = 0;
    q_lat = 0;
    lead_lat_m = 0;
    lead_lon_m = 0;
    phase = 0;
    x_dist = 0;
    y_dist = 0;
    dist = 0;
    head_ang = 0;
    path_ang = 0;
    miss_phase = 0;
    miss_timer = 0;
    control_energy = 0;

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
bool GuidanceLib::Interceptor::Initialize(const string& XML_Parameters, string& message) {
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
    pugi::xml_node Algo = Root.child("Guidance_Lib");
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


    /* Initialize parameters from "Execution" node */
    R_start = Read_Param(stod(Exec_Param.child_value("Start_Distance")), 0, 10000000, error);
    R_stop = Read_Param(stod(Exec_Param.child_value("Stop_Distance")), 0, 5000, error);
    dist_lim = Read_Param(stod(Exec_Param.child_value("Distance_Limit")), 0, 10000000, error);
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
void GuidanceLib::Interceptor::Reinitialize() {
    State = Initial_State;
    Reset();
}

void GuidanceLib::Interceptor::Decompose_R(double& Rx, double& Ry, double& Rz)
{
    Rx = sin(q_lat) * R * cos(q_lon);
    Ry = cos(q_lat) * R * cos(q_lon);
    Rz = R * sin(q_lon);
}
void GuidanceLib::Interceptor::Guidance_Law()
{


}
void GuidanceLib::Interceptor::MPNGuidance(objectstate Target_State) {
    MatrixXd target_pos_M(1, 3), target_vel_M(1, 3);
    MatrixXd seeker_vel_M(1, 3), seeker_pos_M(1, 3);
    MatrixXd accel_M = Eigen::MatrixXd::Zero(1, 3);

    target_vel_M << Target_State.x_vel, Target_State.y_vel, Target_State.z_vel;
    seeker_vel_M << State.x_vel, State.y_vel, State.z_vel;

    double Rx, Ry, Rz, LOS_cos;
    Decompose_R(Rx, Ry, Rz);
    MatrixXd R_r{ {Rx, Ry, Rz } };
    double accelerate_until_t = time_acc;
    double sample_t = Target_State.time - State.time;
    double seeker_max_acc = command_max;

    // Since the seeker starts from velocity == 0, we need time for the seeker to gain velocity
    if (seeker_vel_M.norm() < target_vel_M.norm() && R_r.norm() > accelerate_until_t * sample_t * seeker_vel_M.norm())
        accel_M = R_r / R_r.norm() * seeker_max_acc;
    else
    {
        Eigen::MatrixXd A = seeker_vel_M / seeker_vel_M.norm();
        Eigen::MatrixXd M2 = A * R_r.transpose() / R_r.norm();
        LOS_cos = M2(0, 0);

        if (LOS_cos < 0.9397) // if it is larger than 20 degree
        {
            Eigen::MatrixXd M3 = R_r * A.transpose() * A;

            Eigen::MatrixXd B = Eigen::MatrixXd::Ones(A.rows(), A.cols());
            if (R_r.isApprox(M3))
            {
                Eigen::MatrixXd v_M = target_vel_M * A.transpose() * A;
                if (target_vel_M.isApprox(v_M))
                {
                    std::vector<int> indeces = findNonZeroMinInd(A);
                    int row_ind = indeces[0];
                    int col_ind = indeces[1];

                    Eigen::MatrixXd Z = A * B.transpose();
                    double z = Z(0, 0);
                    B(row_ind, col_ind) = (A(row_ind, col_ind) * B(row_ind, col_ind)
                        - z) / A(row_ind, col_ind);
                    B = B / B.norm();
                }
                else
                    B = (target_vel_M - v_M) / (target_vel_M - v_M).norm();
            }
            else
                B = (R_r - R_r * A.transpose() * A) / (R_r - R_r * A.transpose() * A).norm();
            if (B.norm() == 0)
                accel_M = Eigen::MatrixXd::Zero(B.rows(), B.cols());
            else
                accel_M = B / B.norm() * seeker_max_acc;  // * 9.81 * seeker_max_acc;
        }
        else
        {
            Eigen::MatrixXd lambda_M = R_r / R_r.norm();
            Eigen::MatrixXd M4 = target_vel_M - seeker_vel_M;
            Eigen::MatrixXd Vc_M = M4 * (-1.0) * lambda_M.transpose();
            Eigen::MatrixXd lambda_dot_M = (M4 + Vc_M * lambda_M) / R_r.norm();
            accel_M = Vc_M * lambda_dot_M;
            accel_M = accel_M * n;
        }
    }

    command_x = accel_M(0, 0);
    command_y = accel_M(0, 1);
    command_z = accel_M(0, 2);

    Limit(command_x, command_max, -command_max);
    Limit(command_y, command_max, -command_max);
    Limit(command_z, command_max, -command_max);
}
void GuidanceLib::Interceptor::APNGuidance(objectstate Target_State)
{
    MatrixXd target_pos_M(1, 3), target_vel_M(1, 3);
    MatrixXd seeker_vel_M(1, 3), seeker_pos_M(1, 3);
    MatrixXd acc_M = Eigen::MatrixXd::Zero(1, 3);
    MatrixXd target_acc_M = Eigen::MatrixXd::Zero(1, 3);

    target_vel_M << Target_State.x_vel, Target_State.y_vel, Target_State.z_vel;
    seeker_vel_M << State.x_vel, State.y_vel, State.z_vel;

    double Rx, Ry, Rz, LOS_cos;
    Decompose_R(Rx, Ry, Rz);
    MatrixXd R_r{ {Rx, Ry, Rz } };
    double accelerate_until_t = time_acc;
    double sample_t = Target_State.time - State.time;
    double seeker_max_acc = command_max;

    target_acc_M(0) = (Target_State.x_vel - Last_Target_State.x_vel) / sample_t;
    target_acc_M(1) = (Target_State.y_vel - Last_Target_State.y_vel) / sample_t;
    target_acc_M(2) = (Target_State.z_vel - Last_Target_State.z_vel) / sample_t;

    // Gain enough velocity for the seeker 
    if (seeker_vel_M.norm() < target_vel_M.norm() && R_r.norm() > accelerate_until_t * sample_t * seeker_vel_M.norm())
        acc_M = R_r / R_r.norm() * seeker_max_acc;
    else
    {
        Eigen::MatrixXd A = seeker_vel_M / seeker_vel_M.norm();
        // compute Line of Sight angle 
        Eigen::MatrixXd M2 = A * R_r.transpose() / R_r.norm();
        LOS_cos = M2(0, 0);
        if (LOS_cos < 0.9397) // if it is smaller than 20 degree
        {
            Eigen::MatrixXd M3 = R_r * A.transpose() * A;

            Eigen::MatrixXd B = Eigen::MatrixXd::Ones(A.rows(), A.cols());
            if (R_r.isApprox(M3))
            {
                Eigen::MatrixXd v_M = target_vel_M * A.transpose() * A;
                if (target_vel_M.isApprox(v_M))
                {
                    std::vector<int> indeces = findNonZeroMinInd(A);
                    int row_ind = indeces[0];
                    int col_ind = indeces[1];

                    Eigen::MatrixXd Z = A * B.transpose();
                    double z = Z(0, 0);
                    B(row_ind, col_ind) = (A(row_ind, col_ind) * B(row_ind, col_ind)
                        - z) / A(row_ind, col_ind);
                    B = B / B.norm();
                }
                else
                    B = (target_vel_M - v_M) / (target_vel_M - v_M).norm();
            }
            else
                B = (R_r - R_r * A.transpose() * A) / (R_r - R_r * A.transpose() * A).norm();

            // get the acceleration matrix 
            Eigen::MatrixXd max_acc_M = Eigen::MatrixXd::Zero(1, 3);
            for (int i = 0; i < 3; ++i)
            {
                max_acc_M(0, i) = R_r(0, i) / std::abs(R_r(0, i)) * seeker_max_acc;
            }
            if (B.norm() == 0)
            {
                acc_M = max_acc_M;
            }
            else
            {
                acc_M = B * 9.81 * seeker_max_acc + max_acc_M;
            }

        }
        else
        {
            Eigen::MatrixXd lambda_M = R_r / R_r.norm();
            Eigen::MatrixXd V_r = target_vel_M - seeker_vel_M;
            Eigen::MatrixXd V_c = V_r * (-1.0) * lambda_M.transpose();
            Eigen::MatrixXd lambda_dot_M = (V_r + V_c * lambda_M) / R_r.norm();
            Eigen::MatrixXd acc_M1 = V_c * lambda_dot_M;
            acc_M1 = acc_M1 * n;

            //Eigen::MatrixXd acc_M;
            Eigen::MatrixXd alpha_T = target_vel_M * R_r.transpose() / target_vel_M.norm() / R_r.norm();
            double cos_alpha = alpha_T(0, 0);
            if (cos_alpha > 0)
                acc_M = acc_M1 + target_acc_M * (n / 2);
            else
                acc_M = acc_M1;
        }
    }

    command_x = acc_M(0, 0);
    command_y = acc_M(0, 1);
    command_z = acc_M(0, 2);

    Limit(command_x, command_max, -command_max);
    Limit(command_y, command_max, -command_max);
    Limit(command_z, command_max, -command_max);
}
void GuidanceLib::Interceptor::DGGuidance(objectstate Target_State)
{
    MatrixXd target_pos_M(1, 3), target_vel_M(1, 3);
    MatrixXd seeker_vel_M(1, 3), seeker_pos_M(1, 3);
    MatrixXd acc_M = Eigen::MatrixXd::Zero(1, 3);
    MatrixXd target_acc_M = Eigen::MatrixXd::Zero(1, 3);

    //target_pos_M << Target_State.x_or_lat, Target_State.y_or_lon, Target_State.z_or_alt;
    target_vel_M << Target_State.x_vel, Target_State.y_vel, Target_State.z_vel;
    //seeker_pos_M << State.x_or_lat, State.y_or_lon, State.z_or_alt;
    seeker_vel_M << State.x_vel, State.y_vel, State.z_vel;

    double Rx, Ry, Rz, LOS_cos;
    Decompose_R(Rx, Ry, Rz);
    MatrixXd R_r{ {Rx, Ry, Rz } };
    double accelerate_until_t = time_acc;
    double sample_t = Target_State.time - State.time;
    double seeker_max_acc = command_max;

    target_acc_M(0) = (Target_State.x_vel - Last_Target_State.x_vel) / sample_t;
    target_acc_M(1) = (Target_State.y_vel - Last_Target_State.y_vel) / sample_t;
    target_acc_M(2) = (Target_State.z_vel - Last_Target_State.z_vel) / sample_t;

    if (seeker_vel_M.norm() < target_vel_M.norm() && R_r.norm() > accelerate_until_t * sample_t * seeker_vel_M.norm())
        acc_M = R_r / R_r.norm() * seeker_max_acc;
    else
    {
        Eigen::MatrixXd A = seeker_vel_M / seeker_vel_M.norm();
        // compute Line of Sight angle 
        Eigen::MatrixXd M2 = A * R_r.transpose() / R_r.norm();
        LOS_cos = M2(0, 0);
        double distance = R_r.norm();

        if (LOS_cos < 0.9397) // if it is smaller than 20 degree
        {
            Eigen::MatrixXd M3 = R_r * A.transpose() * A;

            Eigen::MatrixXd B = Eigen::MatrixXd::Ones(A.rows(), A.cols());
            if (R_r.isApprox(M3))
            {
                Eigen::MatrixXd v_M = target_vel_M * A.transpose() * A;
                if (target_vel_M.isApprox(v_M))
                {
                    std::vector<int> indeces = findNonZeroMinInd(A);
                    int row_ind = indeces[0];
                    int col_ind = indeces[1];

                    Eigen::MatrixXd Z = A * B.transpose();
                    double z = Z(0, 0);
                    B(row_ind, col_ind) = (A(row_ind, col_ind) * B(row_ind, col_ind)
                        - z) / A(row_ind, col_ind);
                    B = B / B.norm();
                }
                else
                    B = (target_vel_M - v_M) / (target_vel_M - v_M).norm();
            }
            else
                B = (R_r - R_r * A.transpose() * A) / (R_r - R_r * A.transpose() * A).norm();
            if (B.norm() == 0)
            {
                acc_M = R_r * seeker_max_acc / R_r.norm();
            }
            else
            {
                acc_M = B * 9.81 * seeker_max_acc + R_r * seeker_max_acc / R_r.norm();
            }
        }
        else
        {
            // compute the angles between the V_r and R_r
            //Eigen::MatrixXd R_r = target_pos_M - seeker_pos_M;
            Eigen::MatrixXd eta_M = seeker_vel_M * R_r.transpose() / seeker_vel_M.norm() / R_r.norm();
            double cos_m = eta_M(0, 0);
            Eigen::MatrixXd eta_T = target_vel_M * R_r.transpose() / target_vel_M.norm() / R_r.norm();
            double cos_t = eta_T(0, 0);
            // the second part of the acceleration command 
            Eigen::MatrixXd acc_M2 = target_acc_M * (cos_t / cos_m);

            Eigen::MatrixXd lambda_M = R_r / R_r.norm();
            Eigen::MatrixXd V_r = target_vel_M - seeker_vel_M;
            Eigen::MatrixXd V_c = V_r * (-1.0) * lambda_M.transpose();
            Eigen::MatrixXd lambda_dot_M = (V_r + V_c * lambda_M) / R_r.norm();
            Eigen::MatrixXd acc_M1 = V_c * lambda_dot_M;
            acc_M1 = acc_M1 * (n / cos_m);
            acc_M = acc_M1 + acc_M2;
        }
    }

    command_x = acc_M(0, 0);
    command_y = acc_M(0, 1);
    command_z = acc_M(0, 2);

    Limit(command_x, command_max, -command_max);
    Limit(command_y, command_max, -command_max);
    Limit(command_z, command_max, -command_max);
}
void GuidanceLib::Interceptor::LOSGuidance(objectstate Target_State)
{
    MatrixXd target_pos_M(1, 3), target_vel_M(1, 3);
    MatrixXd seeker_vel_M(1, 3), seeker_pos_M(1, 3);
    MatrixXd acc_M = Eigen::MatrixXd::Zero(1, 3);
    MatrixXd target_acc_M = Eigen::MatrixXd::Zero(1, 3);

    double target_dist, target_azimuth, target_elevation;
    Calc_LOS(Initial_State, Target_State, target_dist, target_azimuth, target_elevation);
    target_pos_M(0) = sin(target_azimuth) * target_dist * cos(target_elevation);
    target_pos_M(1) = cos(target_azimuth) * target_dist * cos(target_elevation);
    target_pos_M(2) = target_dist * sin(target_elevation);
    seeker_pos_M << x_dist, y_dist, State.z_or_alt;

    target_vel_M << Target_State.x_vel, Target_State.y_vel, Target_State.z_vel;
    seeker_vel_M << State.x_vel, State.y_vel, State.z_vel;

    double accelerate_until_t = time_acc;
    double sample_t = Target_State.time - State.time;
    double seeker_max_acc = command_max;

    // convert the matrix to vector
    Eigen::Vector3d origin_pos(0, 0, 0);
    Eigen::Vector3d R_t(target_pos_M(0, 0), target_pos_M(0, 1), target_pos_M(0, 2));
    R_t = R_t - origin_pos;
    Eigen::Vector3d V_t(target_vel_M(0, 0), target_vel_M(0, 1), target_vel_M(0, 2));
    Eigen::Vector3d R_m(seeker_pos_M(0, 0), seeker_pos_M(0, 1), seeker_pos_M(0, 2));
    R_m = R_m - origin_pos;
    Eigen::Vector3d V_m(seeker_vel_M(0, 0), seeker_vel_M(0, 1), seeker_vel_M(0, 2));
    Eigen::Vector3d R_r = R_t - R_m;

    Eigen::Vector3d acc;
    if (V_m.norm() == 0 || R_m.norm() == 0)
        acc = R_t / R_t.norm() * seeker_max_acc;
    else
    {
        Eigen::Vector3d Vec0 = R_t.cross(V_t);
        Eigen::Vector3d V_m_new = Vec0.cross(R_m) / std::pow(R_t.norm(), 2) + R_t / R_t.norm() * V_t.norm();
        acc = (V_m_new - V_m) / sample_t;
    }

    command_x = acc_M(0, 0);
    command_y = acc_M(0, 1);
    command_z = acc_M(0, 2);

    Limit(command_x, command_max, -command_max);
    Limit(command_y, command_max, -command_max);
    Limit(command_z, command_max, -command_max);
}
void GuidanceLib::Interceptor::PursuitGuidance(objectstate Target_State) 
{
    MatrixXd target_pos_M(1, 3), target_vel_M(1, 3);
    MatrixXd seeker_vel_M(1, 3), seeker_pos_M(1, 3);
    MatrixXd acc_M = Eigen::MatrixXd::Zero(1, 3);
    MatrixXd target_acc_M = Eigen::MatrixXd::Zero(1, 3);

    double target_dist, target_azimuth, target_elevation;
    Calc_LOS(Initial_State, Target_State, target_dist, target_azimuth, target_elevation);
    target_pos_M(0) = sin(target_azimuth) * target_dist * cos(target_elevation);
    target_pos_M(1) = cos(target_azimuth) * target_dist * cos(target_elevation);
    target_pos_M(2) = target_dist * sin(target_elevation);
    seeker_pos_M << x_dist, y_dist, State.z_or_alt;

    target_vel_M << Target_State.x_vel, Target_State.y_vel, Target_State.z_vel;
    seeker_vel_M << State.x_vel, State.y_vel, State.z_vel;

    double accelerate_until_t = time_acc;
    double sample_t = Target_State.time - State.time;
    double seeker_max_acc = command_max;

    Eigen::Vector3d R_t(target_pos_M(0, 0), target_pos_M(0, 1), target_pos_M(0, 2));
    Eigen::MatrixXd R_r_M = target_pos_M - seeker_pos_M;
    Eigen::Vector3d R_r(R_r_M(0, 0), R_r_M(0, 1), R_r_M(0, 2));
    Eigen::Vector3d V_t(target_vel_M(0, 0), target_vel_M(0, 1), target_vel_M(0, 2));
    Eigen::Vector3d V_m(seeker_vel_M(0, 0), seeker_vel_M(0, 1), seeker_vel_M(0, 2));
    Eigen::Vector3d Vec_0 = V_m.cross(R_r) / std::pow(R_r.norm(), 2);
    Eigen::Vector3d acc = Vec_0.cross(V_m);
    
    if ((V_m.norm() < V_t.norm()) && (R_r.norm() > accelerate_until_t * sample_t * V_m.norm()))
        acc = R_r / R_r.norm() * seeker_max_acc;
    
    command_x = acc_M(0, 0);
    command_y = acc_M(0, 1);
    command_z = acc_M(0, 2);

    Limit(command_x, command_max, -command_max);
    Limit(command_y, command_max, -command_max);
    Limit(command_z, command_max, -command_max);
}
void GuidanceLib::Interceptor::PurePNGuidance(objectstate Target_State)
{
    MatrixXd target_pos_M(1, 3), target_vel_M(1, 3);
    MatrixXd seeker_vel_M(1, 3), seeker_pos_M(1, 3);
    MatrixXd acc_M = Eigen::MatrixXd::Zero(1, 3);
    MatrixXd target_acc_M = Eigen::MatrixXd::Zero(1, 3);

    double target_dist, target_azimuth, target_elevation;
    Calc_LOS(Initial_State, Target_State, target_dist, target_azimuth, target_elevation);
    target_pos_M(0) = sin(target_azimuth) * target_dist * cos(target_elevation);
    target_pos_M(1) = cos(target_azimuth) * target_dist * cos(target_elevation);
    target_pos_M(2) = target_dist * sin(target_elevation);
    seeker_pos_M << x_dist, y_dist, State.z_or_alt;

    target_vel_M << Target_State.x_vel, Target_State.y_vel, Target_State.z_vel;
    seeker_vel_M << State.x_vel, State.y_vel, State.z_vel;

    double accelerate_until_t = time_acc;
    double sample_t = Target_State.time - State.time;
    double seeker_max_acc = command_max;

    Eigen::MatrixXd V_r_M = target_vel_M - seeker_vel_M;
    Eigen::Vector3d V_r(V_r_M(0, 0), V_r_M(0, 1), V_r_M(0, 2));
    Eigen::MatrixXd R_r_M = target_pos_M - seeker_pos_M;
    Eigen::Vector3d R_r(R_r_M(0, 0), R_r_M(0, 1), R_r_M(0, 2));
    Eigen::Vector3d V_t(target_vel_M(0, 0), target_vel_M(0, 1), target_vel_M(0, 2));
    Eigen::Vector3d V_m(seeker_vel_M(0, 0), seeker_vel_M(0, 1), seeker_vel_M(0, 2));
    Eigen::Vector3d omega = R_r.cross(V_r) / (std::pow(R_r.norm(), 2));
    Eigen::Vector3d acc = V_r.cross(omega) * n;
    
    if (V_m.norm() == 0)
        acc = R_r / R_r.norm() * seeker_max_acc;

    command_x = acc_M(0, 0);
    command_y = acc_M(0, 1);
    command_z = acc_M(0, 2);

    Limit(command_x, command_max, -command_max);
    Limit(command_y, command_max, -command_max);
    Limit(command_z, command_max, -command_max);
}

void GuidanceLib::Interceptor::Course::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] -- Latitude (radians)
        x[1] -- Longitude (radians)
        x[2] -- x distance (m)
        x[3] -- y distance (m)
        x[4] -- Altitude (m)
        x[5] -- x-velocity (m/s)
        x[6] -- y-velocity (m/s)
        x[7] -- z-velocity (m/s) */

    dx[0] = x[6]/(Re+x[4]);
    dx[1] = x[5] / ((Re + x[6]) * cos(x[0]));
    dx[2] = x[5];
    dx[3] = x[6];
    dx[4] = x[7];
    dx[5] = command_x;
    dx[6] = command_y;
    dx[7] = command_z;
}
void GuidanceLib::Interceptor::Check_Param() {
    bool condition1 = (State.z_or_alt < 0) || (abs(State.x_or_lat) >= 89.999999);
    self_destruct = condition1;
}
void GuidanceLib::Interceptor::Track(objectstate Target_State) {
    state_t x = { 0,0,0,0,0,0,0,0 };
    x[0] = deg_to_rad(State.x_or_lat);           /* Initial state values for the current iteration */
    x[1] = deg_to_rad(State.y_or_lon);
    x[2] = x_dist;
    x[3] = y_dist;
    x[4] = State.z_or_alt;
    x[5] = State.x_vel;
    x[6] = State.y_vel;
    x[7] = State.z_vel;

    double t = State.time;
    double t_meas = Target_State.time;
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
        else if (abs(State.time - t) < 0.000001) {
            MPNGuidance(Target_State);
        }
        if (save_csv)
        {
            Rec1({ t, rad_to_deg(x[0]), rad_to_deg(x[1]), x[4], rad_to_deg(q_lat), rad_to_deg(q_lon), command_x, command_y, command_z });
            valid_csv = true;
        }
        Integrator(Course(), x, t, step_size);        /* Solve the ODEs to find system states */   
        dist = sqrt(pow(x[2], 2) + pow(x[3], 2));
    }

    State.id = 0;
    State.time = t;
    State.z_or_alt = x[4];
    State.y_or_lon = rad_to_deg(x[1]);
    State.x_or_lat = rad_to_deg(x[0]);
    if (phase > 0)
    {
        State.x_vel = (std::abs(x[5]) > speed_max) ? ((x[5] / std::abs(x[5])) * speed_max) : x[5];
        State.y_vel = (std::abs(x[6]) > speed_max) ? ((x[6] / std::abs(x[6])) * speed_max) : x[6];
        State.z_vel = (std::abs(x[7]) > speed_max) ? ((x[7] / std::abs(x[7])) * speed_max) : x[7];
    }
    else 
    {
        State.x_vel = 0;
        State.y_vel = 0;
        State.z_vel = 0;
    }
    x_dist = x[2];
    y_dist = x[3];

    double last_Vm = Vm;
    double last_head_ang = head_ang;
    double last_path_ang = path_ang;
    Vm = sqrt(pow(State.x_vel, 2) + pow(State.y_vel, 2) + pow(State.z_vel, 2));
    head_ang = Calc_Head_Ang(State);
    path_ang = Calc_Path_Ang(State);

    double overload_lat = (head_ang - last_head_ang) / (last_Vm * cos(last_path_ang));
    control_energy = control_energy + pow(overload_lat, 2) * step_size;
}

void GuidanceLib::Interceptor::First_Run(objectstate Target_State) {
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
void GuidanceLib::Interceptor::Fine_Calculations(objectstate Last_Interceptor_State, objectstate Target_State) {
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
            Rec1({ i, Fine_Interceptor.x_or_lat, Fine_Interceptor.y_or_lon, Fine_Interceptor.z_or_alt, rad_to_deg(q_lat), rad_to_deg(q_lon), command_x, command_y, command_z });
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
                dist_covered = dist;
                time_taken = State.time;
                energy_i = control_energy;
            }
            speed = Vm;
            altitude = State.z_or_alt;
            accuracy = 1;
            num = 1;

            break;
        }
    }
}

void GuidanceLib::Interceptor::Update_Target_State(objectstate Target_State)
{
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
    
    Last_Interceptor_State = State;
    Track(Target_State);

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
    my_terminated = Is_Hit() || Is_Missed(R_last) || self_destruct || (dist > dist_lim);
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
        dist_covered = dist;
        time_taken = State.time;
        energy_i = control_energy;
    }
    speed = Vm;
    altitude = State.z_or_alt;
    accuracy = 1;
    num = 1;
};

bool GuidanceLib::Interceptor::Get_State(objectstate& xo_state) {
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
void GuidanceLib::Interceptor::Get_Record() {
    if (save_csv && valid_csv)
    {
        Rec1.csv("Output/SMC_States", { "Time", "Latitude", "Longitude", "Altitude", "LOS_Angle_Lat", "LOS_Angle_Lon", "Command_X", "Command_Y", "Command_Z"});
    }

}

bool GuidanceLib::Interceptor::Is_Started() {
    return (phase > 0);
}
bool GuidanceLib::Interceptor::Is_Hit() {
    return (R < R_stop);
}
bool GuidanceLib::Interceptor::Is_Missed(double R_last) {
    bool missed = 0;
    miss_phase = (R < R_last) ? R_last : miss_phase;
    if (R < miss_phase)
    {
        miss_timer++;
    }
    else
    {
        miss_timer = 0;
        miss_phase = 0;
    }
    missed = (miss_timer > 10) ? 1 : 0;
    return missed;
}

bool GuidanceLib::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
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
