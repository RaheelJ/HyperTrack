#pragma once
#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include "i_Interceptor.hpp"
#include "Auxiliary.h"
#include <Eigen/Dense>

using namespace asc;                    /* Namespace of Ascent library */
using namespace std;
using namespace Eigen;

namespace MPC
{
    // Fixed Parameters
    const double J_x = 1.615;			        // kg.m^2
    const double J_y = 247.4367;				// kg.m^2
    const double J_z = 247.4367;				// kg.m^2
    const double L_ = 0.2286;					// m
    const double Lm = 2;						// m
    const double m_z_alpha = -15.2693;
    const double m_z_w = -0.1614;
    const double m_z_delta = -11.8023;
    const double m_y_beta = -15.2693;
    const double m_y_w = -0.1614;
    const double m_y_delta = -11.8023;
    const double m_x_delta = -7.2766;
    const double K_M_y = 0.3;
    const double K_M_z = 0.3;
    const double Re = 6371e3;                   //Radius of the Earth (m)
    const double g = 9.8;
    const double atm_constant = 1.225;		    //Atmospheric Constant (Kg/m3)
    const double Hs = 10.4e3;				    //Scaling Height
    const double S = 0.04087734;			    //Reference Area (m2)
    const double C_z = -19.7995;
    const double C_y = 19.7995;
    const double m = 204.0227;			        // Mass (kg)

    //Near Space Interceptor (NSI) Interface Parameters
    static MatrixXd Bc(3, 5);
    static double Vm{ 1400 }, q{ 0 };           // Speed of the Interceptor, Dynamic Pressure
    static double Fsmax{ 0 };                   // Maximum Thrust by the Reaction Jets (N)
    static double n_y{ 0 }, n_z{ 0 };           // Normal Overload (Lateral and Longitudnal)
    static VectorXd d{ {0, 0, 0, 0, 0, 0} };    // Disturbance Vector
    static VectorXd u{ {0, 0, 0, 0, 0} };       // Inputs              
    /*  u(0) = Aileron Deflection
        u(1) = Rudder Deflection
        u(2) = Elevator Deflection
        u(3) = Eq. Control Surface Deflection in Pitch Axis
        u(4) = Eq. Control Surface Deflection in Yaw Axis
    */

    //Finite Time Disturbance Observer (FTDO) Interface Parameters
    static double L[6] = { 20, 20, 20, 50, 50, 50 };
    static double lambda_0{ 2 }, lambda_1{ 1.5 }, lambda_2{ 1.1 };
    static VectorXd x_{ {0,0,0,0,0,0,0} };
    static Vector3d v{ {0, 0, 0} };
    static VectorXd v_0{ {0,0,0,0,0,0} }, v_1{ {0,0,0,0,0,0} }, v_2{ {0,0,0,0,0,0} };

    class Interceptor :public Interceptor_Template {
    public:

        // Near Space Interceptor
        objectstate Initial_State, State;           //Initial Position, Current Position 
        Vector3d h_x{ {0, 0, 0} }, y{ {0, 0, 0} };  //Outputs              
        double bank_ang = 0;			            //Bank Angle
        double slip_ang = 0;			            //Sideslip Angle
        double attack_ang = 0;			            //Angle of Attack
        double pitch_ang{ 0 };                      //Pitch Angle
        double path_ang{0}, head_ang{0};            //Flight Path and Heading Angle
        double rate_roll{ 0 }, rate_yaw{ 0 }, rate_pitch{ 0 };
        double Ts{ 1e-6 };                          //Sampling time
        struct NSI_Model {				            /* Solves the dynamic state equations of Near Space Interceptor */
            void operator()(const state_t&, state_t&, const double);
        };
        void NSI_Run();                             /* Runs the NSI model given initial conditions and the inputs */

        // Finite Time Disturbance Observer
        double z_0[6]{ 0,0,0,0,0,0 }, z_1[6]{ 0,0,0,0,0,0 }, z_2[6]{ 0,0,0,0,0,0 };
        double Obs_x[6]{ 0,0,0,0,0,0 };
        VectorXd Obs_d{ {0,0,0,0,0,0} };
        double T{ 0 };                              //Current time
        struct FTDO_Model {				            /* Solves the dynamic state equations of FTDO */
            void operator()(const state_t&, state_t&, const double);
        };
        void FTDO_Run();                            /* Runs the FTDO model given initial conditions and the inputs */

        // Non-Linear MPC
        double k0_1{ 465.1317 }, k0_2{ 465.1317 }, k0_3{ 465.1317 };        // Controller gains calculated using                     
        double k1_1{ 34.9424 }, k1_2{ 34.9424 }, k1_3{ 34.9424 };           // predictive period (Tp) and the control order (r)
        Vector3d Ref;                               /* Reference input provided by the guidance law */
        Vector3d Ref_Last, dRef_Last, ddRef_Last;   /* Reference input and its derivative at previous time instant */
        double t_last = 0;                          /* Previous time instant */
        void NMP_Control(double);

        //Dynamic Control Allocation
        double lambda_c = 1;
        VectorXd u_max{ {30, 30, 30, 1, 1} };                       //Upper limit of NSI model inputs (deg)
        VectorXd u_min{ {-30, -30, -30, -1, -1} };                  //Lower limit of NSI model inputs (deg)
        VectorXd u_rate_max{ {300, 300, 300, 200, 200} };           //Upper rate limit of NSI model inputs (deg/s)
        VectorXd u_rate_min{ {-300, -300, -300, -200, -200} };      //Lower rate limit of NSI model inputs (deg/s)
        VectorXd u_T{ { 0,0,0,0,0 } }, u_2T{ { 0,0,0,0,0 } };       //Value of control inputs on previous instants
        VectorXd u_cs{ { 0,0,0,0,0 } };                             //Steady state input values
        void DCA();

        void Track(double);                                         /* Calculate the commands at time t_meas and track the Target */
        bool Is_Hit();
        bool Is_Missed(int);
        bool my_terminated = false;
        bool my_killed = false;
        bool self_destruct = false;
        bool valid_csv = false;
        bool is_first_run = true;
        void Check_Param();
        void First_Run(objectstate);                            /* During the first execution calculate the required quantities */
        void Fine_Calculations(objectstate, objectstate);       /* Performs fine step calculations to detect the collision with greater accuracy */
        void Reset();

        Interceptor() {};
        bool Initialize(const string&, string&);       /* Initializes all the necessary parameters via a "Parameter.xml" file */
        void Reinitialize();
        void Update_Target_State(objectstate);/* This function calculates all the required quantities for tracking from position measurements of Target and Interceptor */
        bool Get_State(objectstate&);
        void Get_Record();
        bool Is_Started();
        bool Is_Terminated(bool&, bool&, int&);
    };
}

