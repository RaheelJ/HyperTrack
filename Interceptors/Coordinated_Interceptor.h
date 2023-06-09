#pragma once
#define _USE_MATH_DEFINES                   /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include "i_Interceptor.hpp"
#include "Auxiliary.h"
#include <Eigen/Dense>

using namespace asc;                       /* Namespace of Ascent library */
using namespace std;
using namespace Eigen;

namespace Coordinated
{
    static double command_L;                /* Normal Overload of the Leader */               
    static double radial_distance;          /* Radial Distance of the Interceptor's plane from centre of the Earth */
    static double head_ang_F;               /* Heading Angle of a Follower */
    static double V_F;                      /* Speed of Follower */
    static double u_x, u_y;                 /* Command of Follower in x and y direction */

    class Interceptor :public Interceptor_Template {
    private:

        // Leader Variables 
        objectstate Initial_State_L{ 0 };       /* Initial State of the Leader */
        double Initial_Velocity_L{ 0 };         /* Initial Speed of the Interceptor */
        objectstate State_L{ 0 };		        /* Current State of the Leader */
        double R_L{ 0 };				        /* LOS Distance between Traget and the Leader */
        double head_ang_L;                      /* Heading Angle of the Leader */
        double q_L{ 0 };                        /* LOS Angle of the Leader */
        Recorder Rec1;                          /* To record the states of the Leader */
        double V_L{ 0 };                        /* Leader's Speed */
        int phase{ 0 };                         /* Phase of Interception */
        double P_L[2]{ 0,0 };                   /* Distance covered by the Leader along x and y-axis */
        double energy{ 0 }, energy_lim{ 0 };    /* Energy consumed by each interceptor and its limit */

        // Parameters of MPNG
        double N{ 0 };                          /* Proportional Gain */
        double gain_const{ 0 };
        double Re{ 0 };                         /* Radius of Earth */

        // Target States
        objectstate Last_Target_State{ 0 };
        double head_ang_T{ 0 };                 /* Heading Angle of Target */
        
        // Follower Variables
        int n{ 0 };                             /* Number of Followers */
        double K{ 0 };                          /* Scale Coefficient */
        vector <objectstate> State_F{ 0 };      /* Current State of Followers */
        MatrixXd P_F{ {0} };                    /* Position of Followers in Local Coordinates */
        MatrixXd P_F_Last{ {0} };
        MatrixXd Initial_P_F{ {0} };            /* Initial Position of Followers in Local Coordinates */
        MatrixXd u_F{ {0} };                    /* Consensus Input */
        MatrixXd W{ {0} };                      /* Node Connection Matrix */
        VectorXd alpha{ 0 };        
        VectorXd beta{ 0 };                     /* Weight Vectors */
        VectorXd gamma{ 0 };
        vector <double> R_F{ 0 };				/* LOS Distance between Traget and the Leader */
        Recorder Rec2;                          /* To record the states of the Leader */
            
        // Measurement Data
        double q_dot{ 0 };
        double lead_ang_L{ 0 };                 /* Interceptor Velocity Leading Angle */
        double lead_ang_T{ 0 };                 /* Target Velocity Leading Angle */
        double t_last{ 0 };
        double q_last{ 0 };
        double a_T{ 0 };                        /* Target Normal Overload */

        /* Parameters to start and stop execution when Target hit */
        double R_stop{ 0 };
        double R_start{ 0 };
        bool save_csv = false;
        bool my_terminated = false;
        bool my_killed = false;
        bool self_destruct = false;

        /* Private Methods */
        void MPNG();
        struct Course_L {				        /* Solves the dynamic state equations of Leader using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
        void CGS(int);                             /* Consensus Protocol */
        struct Course_F {				        /* Solves the dynamic state equations of Followers using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
        void Track_L(double);                     /* Calculate the commands for Leader at time t_meas and track the Target */
        void Track_F(double);                     /* Calculate the commands for Followers at time t_meas and track the Target */
        bool Is_Hit();
        bool Is_Missed(int);
        bool valid_csv = false;
        bool is_first_run = true;
        void Check_Param();
        void First_Run(objectstate);                              /* During the first execution calculate the required quantities */
        void Fine_Calculations_L(objectstate, objectstate);       /* Performs fine step calculations to detect the collision with greater accuracy */
        void Fine_Calculations_F(objectstate, objectstate, int);  /* Performs fine step calculations to detect the collision with greater accuracy */
        void Reset();
        void Calc_Target_Overload(objectstate, double);
        void Initialize_Followers();

    public:
        Interceptor() {};
        bool Initialize(const string&, string&);    /* Initializes all the necessary parameters via a "Parameter.xml" file */
        void Reinitialize();
        void Update_Target_State(objectstate);      /* This function calculates all the required quantities for tracking from position measurements of Target and Interceptor */
        bool Get_State(objectstate&);
        void Get_Record();
        bool Is_Started();
        bool Is_Terminated(bool&, bool&, int&);
    };
}

