#pragma once
#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include "i_Interceptor.hpp"
#include "Auxiliary.h"

using namespace asc;                    /* Namespace of Ascent library */
using namespace std;

namespace NonLinear_Proportional
{
    static double command;               /* Interceptor Normal Acceleration Command */
    static double head_ang;              /* Interceptor Heading Angle */               
    static double radial_distance;       /* Interceptor Radial Distance from Centre of the Earth */

    class Interceptor :public Interceptor_Template {
    private:

        // Interceptor States
        objectstate Initial_State{ 0 };
        double Initial_Velocity{ 0 };
        objectstate State{ 0 };		    /* Interceptor Position */
        double R{ 0 };				    /* Relative Distance */
        double q{ 0 };
        Recorder Rec1;                  /* To record the states of the Interceptor */
        double V_M;                     /* Interceptor Velocity */
        int phase{ 0 };
        double x_dist{ 0 }, y_dist{ 0 };         /* Distance covered along x and y-axis */
        double energy{ 0 }, energy_lim{ 0 };     /* Energy consumed by the interceptor and its limit */

        // Target States
        objectstate Last_Target_State{ 0 };

        // Measurement Data
        double q_ddot{ 0 };
        double q_dot_est{ 0 };
        double V_R{ 0 };
        double t_last{ 0 };
        double q_last{ 0 };
        double R_last{ 0 };
        double q_dot{ 0 };

        // Parameters of NLPG
        double K_P{ 0 }, K_D{ 0 };              /* Proportional and Differential Gains for the Normal Acceleration Command */
        double Re{ 0 };
        double K{ 0 };                          /* Gain of Non-Linear Tracking Differentiator */

        /* Parameters to start and stop execution when Target hit */
        double R_stop{ 0 };
        double R_start{ 0 };
        bool save_csv = false;

        /* Private Methods */
        void NL_Differentiator(state_t&, double);        /* Non-linear Tracking Differentiator */
        void NLPG();                    /* To calculate the Non-Linear Proportional Guidance Command */
        struct Course {				    /* Solves the dynamic state equations of Hypersonic Gliding Vehicle using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
        void Track(double);             /* Calculate the commands at time t_meas and track the Target */
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

      public:
        Interceptor() {};
        ~Interceptor() {};
        bool Initialize(const string&, string&);       /* Initializes all the necessary parameters via the xml file */
        void Reinitialize();
        void Update_Target_State(objectstate);         /* This function calculates all the required quantities for tracking from position measurements of Target and Interceptor */
        bool Get_State(objectstate&);
        void Get_Record();
        bool Is_Started();
        bool Is_Terminated(bool&, bool&, int&);
    };
}

