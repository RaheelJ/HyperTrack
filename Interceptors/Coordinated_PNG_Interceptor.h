#pragma once
#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include "i_Interceptor.hpp"
#include "Auxiliary.h"

using namespace asc;                    /* Namespace of Ascent library */
using namespace std;

namespace Coordinated_PNG
{
    const double Re = 6378100;
    static double aM_lon{ 0 }, aM_lat{ 0 };

    class Interceptor :public Interceptor_Template {
    private:

        // Interceptor States
        vector<objectstate> Initial_State{ 0 };
        double Initial_Velocity{ 0 }, Initial_Altitude{ 0 };
        vector<objectstate> State{ 0 };		                /* Interceptor Position */
        vector<double> LOS_dist_lat{ 0 }, LOS_dist_lon{ 0 }, LOS_dist{ 0 };    /* LOS Distance */
        vector<double> LOS_ang_lat{ 0 }, LOS_ang_lon{ 0 };  /* LOS Angle in Lateral and Longitudnal Plane*/
        vector<double> lead_lat_m{ 0 }, lead_lon_m{ 0 };    /* Leading Angles in Lateral and Longitudnal Planes */
        vector<double> head_ang;                            /* Heading Angle */
        vector<double> path_ang;                            /* Flight Path Angle */
        vector<Recorder> Rec;                                      /* To record the states of the Interceptor */
        double speed_M{ 0 }, energy_lim{ 0 };               /* Interceptor Velocity and its Energy Limit */
        vector<int> phase{ 0 }; 
        vector<double> x_dist{ 0 }, y_dist{ 0 };            /* Distance covered along x, y and z-axis */
        vector<double> energy{ 0 };                         /* Energy consumed by the interceptor*/
        
        // Target States
        objectstate Last_Target_State{ 0 };
        vector<double> lead_lat_t{ 0 }, lead_lon_t{ 0 };    /* Leading Angles of Target in Lateral and Longitudnal Planes */
        
        // Measurement Data
        vector<double> LOS_ang_lat_dot{ 0 }, LOS_ang_lon_dot{ 0 };
        vector<double> LOS_dist_lat_dot{ 0 }, LOS_dist_dot{ 0 };
        vector<double> LOS_dist_lon_dot{ 0 };

        // Parameters of Cooperative Coverage Strategy
        double a_M_max{ 0 }, aM_max{ 0 };                   /* Coverage and Max Interceptor Acceleration */
        double aT_max{ 0 };                                 /* Max Target Acceleration */
        double cov_ratio{ 0.5 };
        vector<Point_2D> CoC{ 0 };                          /* Centre of each Small Circle */
        int num_M{ 1 };                                     /* Total Number of interceptors given by the Coverage Strategy */

        // Parameters of PNG Law with Cooperative Bias
        double N{ 0 }, p{ 0 };                                  /* Navigation Ratio and Weighting Factor for Dynamic Adjustment */
        vector<double> lead_S_lat_t{ 0 }, lead_S_lon_t{ 0 };    /* Target Leading Angles for Interception in Standard Trajectory */
        vector<double> aT_S_lon{ 0 }, aT_S_lat{ 0 };            /* Target Accelerations for Standard Interceptor Trajectory */
        vector<double> B_lat{ 0 }, B_lon{ 0 };                  /* Cooperative Biases in Lateral and Longitudnal Plane */
        vector<double> command_lat;                             /* Normal Overload Command in Lateral Plane */
        vector<double> command_lon;                             /* Normal Overload Command in Longitudnal Plane */

        /* Parameters to start and stop execution when Target hit */
        double dist_stop{ 0 };
        double dist_start{ 0 };
        bool save_csv = false;
       
        // Flags to Control the Execution
        bool my_terminated = false;
        bool my_killed = false;
        bool self_destruct = false;
        bool valid_csv = false;
        bool is_first_run = true;
        
        /* Private Methods */
        double sin_DOmF(int);
        void Calc_Bias(int);
        void Adjust_Bias(int);
        void CPNG(int);
        struct Course {				    /* Solves the dynamic state equations of Hypersonic Gliding Vehicle using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
        void Coop_Coverage(int);                                        /* Run Algorithm to Find the Number of Interceptors and their Coverage Origins */
        void Track(double, int);                                        /* Calculate the commands at time t_meas and track the Target */        
        void Check_Param(int);
        void First_Run(objectstate);                                    /* During the first execution calculate the required quantities */
        void Fine_Calculations(objectstate, objectstate, int);          /* Performs fine step calculations to detect the collision with greater accuracy */
        bool Is_Hit(int);
        bool Is_Missed(int, int);
        void Reset();

    public:
        Interceptor() {};
        bool Initialize(const string&, string&);        /* Initializes all the necessary parameters via a "Parameter.xml" file */
        void Reinitialize();
        void Update_Target_State(objectstate);          /* This function calculates all the required quantities for tracking from position measurements of Target and Interceptor */
        bool Get_State(objectstate&);
        void Get_Record();
        bool Is_Started();
        bool Is_Terminated(bool&, bool&, int&);
    };
}

