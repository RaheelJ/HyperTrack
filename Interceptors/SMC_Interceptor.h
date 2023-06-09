#pragma once
#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include "i_Interceptor.hpp"
#include "Auxiliary.h"

using namespace asc;                    /* Namespace of Ascent library */
using namespace std;

namespace SMC
{
    const double Re = 6378100;
    static double command_lat;          /* Normal Overload Command in Lateral Plane */
    static double command_lon;          /* Normal Overload Command in Longitudnal Plane */
    static double head_ang;             /* Heading Angle */
    static double path_ang;             /* Flight Path Angle */

    class Interceptor :public Interceptor_Template {
    private:

        // Interceptor States
        objectstate Initial_State{ 0 };
        double Initial_Velocity{ 0 }, Initial_Altitude{ 0 };
        objectstate State{ 0 };		                /* Interceptor Position */
        double R{ 0 };				                /* LOS Distance */
        double q_lat{ 0 }, q_lon{ 0 };              /* LOS Angle in Lateral and Longitudnal Plane*/
        double lead_lat_m{ 0 }, lead_lon_m{ 0 };    /* Leading Angles in Lateral and Longitudnal Planes */
        Recorder Rec1;                              /* To record the states of the Interceptor */
        double Vm{ 0 };                             /* Interceptor Velocity */
        int phase{ 0 }; 
        double x_dist{ 0 }, y_dist{ 0 };            /* Distance covered along x, y and z-axis */
        double energy{ 0 }, energy_lim{ 0 };        /* Energy consumed by the interceptor and its limit */
        
        // Target States
        objectstate Last_Target_State{ 0 };
        double lead_lat_t{ 0 }, lead_lon_t{ 0 };    /* Leading Angles of Target in Lateral and Longitudnal Planes */
        
        // Measurement Data
        double q_lat_dot{ 0 }, q_lon_dot{ 0 };

        // Parameters of SMC Guidance Law
        double alpha{ 0 }, p{ 0 };
        double r0{ 0 }, r1{ 0 }, r2{ 0 };
        double n1{ 0 }, n2{ 0 };
        double k1{ 0 }, k2{ 0 };

        /* Parameters to start and stop execution when Target hit */
        double R_stop{ 0 };
        double R_start{ 0 };
        bool save_csv = false;
       
        // Flags to Control the Execution
        bool my_terminated = false;
        bool my_killed = false;
        bool self_destruct = false;
        bool valid_csv = false;
        bool is_first_run = true;
        
        /* Private Methods */
        void SMCGuidance();
        struct Course {				    /* Solves the dynamic state equations of Hypersonic Gliding Vehicle using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
        void Track(double);                     /* Calculate the commands at time t_meas and track the Target */        
        void Check_Param();
        void First_Run(objectstate);                            /* During the first execution calculate the required quantities */
        void Fine_Calculations(objectstate, objectstate);       /* Performs fine step calculations to detect the collision with greater accuracy */
        bool Is_Hit();
        bool Is_Missed(int);
        void Reset();

    public:
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

