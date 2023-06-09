#pragma once
#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include "i_target.hpp"
#include "Auxiliary.h"

using namespace asc;                    /* Namespace of Ascent library */
using namespace std;

namespace TargetLib
{
    const double Re = 6378100;
    static double command_x;          /* Normal Overload Command in Lateral Plane */
    static double command_y;          /* Normal Overload Command in Longitudnal Plane */
    static double command_z;          /* Normal Overload Command in Longitudnal Plane */

    class Target :public Target_Template {
    private:

        // Target States
        std::string Name;
        double starttime{ 0 };
        TargetState Initial_State{ 0 };
        double Initial_Speed{ 0 }, Initial_Altitude{ 0 };
        TargetState State{ 0 };		                /* Interceptor Position */
        double R{ 0 };				                /* LOS Distance */
        double q_lat{ 0 }, q_lon{ 0 };              /* LOS Angle in Lateral and Longitudnal Plane*/
        double lead_lat_m{ 0 }, lead_lon_m{ 0 };    /* Leading Angles in Lateral and Longitudnal Planes */
        Recorder Rec1;                              /* To record the states of the Interceptor */
        double Vm{ 0 };                             /* Interceptor Speed */
        int phase{ 0 }; 
        double x_dist{ 0 }, y_dist{ 0 };            /* Distance covered along x, y and z-axis */
        double dist{ 0 }, dist_lim{ 0 };            /* Distance covered by the interceptor and its limit */
        double head_ang{ 0 }, path_ang{0};          /* Heading Angle */
        int miss_timer{ 0 }, miss_phase{ 0 };
        int selected_guidance{ 0 };
        double control_energy{ 0 };     /* Accumulated maneuver effort or control effort for the interception calculated as:  */
                                        /* 2-Norm of the normal overload curve in the lateral plane*/

        // Destination States
        TargetState Last_Target_State{ 0 };
        double lead_lat_t{ 0 }, lead_lon_t{ 0 };    /* Leading Angles of Target in Lateral and Longitudnal Planes */
        std::vector<TargetState> WayPoints;
        int num_wp{ 0 };
        std::vector<TargetState> Destinations{ 0 };
        bool more_data{ true };
        TargetState Current_Destination{ 0 };
        

        // Measurement Data
        double q_lat_dot{ 0 }, q_lon_dot{ 0 };

        // Parameters of Guidance Laws
        double command_max{ 0 };
        double time_acc{ 0 };                       /* The Time to Accelerate */
        double n{ 0 };
        double speed_max{ 0 };
        double turn_max{ 0 };

        /* Parameters to start and stop execution when destination reached */
        double R_stop{ 0 };
        double R_start{ 20e8 };
        bool save_csv = false;
       
        // Flags to Control the Execution
        bool my_terminated = false;
        bool is_reached = false;
        bool my_killed = false;
        bool self_destruct = false;
        bool valid_csv = false;
        bool is_first_run = true;
        
        /* Private Methods */
        void Decompose_R(double&, double&, double&);
        void Guidance_Law(TargetState);
        void MPNGuidance(TargetState);
        void APNGuidance(TargetState);
        void DGGuidance(TargetState);
        void LOSGuidance(TargetState);
        void PursuitGuidance(TargetState);
        void PurePNGuidance(TargetState);
        struct Course {				    /* Solves the dynamic state equations of Hypersonic Gliding Vehicle using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
        void Track(TargetState);                                /* Calculate the commands at time t_meas and track the Target */
        void Check_Param();
        void First_Run(TargetState);                            /* During the first execution calculate the required quantities */
        void Fine_Calculations(TargetState, TargetState);       /* Performs fine step calculations to detect the collision with greater accuracy */
        bool Is_Hit();
        bool Is_Missed(double);
        void Reset();
        bool Get_Destination(TargetState&);

        /* Public Methods from i_Interceptor */
        void Update_Target_State(TargetState);                  /* This function calculates all the required quantities for tracking from position measurements of Target and Interceptor */
        bool Get_State(TargetState&);
        void Get_Record();
        bool Is_Started();
        bool Is_Terminated(bool&, bool&, int&);

    public:
        Target() {};
        bool 	initialize(const std::string&, std::string&);
        bool 	reinitialize();
        bool 	get_state(double, std::vector<TargetState>&);
        bool    get_target_property(int xi_target_id, TargetProperty&);
        bool 	update_target(std::vector<TargetChange>&);
        //void 	get_performance_measures(TargetPerfMeasure&);
        int     get_expected_number_of_targets();
        double  get_start_time();
        bool    has_more_data();
    };
}

