#pragma once
#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include <string>
#include <ascent/Ascent.h>
#include <libInterpolate/Interpolate.hpp>
#include "i_Interceptor.hpp"
#include "Auxiliary.h"

using namespace asc;                    /* Namespace of Ascent library */
using namespace std;

namespace Hypersonic
{
    static double D;                    /* Drag Acceleration */
    static double g;                    /* Gravitational Acceleration */
    static double L;                    /* Lift of Interceptor */
    static double bank_ang;             /* Banking Angle */
    static double n_long;               /* Command of Reaction-Jet in Longitudnal Plane */
    static double n_lat;                /* Command of Reaction-Jet in Lateral Plane */

    class Interceptor :public Interceptor_Template {
    private:
        // Interceptor States
        objectstate Initial_State{ 0 };
        objectstate State{ 0 };	            /* Interceptor State */
        double head_ang{ 0 };	            /* Heading Angle */
        double path_ang{ 0 };               /* Flight Path Angle */
        double attack_ang{ 0 };			    /* Angle of Attack */
        double Vm{ 0 };                     /* Interceptor Speed */
        double x_dist{ 0 }, y_dist{ 0 };    /* Distance covered along x and y-axis */
        double energy{ 0 }, energy_lim{ 0 };     /* Energy consumed by the interceptor and its limit */
        Recorder Rec1;                      /* To record the states of the Interceptor */

        // Target State
        objectstate Last_Target_State{ 0 };

        // Measurement Data
        double delta_V{ 0 };		    /* Velocity Difference w.r.t LOS */
        double s{ 0 };				    /* Relative Distance */
        double delta_h{ 0 };            /* Altitude Difference */
        double q_dot_long{ 0 };         /* LOS Angle Rate in Longitudnal Plane */
        double q_long{ 0 };
        double q_dot_lat{ 0 };          /* LOS Angle Rate in Lateral Plane */
        double q_lat{ 0 };

        // Parameters for Boost Phase
        double boost_vel{ 0 }, boost_alt{ 0 };

        // Parameters for Aiming Phase
        double K1{ 0 }, K2{ 0 };

        // Parameters for calculation of Lift Coefficient
        double lift_coeff{ 0 };			 /* Lift Coefficient */
        double k1{ 0 }, k2{ 0 };         /* Constants for modified QEGC */
        double m{ 0 };                   /* Mass of Interceptor */
        double Sa{ 0 };                  /* Reference Area */
        double Re{ 0 };                  /* Radius of Earth */
        double atm_constant{ 0 };		 /* Atmospheric Density at Sea Level */
        double Hs{ 0 };					 /* Scale Altitude */

        // Interpolator for Angle of Attack calculation
        //_2D::ThinPlateSplineInterpolator<double> Interpolate;
        //_2D::BilinearInterpolator<double> Interpolate;

        /* Parameters for Velocity Error Corridor and Bank Reference*/
        double slope{ 0 }, bias{ 0 }, dead_band{ 0 };
        double gain{ 0 };
        int phase{ 0 };                                      /* 0 = Not Tracking, 1 = Boost, 2 = Glide, 3 = Aim, 4 = Hit & Kill */
        double s0{ 0 }, s1{ 0 }, s2{ 0 }, s3{ 0 };           /* Ranges for Relative Distance */

        /* Parameters for Reaction-Jet */
        double N_long{ 0 }, N_lat{ 0 };           /* Constant Loads in Longitudinal and Lateral Planes */
        double delta_q{ 0 }, delta_Q{ 0 };        /* ON and OFF Thresholds */

        /* Parameters to stop execution when Target hit */
        double hit_radius{ 0 };
        bool is_first_run = true;
        bool save_csv = false;

		double q_lat_last = 0;
		double q_long_last = 0;
		double t_last = 0;
        
        /* Private Methods */
        void Boost();                   /* Provides initial speed boost and altitude to the interceptor */
        double V_Lim(); 			    /* Returns  Velocity Bound for Banking Commands based on s */
        double Bank_Ref(); 			    /* Returns reference for Banking Command based on s */
        void Calc_Bank();			    /* Calulates Banking Command for Gliding Phase*/
        double Atm_Den();   		    /* Returns Atmospheric Density */
        void Calc_Lift(); 			    /* Calculate Lift using Modified QEGC and then Calulate Lift Coefficient */
        void Calc_AoA();    			/* Calculate Angle of Attack from Lift Coefficient */
        void Glide();
        void Aim(); 				    /* Calculate Banking and Lift commands using PNG in the Aiming mode */
        void Jet_Law(double, double, double&); 	/* Switching Control Law of the Reaction-Jet */
        void Kill(); 				    /* Change Velocity direction using Reaction Jet in Hit to Kill Phase */
        struct Course {				    /* Solves the dynamic state equations of Hypersonic Gliding Vehicle using input commands  */
            void operator()(const state_t&, state_t&, const double);
        };
		bool Is_Hit();
        bool Is_Missed(int);
		void Track(double);                                     /* Calculate the commands at time t_meas and track the Target */
        void Reset();
        void Fine_Calculations(objectstate, objectstate);       /* Performs fine step calculations to detect the collision with greater accuracy */
        void First_Run(objectstate);                            /* During the first execution calculate the required quantities */
        void Check_Param();                                     /* Check system parameters for any abnormal conditions */
        bool my_terminated = false;
        bool my_killed = false;
        bool self_destruct = false;
        bool valid_csv = false;

    public:
        Interceptor() {};
        bool Initialize(const string&, string&);                /* Initializes all the necessary parameters via a "Parameter.xml" file */
        void Reinitialize();
        void Update_Target_State(objectstate);                  /* This function calculates all the required quantities for tracking from position measurements of Target and Interceptor */        
		bool Get_State(objectstate&);
        void Get_Record();
        bool Is_Started();       
		bool Is_Terminated(bool&, bool&, int&);        
    };
}

