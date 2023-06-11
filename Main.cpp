#include<Interceptor/Fractional_Interceptor.h>
#include<Interceptor/Hypersonic_Interceptor.h>
#include<Interceptor/Coordinated_PNG_Interceptor.h>
#include<Interceptor/Coordinated_TPNG_Interceptor.h>
#include<Interceptor/NonLinear_Proportional_Interceptor.h>
#include<Interceptor/MPC_Interceptor.h>
#include<iostream>
#include<Interceptor/Auxiliary.h>
#include<Interceptor/i_Interceptor.hpp>
#include<Interceptor/i_target.hpp>
#include<Interceptor/Target_Library.h>
#include<Interceptor/Guidance_Library.h>

int main() {
    //iTarget T1;
    //metrics M1;
    double t_meas = 0;                                      /* Start time of Target flight */
    iInterceptor I1;                                        /* Initialize the Interceptor */
    objectstate Target_State;
    std::vector<objectstate> Target_Vector(1);
    objectstate Interceptor_State;
    Interceptor_State.id = 0;
    bool tracking, init, interceptor_terminate, target_hit;
    int target_id;
    string message;
    InterceptorPerfMeasure M1;

    init = I1.initialize("Input/Case_Study_BGV/NonLinearProportional_Interceptor_Halifax_d12.xml", message);
    if (init) {}
    else {
        cout << message;
        return 0;
    }

    /* Loop for tracking */
    for (t_meas = 400; t_meas <= 1000; t_meas = t_meas + 1) {
        Read_CSV(t_meas, "Input/Case_Study_BGV/BGV_Hypersonic_Moscow.csv", Target_State);           /* Read the Target state */
        Target_Vector[0] = Target_State;
        cout << Target_State.x_or_lat << '\t' << Target_State.y_or_lon << '\t' << Target_State.z_or_alt << '\t' << Target_State.time << "\n";

        I1.update_target_state(Target_Vector);                            /* Follow the Target and save the Trajectory */

        tracking = I1.is_started(t_meas);
        if (tracking)
        {
            cout << "In Pursuit" << "\n";
        }

        I1.get_state(Interceptor_State);

        cout << Interceptor_State.x_or_lat << '\t' << Interceptor_State.y_or_lon << '\t' << Interceptor_State.z_or_alt << '\t' << Interceptor_State.time << "\n";
        cout << Interceptor_State.x_vel << '\t' << Interceptor_State.y_vel << '\t' << Interceptor_State.z_vel << "\n \n";

        I1.is_terminated(interceptor_terminate, target_hit, target_id);
        if (target_hit) {
            cout << '\n' << "Target Hit !!" << '\n';
            break;
        }
        else if (interceptor_terminate) {
            cout << '\n' << "Target Missed !!" << '\n';
            break;
        }
    }
    I1.get_csv();
    I1.get_performance_measures(M1);
    cout << "Distance Covered: " << M1.measures[1].value << '\n';
    cout << "Time Taken: " << M1.measures[2].value << '\n';
    cout << "Energy Expended: " << M1.measures[5].value << '\n';
    cout << "Accuracy: " << M1.measures[4].value << '\n';
    cout << "Number Expended: " << M1.measures[3].value << '\n';
    cout << "Distance to Target: " << M1.measures[0].value << "\n\n";

    return 0;
}

//int main()
//{
//	iHypersonicTarget T1;
//	int total_targets;
//	bool more_data, initialized, reinitialized, valid_state, got_property, target_updated;
//	double start_time;
//	std::string message;
//	std::vector<TargetState> TState;
//	TargetProperty TProperty;
//	std::vector<TargetChange> TChange;
//
//	initialized = T1.initialize("Input/Hypersonic_Target.xml", message);
//	total_targets = T1.get_expected_number_of_targets();
//	got_property = T1.get_target_property(1, TProperty);
//
//	for (double i = 0; i < 4000; i++)
//	{
//		if (T1.has_more_data())
//		{
//			valid_state = T1.get_state(i, TState);
//			if (valid_state)
//			{
//				cout << TState[0].x_or_lat << '\t' << TState[0].y_or_lon << '\t' << TState[0].z_or_alt << '\t' << TState[0].time << "\n";
//				cout << TState[0].x_vel << '\t' << TState[0].y_vel << '\t' << TState[0].z_vel << "\n \n";
//				TState.pop_back();
//			}
//		}
//		else
//			break;
//
//		/*if (i == 1400)
//		{
//			TargetChange temp_change;
//			temp_change.id = 1;
//			temp_change.terminate = true;
//			TChange.push_back(temp_change);
//			T1.update_target(TChange);
//		}*/
//	}
//}

//int main()
//{
//    double t_meas = 0;
//    int i, j;
//    MPC::Interceptor I1;
//    asc::Recorder R1;
//
//    I1.Reinitialize();
//    for (i = 0; i < 30; i++)
//    {
//        t_meas = i * 0.1;
//        cout << "Loop: " << i << "\t" << "Time: " << t_meas << std::endl << std::endl;
//       
//        if (t_meas <= 1)
//        {
//            I1.Ref << 0, 5, 5;
//        }
//        else if (t_meas <= 2)
//        {
//            I1.Ref << 0, 10, 10;
//        }
//        else
//        {
//            I1.Ref << 0, 5, 5;
//        }
//        cout << "Reference: " << I1.Ref.transpose() << std::endl << std::endl;
//
//        for (j = t_meas; j < t_meas + 0.1; j = j + I1.Ts)
//        {
//            MPC::d << 0.2 * sin(2 * j), 0.2 * sin(2 * j), 0.2 * sin(2 * j), sin(2 * j), sin(2 * j), sin(2 * j);
//            cout << "Disturbance: " << MPC::d.transpose() << std::endl << std::endl;
//
//            I1.FTDO_Run();
//            I1.NMP_Control(t_meas);
//            I1.DCA();
//            I1.NSI_Run();
//        }
//
//        R1({ MPC::v(0), I1.Ref(0), MPC::v(1), I1.Ref(1), MPC::v(2), I1.Ref(2), t_meas });
//    }
//
//    R1.csv("Output/MPC_States", { "Bank_Angle", "Bank_Angle_Ref", "Overload_Longitudnal", "Overload_Longitudnal_Reference", "Overload_Lateral", "Overload_Longitudnal_Reference" });
//}


