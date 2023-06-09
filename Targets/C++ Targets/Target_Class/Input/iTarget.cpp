#include<iTarget.hpp>
#include<Matlab_Cpp_Interface.h>

// Constructor
iTarget::iTarget()
{
    my_last_update_time = -1;

    dist_destination.short_name = "dist_destination";
    dist_destination.long_name = "Minimum distance between the target and the destination";

    time_destination.short_name = "time_destination";
    time_destination.long_name = "Time required by the target to reach the destination";
}

iTarget::~iTarget()
{
    
}
        
// Initialization and Trajectory Optimization
void iTarget::initialize(int& status, std::string& message, project_detail_struct project_details, const char* config)
{
    my_project_details = project_details;

    // Declare outputs for the function
    mwArray solution_out;

    // Define inputs for the function;
    int FID = 1;        // 1 = Initialize
    mwArray solution_in;
    double current_time = 0;

    // Call the function
    run_main(status, message, solution_out, &my_target_start_states, &my_target_states, FID, config, solution_in, current_time);

    my_status = status;
    my_solution = solution_out;
    my_num_targets = 1;
}

// Reinitialze to the initial states
void iTarget::reinitialize()
{
    my_last_update_time = -1;
    my_target_states = my_target_start_states;
}

int iTarget::get_expected_num_of_targets()
{
    return my_num_targets;
}

double iTarget::get_start_time()
{
    return 0.0;
}

void iTarget::get_targets(target_struct* targets, double current_time)
{
    int status_mat  = 0;            // Status returned by Matlab dll
    int dim = 0;                    // Dimensions of my_solution
    double tf = 0;                  // Total time required to reach destination 
    mwArray temp_Array;
    object_state_struct dest;       // Destination Position
    object_state_struct tar;        // Current Target Position
    double dist, temp;              // Calculated distance between target and destination in the lateral plane

    // Declare outputs for the function
    int status = 0;
    std::string message = "Nil";
    mwArray solution_out;
    target_struct start_states;

    // Define inputs for the function;
    int FID = 2;        // 2 = Get Targets

    // Call the function
    if (my_status == 0 && current_time >= 0)
    {        
        // Calculate the performance metrics
        dim = my_solution.NumberOfElements();
        tf = my_solution("tf", dim);

        if (current_time >= tf)
        {
            return;
        }

        run_main(status, message, solution_out, &start_states, targets, FID, "Nil", my_solution, current_time);

        temp_Array = my_solution("X", dim);
        dest.pos_x_or_lat = rad_to_deg(temp_Array(78, 3));
        dest.pos_y_or_long = rad_to_deg(temp_Array(78, 2));
        tar.pos_x_or_lat = rad_to_deg(targets->pos_x_or_lat);
        tar.pos_y_or_long = rad_to_deg(targets->pos_y_or_long);
        Calc_Distance(tar, dest, dist, temp);
        
        dist_destination.value = sqrt(pow(dist, 2) + pow(targets[0].pos_z_or_alt, 2));
        time_destination.value = tf - current_time;
    }
    targets[0].pos_in_geo = 1;
    my_target_states =  targets[0];
}

void iTarget::get_performance_measures(TargetPerfMeasure& xo_mop)
{
    xo_mop.measures.resize(2);

    xo_mop.measures[0] = dist_destination;
    xo_mop.measures[1] = time_destination;
}


            