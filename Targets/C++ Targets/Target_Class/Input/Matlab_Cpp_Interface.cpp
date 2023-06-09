#define _USE_MATH_DEFINES			/* Adds math constants */ 

#include <iostream>
#include <Matlab_Cpp_Interface.h>
#include <GeographicLib/Geodesic.hpp>
#include <cmath>             

using namespace GeographicLib;

double rad_to_deg(double in) {
	return 180 * fmod((in), 2 * M_PI) / M_PI;
}
double deg_to_rad(double in) {
	return M_PI * fmod((in), 360) / 180;
}

bool Calc_Distance(object_state_struct Target, object_state_struct Destination, double& xo_dist, double& xo_angle) {
	Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
	double angle2;
	geod.Inverse(Target.pos_x_or_lat, Target.pos_y_or_long, Destination.pos_x_or_lat, Destination.pos_y_or_long, xo_dist, xo_angle, angle2);
	return true;
}

void Matlab_dllRun(int& status, std::string& message, mwArray& solution_out, target_struct* start_states, target_struct* targets,
				   int FID, const char* config_file, mwArray solution_in, double current_time)
{

	// Input:FID
	mxDouble FIDInData = (double) FID;
	mwArray FIDIn(FIDInData);

	// Input:config_file
	const char* config_fileInData = config_file;
	mwArray config_fileIn(config_fileInData);

	// Input:solution_in
	mwArray solutionIn = solution_in;

	// Input:current_time
	mxDouble current_timeInData = current_time;
	mwArray current_timeIn(current_timeInData);

	// Outputs of Matlab_dll function
	mwArray statusOut;
	mwArray messageOut;
	mwArray solutionOut;
	mwArray start_statesOut;
	mwArray targetsOut;

	try
	{
		Matlab_dll(5, statusOut, messageOut, solutionOut, start_statesOut, targetsOut, FIDIn, config_fileIn, solutionIn, current_timeIn);
	}
	catch (const mwException& e)
	{
		std::cerr << e.what() << std::endl;
		return;
	}
	catch (...)
	{
		std::cerr << "Unexpected error thrown" << std::endl;
		return;
	}

	status = (int) statusOut;
	message = (std::string)(messageOut.ToString());
	solution_out = solutionOut;

	targets[0].pos_in_geo = 1;
	targets[0].pos_x_or_lat = targetsOut("pos_x_or_lat", 1);
	targets[0].pos_y_or_long = targetsOut("pos_y_or_long", 1);
	targets[0].pos_z_or_alt = targetsOut("pos_z_or_alt", 1);
	targets[0].vel_x = targetsOut("vel_x", 1);
	targets[0].vel_y = targetsOut("vel_y", 1);
	targets[0].vel_z = targetsOut("vel_z", 1);

	start_states[0].pos_in_geo = 1;
	start_states[0].pos_x_or_lat = start_statesOut("pos_x_or_lat", 1);
	start_states[0].pos_y_or_long = start_statesOut("pos_y_or_long", 1);
	start_states[0].pos_z_or_alt = start_statesOut("pos_z_or_alt", 1);
	start_states[0].vel_x = start_statesOut("vel_x", 1);
	start_states[0].vel_y = start_statesOut("vel_y", 1);
	start_states[0].vel_z = start_statesOut("vel_z", 1);
}

int run_main(int& status, std::string& message, mwArray& solution_out, target_struct* start_states, target_struct* targets,
			 int FID, const char* config_file, mwArray solution_in, double current_time)
{
	if (!Matlab_dllInitialize())
	{
		std::cerr << "Could not initialize the library properly" << std::endl;
		return -2;
	}
	else
	{
		Matlab_dllRun(status, message, solution_out, start_states, targets, FID, config_file, solution_in, current_time);
		// Call the application and library termination routine
		Matlab_dllTerminate();
	}
	return 0;
}


