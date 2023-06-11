#include <iTarget.hpp>
#include <Matlab_Cpp_Interface.h>

int main()
{
	// Prepare to run Matlab runtime. Should be defined at the start of main function to use mwArray data structure. 
	if (!mclInitializeApplication(nullptr, 0))
	{
		std::cerr << "Could not initialize the application." << std::endl;
	}

	iTarget Tar1;
	int i, status;
	std::string message;
	project_detail_struct proj_details;
	double time = 0, time_step = 100;
	target_struct tar_position;
	TargetPerfMeasure xo_mop;

	std::cout << "Initializing ... " << std::endl;
	Tar1.initialize(status, message, proj_details, "Input/Output_ABH_1.xml");
	if (status == 0)
	{
		std::cout << message << std::endl;
		std::cout << "Target Initialized !" << std::endl;
	}
	else
	{
		std::cout << message << std::endl;
		std::cout << "Target not initialized !" << std::endl;
		return -1;
	}

	for (i = 0; i <= 50; i++)
	{
		time = i * time_step;
		Tar1.get_targets(&tar_position, time);
		std::cout << std::endl << "Target: " << rad_to_deg(tar_position.pos_x_or_lat) << "  " << rad_to_deg(tar_position.pos_y_or_long) << "  " << tar_position.pos_z_or_alt << std::endl;

		Tar1.get_performance_measures(xo_mop);
		std::cout << "Time to Destination: " << xo_mop.measures[1].value << std::endl;
		std::cout << "Distance to Destination: " << xo_mop.measures[0].value << std::endl;
	}

	// Call mclTerminateApplication at the end of the application to shut down all MATLAB Runtime instances.
	mclTerminateApplication();

	return 0;
}