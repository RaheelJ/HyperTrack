#pragma

#include <stdlib.h>
#include <string>

struct project_detail_struct
{
	double reference_time;					// Reference time in seconds from 1970/01/01
    double lrf_latitude;					// Local reference frame latitude
    double lrf_longitude;					// Local reference frame longitude
    double lrf_altitude;					// Local reference frame altitude

	project_detail_struct()
	{
		reference_time = 0;
		lrf_latitude = 0;
		lrf_longitude = 0;
		lrf_altitude = 0;
	}
};

struct target_transmitters_struct
{
	bool ais;							//0: no AIS transmiiter, 1: has AIS transmitter
    bool ads_b;							//0: no ADS-B transmiiter, 1: has ADS-B transmitter
    bool esm;							//0: no ESM transmiiter, 1: has ESM transmitter

	target_transmitters_struct() 
	{
		ais = 0;
		ads_b = 0;
		esm = 0;
	}
};

struct target_struct
{
	int id;						//target ID
	char* name;					//target name
	char* type;					//target type; options: aircraft,ship,boat,satellite,submarine,helicopter,car,truck,person,missile,other
    bool is_alive;				//1: alive, 0: dead 
    double time;				//time from reference time
    bool pos_in_geo;			//specifies if the target state in GEO or LRF (local reference frame). "pos_in_geo = 0" means LRF, "pos_in_geo = 1" means GEO
    double pos_x_or_lat;		//target x (in m) if "pos_in_geo = 0, otherwise target latitude (in deg)
    double pos_y_or_long;		//target y (in m) if "pos_in_geo = 0, otherwise target longitude (in deg)
    double pos_z_or_alt;		//target z (in m) if "pos_in_geo = 0, otherwise target altitude (in m)
    double vel_x;				//target velocity in the x/east direction (in m/s)
    double vel_y;				//target velocity in the y/north direction (in m/s)
    double vel_z;				//target velocity in the z/up direction (in m/s)
    target_transmitters_struct	transmitters;

	target_struct()
	{
		id = 0;
		name = NULL;
		type = NULL;
		is_alive = true;
		time = 0;
		pos_in_geo = 1;
		pos_x_or_lat = 0;
		pos_y_or_long = 0;
		pos_z_or_alt = 0;
		vel_x = 0;
		vel_y = 0;
		vel_z = 0;
	}
};

struct object_state_struct
{
	double time;					//time from reference time
    bool pos_in_geo;				//specifies if the object position in GEO or LRF (local reference frame). "pos_in_geo = 0" means LRF, "pos_in_geo = 1" means GEO
    double pos_x_or_lat;			//position x (in m) if "pos_in_geo = 0, otherwise object latitude (in deg)
    double pos_y_or_long;			//position y (in m) if "pos_in_geo = 0, otherwise object longitude (in deg)
    double pos_z_or_alt;			//position z (in m) if "pos_in_geo = 0, otherwise object altitude (in m)
    double vel_x;					//object velocity in the x/east direction (in m/s)
    double vel_y;					//object velocity in the y/north direction (in m/s)
    double vel_z;					//object velocity in the z/up direction (in m/s)

	object_state_struct()
	{
		time = 0;
		pos_in_geo = 1;
		pos_x_or_lat = 0;
		pos_y_or_long = 0;
		pos_z_or_alt = 0;
		vel_x = 0;
		vel_y = 0;
		vel_z = 0;
	}
};

struct waypoint_struct
{
	double latitude;
	double longitude;
	double altitude;
    double rest_time;
    double leaving_speed;

	waypoint_struct()
	{
		latitude = 0;
		longitude = 0;
		altitude = 0;
		rest_time = 0;
		leaving_speed = 0;
	}
};

struct target_change_struct
{
	int id;
	bool internal_target;		//0: target comes fom outside, 1: target created here
	char* name;
	bool start;
	bool terminate;
	bool change_start_time;
    bool change_end_time;
	double start_time;
	double end_time;
	bool change_heading;
	double new_heading;
	double new_vertical_heading;
	bool change_speed;
	double new_speed;
	bool change_currrent_state;
	object_state_struct current_state;
	bool change_waypoints;
	waypoint_struct* waypoints;

	target_change_struct()
	{
		id = 0;
		internal_target = 0;
		name = NULL;
		start = 0;
		terminate = 0;
		change_start_time = 0;
		change_end_time = 0;
		start_time = 0;
		end_time = 0;
		change_heading = 0;
		new_heading = 0;
		new_vertical_heading = 0;
		change_speed = 0;
		new_speed = 0;
		change_currrent_state = 0;
		change_waypoints = 0;
		waypoints = NULL;
	}
};


class iTruthGenerator
{
public:
	iTruthGenerator() {};
	
 /* project_detail_struct 				project_det;
	target_transmitters_struct*			targets_trans;
	target_struct*						targets;
	object_state_struct*				objects;
	waypoint_struct*					waypoints;
	target_change_struct*				targets_change; */

	// Initialize the truth-generator with given configuration file 
	// project_details is a struct defined under, create project_detail_struct()        
	// config file is the project file (xml, m, or any user specified file). It can be empty too.        
	// return status(1: success, 0: false), and error_message (if status = 0) 
	void initialize(int& status, std::string& message, project_detail_struct project_details, const char* config) {};

	// Reinitialize the truth-generator
	// Reinitialize is same as running the project again. So, except parameter setting, all the data may need to be cleared.
	void reinitialize() {};

	// Get the expected (or exact) number of targets, used for memory allocation
	int get_expected_num_of_targets() {};
	
	// Get the start time of the first target
	double get_start_time() {};
	
	// Get the target states of all the targets at the given time return array of targets as an array of struct defined under create_target_struct()
	void get_targets(target_struct* targets, double current_time) {};
	
	// Apply the changes, changes is an array of struct defined as target_change_struct        
	void apply_changes(target_change_struct* changes) {};
};