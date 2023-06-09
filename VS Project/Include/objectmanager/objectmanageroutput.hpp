#ifndef __OBJECTMANAGEROUTPUT_H__
#define __OBJECTMANAGEROUTPUT_H__

#pragma once

#include<map>

//BEGIN_INTEGRATION_NAMESPACE
typedef std::string string_type;
typedef double value_type;
typedef  int etarget_type;

struct platformstate
{
	value_type                 time;
	bool                       state_in_geo;	
	
	value_type                 x_or_lat;
	value_type                 y_or_lon;
	value_type                 z_or_alt;

	value_type                 x_vel;
	value_type                 y_vel;
	value_type                 z_vel;	

	//objectattitude             attitude;

	platformstate()
		: time(),
		state_in_geo(),		
		x_or_lat(),
		y_or_lon(),
		z_or_alt(),
		x_vel(),
		y_vel(),
		z_vel()
	{
	}
};

struct omposition
{
	value_type         latitude;
	value_type         longitude;
	value_type         altitude;
};

struct omvelocity
{
	value_type         east;
	value_type         north;
	value_type         up;
};

struct omwaypoint
{
	omposition        point;
	value_type        rest_time = 0;
	value_type        leaving_speed;
};

struct omnewsensor
{
	//sensordetails      sensor_details;
	bool               on_platform;
	string_type        platform_name;
	bool               sensor_on_off_with_platform = false;
};

struct omsensoronoff
{	
	bool                             on;
	value_type                       time;
};

struct omsensorchange
{
	bool               change_sampling_time = false;
	value_type         new_sampling_time;

	bool               set_next_sample_time = false;
	value_type         next_sample_time;

	bool               change_status = false; // change on/off
	bool               on_sensor;
	value_type         status_change_time;

	bool               change_search_pattern = false;
	//get search pattern    

	bool               change_coverage = false;
	//std::vector<measurementcoverage> new_coverage;
};

struct omsensorchanges
{
	std::map<string_type, omnewsensor>         new_sensors;
	std::map<string_type, omsensorchange>      changes;	
};


/*
struct omplatformpath
{
	value_type                  start_time;
	std::vector<omwaypoint>     waypoints;
};

struct omplatformchange
{
	bool                        change_path = false;
	omplatformpath              path;
};*/

struct omtargetchange
{
	int            id;
	string_type    name;
	bool           start = false;
	bool           terminate = false;
	bool           change_start_time = false;
	bool           change_end_time = false;
	value_type     start_time;
	value_type     end_time;
	bool           change_heading = false;
	value_type     new_heading = 0.0;
	value_type     new_vertical_heading = 0.0;
	bool           change_speed = false;
	value_type     new_speed;
	bool           change_currrent_state = false;
	platformstate  new_future_states;
	bool           change_waypoints;
	std::vector<omwaypoint>     waypoints;
};

struct omnewtarget
{
	int                         id;
	string_type                 name;
	etarget_type                type;

	value_type                  start_time;
	std::vector<omwaypoint>     path;
};

struct omplatformchanges
{
	std::map<string_type, omnewtarget>        new_platforms;
	std::map<string_type, omtargetchange>     changes;
};

struct omtargetchanges
{
	std::vector<omnewtarget>         new_targets;
	std::vector<omtargetchange>      changes;

	std::map<int, string_type>    target_display_message;
};

struct omsubmanagerchanges
{
	string_type                submanager_name;
		
	std::vector<int>           assigned_tracks;
	std::vector<omposition>    assigned_region;//polygon
};

struct omaoichanges
{
	string_type       name;
	string_type       description;
	value_type        start_time;
	value_type        end_time;
	//sensorcoverage    regions;	
};

struct omoutput
{
	std::vector<omsensorchanges>       sensor_changes;
	std::vector<omplatformchanges>     platform_changes;
	std::vector<omtargetchanges>       target_changes;
	std::vector<omaoichanges>          AOI_changes;

	std::vector<omsubmanagerchanges>   submanager_changes;
};

//END_INTEGRATION_NAMESPACE

#endif // __OBJECTMANAGEROUTPUT_H__
