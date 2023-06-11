#ifndef __IOBJECTMANAGER_H__
#define __IOBJECTMANAGER_H__

#pragma once

//#include <integration/entity/tracker/track.hpp>
//#include <integration/entity/objectmanager/objectmanageroutput.hpp>
#include "objectmanager/objectmanageroutput.hpp"
#include "objectmanager/track.hpp"

typedef std::string string_type;
typedef double value_type;

struct sensorstate
{
	bool                       state_in_geo;

	value_type                 time;

	value_type                 x_or_lat;
	value_type                 y_or_lon;
	value_type                 z_or_alt;

	value_type                 x_velocity;
	value_type                 y_velocity;
	value_type                 z_velocity;

//	objectattitude             attitude;

	sensorstate()
		: state_in_geo(false),
		time(0),
		x_or_lat(0),
		y_or_lon(0),
		z_or_alt(0),
		x_velocity(0),
		y_velocity(0),
		z_velocity(0)
	{
	}
};

//BEGIN_INTEGRATION_NAMESPACE

struct iobjectmanagerinputparams
{
	string_type             name;
    string_type				config_file;
    int						coordinate_dimension;
/*	time_t                  reference_time;
	edistanceunit           distance_unit;
	etimeunit               time_unit;
	eangleunit              angle_unit;*/
	value_type              lrf_latitude;
	value_type              lrf_longitude;
	value_type              lrf_altitude;

	/*std::vector<sensordetails>   sensors;
	std::vector<platformdetails> platforms;
	std::vector<trackerdetails>  trackers;*/
};

struct iobjectmanager
{    
	typedef iobjectmanagerinputparams  iobjectmanagerinputparams_t;
	typedef sensorstate sensorstate_t;	
	typedef omoutput                          omoutput_t;

    typedef std::vector<track> tracks_t;

	virtual ~iobjectmanager()
    {
    }
	
	virtual bool initialize(const iobjectmanagerinputparams_t& xi_input_params, string_type& xo_message) = 0;
    virtual void uninitialize() = 0;
    virtual bool reinitialize() = 0;
	//virtual bool reinitialize(const iobjectmanagerinputparams_t& xi_input_params, string_type& xo_message) = 0;
	virtual bool prepare_for_new_run() = 0;
    //virtual bool change_parameters(string_type xi_config_file, string_type& xo_message) = 0;    

	//virtual void on_run_thread_open() {};
	//virtual void on_run_thread_close() {};
	
	virtual bool update_sensor_state(const string_type& xi_sensor_name, value_type xi_current_time, sensorstate_t& xi_state) { return true; };	
	//virtual bool add_targets(value_type xi_current_time, const targets_t& xi_targets) { return true; };
	virtual bool add_tracks(const string_type& xi_tracker_name, value_type xi_current_time, const tracks_t& xi_tracks) { return true; };
	//virtual bool add_measurements(const string_type& xi_sensor_name, const measurements_t &xi_measurements) { return true; };

	/*virtual bool manage_objects_offline(omoutput_t& xo_output)
	{
		return false;
	};*/
	virtual bool manage_objects_online(value_type xi_current_time, omoutput_t& xo_output)
	{
		return false;
	};	
};

//END_INTEGRATION_NAMESPACE

#endif // __IOBJECTMANAGER_H__
