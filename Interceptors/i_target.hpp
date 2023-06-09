/////////////////////////////////////////////////////////////////////////
#ifndef _IHYPERSONIC_TARGET_H
#define _IHYPERSONIC_TARGET_H

#if defined(_WINDOWS) && defined(_APP_MODE_DLL)
#define DLLDIR __declspec(dllexport)
#else
#define DLLDIR
#endif

#include <string>
#include <vector>

//struct PerfMeasure
//{
//	std::string        short_name;
//	std::string        long_name;
//	double             value;
//};
//
//struct TargetPerfMeasure
//{
//	std::vector<PerfMeasure> measures;
//};

struct TargetChange
{
	int 	id;
	bool 	terminate;
};

enum etype {
	etype_missile,
	etype_hypersonic_vehicle
};

struct TargetProperty
{
	int 	     id;
	std::string  name;
	etype        type;
};

struct TargetState
{
	int		id;

	double	time;
	bool    state_in_geo;

	double  x_or_lat;
	double  y_or_lon;
	double  z_or_alt;

	double  x_vel;
	double  y_vel;
	double  z_vel;

	bool    is_alive = true;
};

class Target_Template; 

class DLLDIR iHypersonicTarget
{
private:

public:
	iHypersonicTarget();
	~iHypersonicTarget();

	bool 	initialize(const std::string& xi_init_file_name, std::string& xo_message);
	bool 	reinitialize();
	
	bool 	get_state(double xi_time, std::vector<TargetState>& xo_state);
	bool    get_target_property(int xi_target_id, TargetProperty& xo_property);
	bool 	update_target(std::vector<TargetChange>& xo_changes);
	//void 	get_performance_measures(TargetPerfMeasure& xo_mop);

	int     get_expected_number_of_targets();
	double  get_start_time();

	bool	has_more_data();

private:
	Target_Template* my_ptr;	
};

#endif // !_IHYPERSONIC_TARGET_H

