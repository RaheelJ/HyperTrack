#pragma once

#include<iTruthGenerator.hpp>
#include<Matlab_dll.h>

struct PerfMeasure
{
	std::string        short_name;
	std::string        long_name;
	double             value;
};


struct TargetPerfMeasure
{
	std::vector<PerfMeasure> measures;
};

class iTarget :public iTruthGenerator
{
public:
	iTarget();
	~iTarget();

private:
	double my_last_update_time{ 0 };
	int my_num_targets{ 0 };
	project_detail_struct my_project_details;
	target_struct my_target_states;
	target_struct my_target_start_states;
	mwArray my_solution;
	int my_status;
	PerfMeasure dist_destination;
	PerfMeasure time_destination;

public:
	void initialize(int& status, std::string& message, project_detail_struct project_details, const char* config);
	void reinitialize();
	int get_expected_num_of_targets();
	double get_start_time();
	void get_targets(target_struct* targets, double current_time);
	void get_performance_measures(TargetPerfMeasure& xo_mop);
};