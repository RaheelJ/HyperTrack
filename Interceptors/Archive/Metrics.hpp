#pragma once

#include "i_Interceptor.hpp"

class iTarget /* iTarget class requires C++ implementation based on MATLAB generated dll */
{
public:
	iTarget() {};
	~iTarget() {};
};

class metrics
{
	private:
		double g_acc(double alt);

	public:
		double dist_LOS;				/* Minimum distance between target and interceptor(s) */
		double dist_t;					/* Minimum distance between target and intended destination */ 
		double time_t;					/* Time to destination (by target) */
		double dist_i;					/* Distance traveled by interceptor(s) */
		double time_i;					/* Travel time by interceptor */
		double energy_i;				/* Energy expended by interceptors per unit mass */
		double num_i;					/* Number of expended interceptor(s) */ 
		double kill_rate_i;				/* Kill rate based on previous data */

		metrics();
		~metrics();

		void calculate_all(iInterceptor& interceptor, iTarget& target);

};
