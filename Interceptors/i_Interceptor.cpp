#include "Hypersonic_Interceptor.h"
#include "Fractional_Interceptor.h"
#include "NonLinear_Proportional_Interceptor.h"
#include "Coordinated_Interceptor.h"
#include "Coordinated_PNG_Interceptor.h"
#include "Coordinated_TPNG_Interceptor.h"
#include "SMC_Interceptor.h"
#include "xml/pugixml.hpp"

iInterceptor::iInterceptor():my_ptr(NULL)
{
}

iInterceptor::~iInterceptor()
{
    if (my_ptr) {
        delete my_ptr;
    }
}

bool iInterceptor::initialize(const std::string& xi_init_file_name, std::string& xo_message)
{
    pugi::xml_document Doc;
    pugi::xml_parse_result Result = Doc.load_file(xi_init_file_name.c_str());	
    if (Result) {}
    else {
        xo_message = "File not found or formatting error!";
        return false;
    }
    pugi::xml_node Root = Doc.child("Parameters");                              /* Accessing Root Node */
    if (Root) {}
    else {
        xo_message = "Root Node (Parameters) not found!";
        return false;
    }

    string my_algo = Root.child("Selection_Execution").child_value("Method");
    int method = 1 * (int)(my_algo == "Hypersonic") + 2 * (int)(my_algo == "Fractional") + 3 * (int)(my_algo == "NonLinear_Proportional") + 4 * (int)(my_algo == "Coordinated") +
                 5 * (int)(my_algo == "Adaptive_SMC") + 6 * (int)(my_algo == "Coordinated_PNG") + 7 * (int)(my_algo == "Coordinated_TPNG");
    switch (method) {
        case 1:
            my_ptr = new Hypersonic::Interceptor();
            break;
        case 2:
            my_ptr = new Fractional::Interceptor;
            break;
        case 3:
            my_ptr = new NonLinear_Proportional::Interceptor;
            break;
        case 4:
            my_ptr = new Coordinated::Interceptor;
            break;
        case 5:
            my_ptr = new SMC::Interceptor();
            break;
        case 6:
            my_ptr = new Coordinated_PNG::Interceptor();
            break;
        case 7:
            my_ptr = new Coordinated_TPNG::Interceptor();
            break;
        default:
            xo_message = "Invalid Method!";
            return false;
            break;
    }

    if (!my_ptr) {
        return false;
    }
    return (*my_ptr).Initialize(xi_init_file_name, xo_message);		
}

bool iInterceptor::reinitialize()
{
    if (!my_ptr) {
        return false;
    }
	my_assigned_target_id = -1;
    (*my_ptr).Reinitialize();
	
	return true;
}

void iInterceptor::update_target_state(const std::vector<objectstate>& xi_target)
{
	bool no_target = true;
	if (!xi_target.empty())
	{
		if (my_assigned_target_id < 0)
		{
			my_assigned_target_id = xi_target[0].id;
		}

		for (const auto& t : xi_target)
		{
			if (t.id == my_assigned_target_id)
			{				
				(*my_ptr).Update_Target_State(t);
				no_target = false;
				break;
			}
		}
	}

	if (no_target)
	{
		objectstate no_target{ 0 };
		no_target.id = -1;
		no_target.time = -1;
		(*my_ptr).Update_Target_State(no_target);
	}

    /*if (my_ptr) {
        if (!xi_target.empty())
        {
            (*my_ptr).Update_Target_State(xi_target[0]);
        }
        else
        {
            objectstate no_target{ 0 };
            no_target.id = -1;
            no_target.time = -1;
            (*my_ptr).Update_Target_State(no_target);
        }
    }*/
}

bool iInterceptor::get_state(objectstate& xo_state)
{
    if (!my_ptr)
    {
        return false;
    }
	return (*my_ptr).Get_State(xo_state);
}

bool iInterceptor::is_started(double xi_current_time)
{
    if (!my_ptr)
    {
        return false;
    }
	return (*my_ptr).Is_Started();
}

bool iInterceptor::is_terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
{	
    if (!my_ptr)
    {
        return false;
    }
	return (*my_ptr).Is_Terminated(xo_interceptor, xo_target, xo_target_id);		
}

void iInterceptor::get_csv() 
{
    if (my_ptr) {
        (*my_ptr).Get_Record();
    }
}

// void iInterceptor::get_details(double* details)
// {
    // if (my_ptr)
    // {
        // details[0] = (*my_ptr).dist_LOS;
        // details[1] = (*my_ptr).speed;
        // details[2] = (*my_ptr).altitude;
        // details[3] = (*my_ptr).dist_covered;
        // details[4] = (*my_ptr).time_taken;
        // details[5] = (*my_ptr).accuracy;
        // details[6] = (*my_ptr).num;
    // }
// }

void iInterceptor::get_performance_measures(InterceptorPerfMeasure& xo_mop)
{
	double g0 = 9.8066498;          /* Gravitational acceleration at sea level (m/s) */
    double R = 6.378166e6;          /* Radius of Earth (m) */
    double tot_g, energy_inst;
	double speed = (*my_ptr).speed;
	double alt = (*my_ptr).altitude;
	
	xo_mop.measures.resize(7);

	xo_mop.measures[0].short_name = "dist_LOS";
	xo_mop.measures[0].long_name = "Minimum distance between the interceptor and the target";
	xo_mop.measures[0].value = (*my_ptr).dist_LOS;
	
	xo_mop.measures[1].short_name = "dist_covered";
	xo_mop.measures[1].long_name = "Distance covered by the interceptor";
	xo_mop.measures[1].value = (*my_ptr).dist_covered;
	
	xo_mop.measures[2].short_name = "time_taken";
	xo_mop.measures[2].long_name = "Time taken to cover the distance";
	xo_mop.measures[2].value = (*my_ptr).time_taken;
	
	xo_mop.measures[3].short_name = "num";
	xo_mop.measures[3].long_name = "Number of interceptors";
	xo_mop.measures[3].value = (*my_ptr).num;
	
	xo_mop.measures[4].short_name = "accuracy";
	xo_mop.measures[4].long_name = "Kill rate of the interceptor";
	xo_mop.measures[4].value = (*my_ptr).accuracy;
	
    tot_g = g0 * (alt - alt * alt / R);						 /* Total Gravitational Acceleration */
    energy_inst = speed * speed / 2 + tot_g * alt;           /* Sum of potential + kinetic energy (J) per unit mass */
	xo_mop.measures[5].short_name = "energy_inst";
	xo_mop.measures[5].long_name = "Energy expended by the interceptor";
	xo_mop.measures[5].value = energy_inst;

    xo_mop.measures[6].short_name = "energy";
    xo_mop.measures[6].long_name = "Maneuver Effort of the Interceptor";
    xo_mop.measures[6].value = (*my_ptr).energy_i;
}
