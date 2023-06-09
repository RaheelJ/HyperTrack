#include "i_target.hpp"
#include "xml/pugixml.hpp"
#include "Target_Library.h"

iHypersonicTarget::iHypersonicTarget() :my_ptr(NULL)
{
}

iHypersonicTarget::~iHypersonicTarget()
{
    if (my_ptr) {
        delete my_ptr;
    }
}

bool iHypersonicTarget::initialize(const std::string& xi_init_file_name, std::string& xo_message)
{
    my_ptr = new TargetLib::Target();
    if (!my_ptr) {
        return false;
    }
    return (*my_ptr).initialize(xi_init_file_name, xo_message);
}

bool iHypersonicTarget::reinitialize()
{
    if (!my_ptr) {
        return false;
    }
    (*my_ptr).reinitialize();

    return true;
}

bool iHypersonicTarget::get_state(double xi_time, std::vector<TargetState>& xo_state)
{
    if (!my_ptr)
    {
        return false;
    }
    return (*my_ptr).get_state(xi_time, xo_state);
}

bool iHypersonicTarget::get_target_property(int xi_target_id, TargetProperty& xo_property)
{
    if (!my_ptr)
    {
        return false;
    }
    return (*my_ptr).get_target_property(xi_target_id, xo_property);
}

bool iHypersonicTarget::update_target(std::vector<TargetChange>& xo_changes)
{
    if (!my_ptr)
    {
        return false;
    }
    return (*my_ptr).update_target(xo_changes);
}

int iHypersonicTarget::get_expected_number_of_targets() {
    if (!my_ptr)
    {
        return false;
    }
    return (*my_ptr).get_expected_number_of_targets();
}

double iHypersonicTarget::get_start_time()
{
    if (!my_ptr)
    {
        return false;
    }
    return (*my_ptr).get_start_time();
}

bool iHypersonicTarget::has_more_data()
{
    if (!my_ptr)
    {
        return false;
    }
    return (*my_ptr).has_more_data();
}
