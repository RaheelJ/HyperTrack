#pragma once

#include <iTarget.hpp>
#include <Matlab_dll.h>

double rad_to_deg(double in);
double deg_to_rad(double in);
bool Calc_Distance(object_state_struct, object_state_struct, double&, double&);
void Matlab_dllRun(int&, std::string&, mwArray&, target_struct*, target_struct*, int, const char*, mwArray, double);
int run_main(int&, std::string&, mwArray&, target_struct*, target_struct*, int, const char*, mwArray, double);