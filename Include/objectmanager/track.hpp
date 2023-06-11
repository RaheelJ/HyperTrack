#ifndef __TRACK_H__
#define __TRACK_H__

#pragma once


//BEGIN_INTEGRATION_NAMESPACE
typedef std::string string_type;
typedef double value_type;

struct track
{
    int     id;
	int     display_id; 

    value_type time;

    bool    pos_in_geo;

    value_type pos_x_or_lat;
    value_type pos_y_or_long;
    value_type pos_z_or_alt;

    value_type vel_x;
    value_type vel_y;
    value_type vel_z;

    value_type pos_var_x;
    value_type pos_var_y;
    value_type pos_var_z;

    value_type pos_cov_xy;
    value_type pos_cov_xz;
    value_type pos_cov_yz;

    value_type vel_var_x;
    value_type vel_var_y;
    value_type vel_var_z;

    value_type vel_cov_xy;
    value_type vel_cov_xz;
    value_type vel_cov_yz;

	value_type cov_pos_x_vel_x;
	value_type cov_pos_x_vel_y;
	value_type cov_pos_x_vel_z;

	value_type cov_pos_y_vel_x;
	value_type cov_pos_y_vel_y;
	value_type cov_pos_y_vel_z;

	value_type cov_pos_z_vel_x;
	value_type cov_pos_z_vel_y;
	value_type cov_pos_z_vel_z;	
	

    track()
        : id(),          
		  display_id(),
          time(),
          pos_in_geo(),
          pos_x_or_lat(),
          pos_y_or_long(),
          pos_z_or_alt(),
          vel_x(),
          vel_y(),
          vel_z(),
          pos_var_x(),
          pos_var_y(),
          pos_var_z(),
          pos_cov_xy(),
          pos_cov_xz(),
          pos_cov_yz(),
          vel_var_x(),
          vel_var_y(),
          vel_var_z(),
          vel_cov_xy(),
          vel_cov_xz(),
          vel_cov_yz(),
		  cov_pos_x_vel_x(),
		  cov_pos_x_vel_y(),
		  cov_pos_x_vel_z(),		
		  cov_pos_y_vel_x(),
		  cov_pos_y_vel_y(),
		  cov_pos_y_vel_z(),
		  cov_pos_z_vel_x(),
		  cov_pos_z_vel_y(),
		  cov_pos_z_vel_z()
    {
    }
};

//END_INTEGRATION_NAMESPACE

#endif // __TRACK_H__
