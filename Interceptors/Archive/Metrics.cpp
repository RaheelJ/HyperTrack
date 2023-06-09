#include "Metrics.hpp"
#include "xml/pugixml.hpp"

metrics::metrics()
{
    dist_LOS = 0;
    time_t = 0;
    dist_t = 0;
    dist_i = 0;
    time_i = 0;
    energy_i = 0;
    num_i = 0;
    kill_rate_i = 0;
}

metrics::~metrics()
{
}

double metrics::g_acc(double alt)   /* Calculate the integrrated gravitational acceleration as a function of altitude */
{
    double g0 = .8066498;           /* Gravitational acceleration at sea level (m/s) */
    double R = 6.378166e6;          /* Radius of Earth (m) */
    double tot_g;

    tot_g = g0 * (alt - alt * alt / R);
    return tot_g;
}

void metrics::calculate_all(iInterceptor& interceptor, iTarget& target)
{
    double details[7] = { 0 };
    double speed, altitude, dist_covered, time_taken, accuracy, num;
    
    //interceptor.get_details(details);
    speed = details[1];
    altitude = details[2];
    dist_covered = details[3];
    time_taken = details[4];
    accuracy = details[5];
    num = details[6];

    dist_LOS = details[0];
    dist_i = dist_covered;
    time_i = time_taken;
    num_i = num;
    kill_rate_i = accuracy;

    energy_i = speed * speed / 2 + g_acc(altitude) * altitude;              /* Sum of potential + kinetic energy (J) */
}
