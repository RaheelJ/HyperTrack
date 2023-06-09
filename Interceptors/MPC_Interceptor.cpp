#define _USE_MATH_DEFINES               /* Adds math constants */ 

#include "MPC_Interceptor.h"
#include <iostream>
#include "xml/pugixml.hpp"
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

void MPC::Interceptor::Reset() {
}
bool MPC::Interceptor::Initialize(const string& XML_Parameters, string& message) {
    State.x_or_lat = 0;
    State.y_or_lon = 0;
    State.z_or_alt = 20e3;
	return true;
}
void MPC::Interceptor::Reinitialize() {
    State.time = 0;
    State.x_or_lat = 0;
    State.y_or_lon = 0;
    State.z_or_alt = 32e3;

    bank_ang = deg_to_rad(5);
    attack_ang = deg_to_rad(5);
    slip_ang = deg_to_rad(-5);

    z_0[0] = bank_ang;
    z_0[1] = slip_ang;
    z_0[2] = attack_ang;

    double atm_den{ 0 };                                    // Atmospheric Density
    atm_den = atm_constant * exp(-State.z_or_alt / Hs);
    q = atm_den * pow(Vm, 2) / 2;

    h_x(0) = bank_ang;
    h_x(1) = q * S * C_z * slip_ang / m;
    h_x(2) = q * S * C_y * attack_ang / m;
    y = h_x;

    int i;
    for (i = 0; i < 3; i++)
    {
        u_max(i) = deg_to_rad(u_max(i));
        u_min(i) = deg_to_rad(u_min(i));
        u_rate_max(i) = deg_to_rad(u_rate_max(i));
        u_rate_min(i) = deg_to_rad(u_rate_min(i));
    }
}

void MPC::Interceptor::NSI_Model::operator()(const state_t& x, state_t& dx, const double) {
    /*  x[0] = bank_ang
        x[1] = slip_ang
        x[2] = attack_ang
        x[3] = rate_roll
        x[4] = rate_yaw
        x[5] = rate_pitch
        x[6] = pitch_ang
        x[7] = longitude
        x[8] = latitude
        x[9] = altitude
        x[10] = path_ang
        x[11] = head_ang
    */

    // Variable Parameters
    VectorXd vv{ {0, 0, 0} };            // Virtual Control Moment (Mx, My, Mz)
    vv = Bc * u;

    VectorXd f_x{ {0,0,0,0,0,0} };
    f_x(0) = x[5] * sin(x[0]) * tan(x[6]) - x[4] * cos(x[0]) * tan(x[6]) + x[3];
    f_x(1) = x[3] * sin(x[2]) * x[4] * cos(x[2]) + n_y * sin(x[2]) * sin(x[1]) / Vm + n_z * cos(x[1]) / Vm;
    f_x(2) = x[4] * tan(x[1]) * sin(x[2]) + x[5] - (n_y * cos(x[2])) / (Vm * cos(x[1])) - x[3] * tan(x[1]) * cos(x[2]);
    f_x(3) = ((J_y - J_z) / J_x) * x[4] * x[5];
    f_x(4) = ((J_z - J_x) / J_y) * x[3] * x[5] + q * S * L_ * m_y_beta * x[1] / J_y + q * S * L_ * m_y_w * x[4] / J_y;
    f_x(5) = ((J_x - J_y) / J_z) * x[3] * x[4] + q * S * L_ * m_z_alpha * x[2] / J_z + q * S * L_ * m_z_w * x[5] / J_z;

    MatrixXd g1_x = MatrixXd::Zero(6, 3);
    g1_x(3, 0) = 1 / J_x;
    g1_x(4, 1) = 1 / J_y;
    g1_x(5, 2) = 1 / J_z;

    MatrixXd g2_x = MatrixXd::Identity(6, 6);
    
    VectorXd NSI_Dynamics(6);
    NSI_Dynamics = f_x + g1_x * vv + g2_x * d;

    dx[0] = NSI_Dynamics(0);
    dx[1] = NSI_Dynamics(1);
    dx[2] = NSI_Dynamics(2);
    dx[3] = NSI_Dynamics(3);
    dx[4] = NSI_Dynamics(4);
    dx[5] = NSI_Dynamics(5);

    dx[6] = x[4] * sin(x[0]) + x[5] * cos(x[0]);
    dx[7] = Vm * cos(x[10]) * sin(x[11]) / ((x[9] + Re) * cos(x[8]));
    dx[8] = Vm * cos(x[10]) * cos(x[11]) / (x[9] + Re);
    dx[9] = Vm * sin(x[10]);
    dx[10] = x[5] - dx[2];
    dx[11] = x[4] + dx[1];
}
void MPC::Interceptor::NSI_Run() 
{
    state_t x = { 0,0,0,0,0,0,0,0,0,0,0,0 };        //Initial States for this iteration
    x[0] = bank_ang;
    x[1] = slip_ang;
    x[2] = attack_ang;
    x[3] = rate_roll;
    x[4] = rate_yaw;
    x[5] = rate_pitch;
    x[6] = pitch_ang;
    x[7] = State.y_or_lon;
    x[8] = State.x_or_lat;
    x[9] = State.z_or_alt;
    x[10] = path_ang;
    x[11] = head_ang;

    double t = State.time;
    double step_size = Ts;

    RK4 Integrator;
    Integrator(NSI_Model(), x, t, step_size);      /* Solve the ODEs to find system states */
    
    // States
    bank_ang = x[0];
    slip_ang = x[1];
    attack_ang = x[2];
    rate_roll = x[3];
    rate_yaw = x[4];
    rate_pitch = x[5];
    pitch_ang = x[6];
    State.y_or_lon = x[7];
    State.x_or_lat = x[8];
    State.z_or_alt = x[9];
    path_ang = x[10];
    head_ang = x[11];

    // Outputs, h(x)
    h_x(0) = bank_ang;
    h_x(1) = q * S * C_z * slip_ang / m;
    h_x(2) = q * S * C_y * attack_ang / m;
    y = h_x;
    cout << "Interceptor Output: " << y.transpose() << std::endl << std::endl;
}

void MPC::Interceptor::FTDO_Model::operator()(const state_t& x, state_t& dx, const double)
{
    VectorXd f_x{ {0,0,0,0,0,0} };
    f_x(0) = x_[5] * sin(x_[0]) * tan(x_[6]) - x_[4] * cos(x_[0]) * tan(x_[6]) + x_[3];
    f_x(1) = x_[3] * sin(x_[2]) * x_[4] * cos(x_[2]) + n_y * sin(x_[2]) * sin(x_[1]) / Vm + n_z * cos(x_[1]) / Vm;
    f_x(2) = x_[4] * tan(x_[1]) * sin(x_[2]) + x_[5] - (n_y * cos(x_[2])) / (Vm * cos(x_[1])) - x_[3] * tan(x_[1]) * cos(x_[2]);
    f_x(3) = ((J_y - J_z) / J_x) * x_[4] * x_[5];
    f_x(4) = ((J_z - J_x) / J_y) * x_[3] * x_[5] + q * S * L_ * m_y_beta * x_[1] / J_y + q * S * L_ * m_y_w * x_[4] / J_y;
    f_x(5) = ((J_x - J_y) / J_z) * x_[3] * x_[4] + q * S * L_ * m_z_alpha * x_[2] / J_z + q * S * L_ * m_z_w * x_[5] / J_z;

    MatrixXd g1_x = MatrixXd::Zero(6, 3);
    g1_x(3, 0) = 1 / J_x;
    g1_x(4, 1) = 1 / J_y;
    g1_x(5, 2) = 1 / J_z;

    /*  x[0] ... x[5] = z_0[0] ... z_0[5]
        x[6] ... x[11] = z_1[0] ... z_1[5]
        x[12] ... x[17] = z_2[0] ... z_2[5]
    */

    VectorXd L_0{ {0,0,0,0,0,0} }, L_1{ {0,0,0,0,0,0} }, L_2{ {0,0,0,0,0,0} };
    VectorXd Obs_Dynamics_1(6), Obs_Dynamics_2(6), Obs_Dynamics_3(6);
    int i;
    for (i = 0; i < 6; i++)
    {
        L_0(i) = pow(L[i], 1 / 3) * pow(abs(x[i] - x_(i)), 2 / 3) * signum(x[i] - x_(i));
        L_1(i) = pow(L[i], 1 / 2) * pow(abs(x[i + 6] - v_0(i)), 1 / 2) * signum(x[i + 6] - v_0(i));
        L_2(i) = L[i] * signum(x[i + 12] - v_1(i));

        v_0(i) = -lambda_0 * L_0(i) + x[i+6];
        v_1(i) = -lambda_1 * L_1(i) + x[i+12];
        v_2(i) = -lambda_2 * L_2(i);
    }

    Obs_Dynamics_1 = v_0 + f_x + g1_x * v;
    Obs_Dynamics_2 = v_1;
    Obs_Dynamics_3 = v_2;

    for (i = 0; i < 6; i++)
    {
        dx[i] = Obs_Dynamics_1(i);
        dx[i + 6] = Obs_Dynamics_2(i);
        dx[i + 12] = Obs_Dynamics_3(i);
    }
}
void MPC::Interceptor::FTDO_Run()
{
    int i;
    state_t x;
    x.resize(18, 0);

    double atm_den{ 0 };                                    // Atmospheric Density
    atm_den = atm_constant * exp(-State.z_or_alt / Hs);
    q = atm_den * pow(Vm, 2) / 2;
    n_y = q * S * C_y * attack_ang / m;
    n_z = q * S * C_z * slip_ang / m;

    Matrix3d Bc_1 = Matrix3d::Zero();
    Bc_1(0, 0) = q * S * L_ * m_x_delta;
    Bc_1(1, 1) = q * S * L_ * m_y_delta;
    Bc_1(2, 2) = q * S * L_ * m_z_delta;

    MatrixXd Bc_2 = MatrixXd::Zero(3, 2);
    Fsmax = 4 * 9.8 * m;
    Bc_2(1, 1) = (1 + K_M_y) * Fsmax * Lm;
    Bc_2(2, 0) = (1 + K_M_z) * Fsmax * Lm;

    Bc << Bc_1, Bc_2;
    cout << "Bc: " << std::endl << Bc << std::endl << std::endl;

    for (i = 0; i < 6; i++)
    {
        x[i] = z_0[i];
        x[i + 6] = z_1[i];
        x[i + 12] = z_2[i];
    }
    x_(0) = bank_ang;
    x_(1) = slip_ang;
    x_(2) = attack_ang;
    x_(3) = rate_roll;
    x_(4) = rate_yaw;
    x_(5) = rate_pitch;
    x_(6) = pitch_ang;

    double t = T;
    double step_size = Ts;

    RK4 Integrator;
    Integrator(FTDO_Model(), x, t, step_size);      /* Solve the ODEs to find system states */

    for (i = 0; i < 6; i++)
    {
        z_0[i] = x[i];
        z_1[i] = x[i + 6];
        z_2[i] = x[i + 12];

        Obs_d(i) = z_1[i];
        Obs_x[i] = z_0[i];
    }
    cout << "Observed Disturbance: " << Obs_d.transpose() << std::endl << std::endl;
}

void MPC::Interceptor::NMP_Control(double t_meas)
{
    //Lie Derivative Expressions of Roll Angle Subsystem
    double f_gamma = rate_pitch * sin(bank_ang) * tan(pitch_ang) - rate_yaw * cos(bank_ang) * tan(pitch_ang) + rate_roll;
    double f_beta = rate_roll * sin(attack_ang) + rate_yaw * cos(attack_ang) + n_y * sin(attack_ang) * sin(slip_ang) / Vm + n_z * cos(slip_ang) / Vm;
    double f_alpha = rate_yaw * tan(slip_ang) * sin(attack_ang) + rate_pitch - (n_y * cos(attack_ang)) / (Vm * cos(slip_ang)) - rate_roll * tan(slip_ang) * cos(attack_ang);
    double f_wx = ((J_y - J_z) / J_x) * rate_yaw * rate_pitch;
    double f_wy = ((J_z - J_x) / J_y) * rate_roll * rate_pitch + q * S * L_ * m_y_beta * slip_ang / J_y + q * S * L_ * m_y_w * rate_yaw / J_y;
    double f_wz = ((J_x - J_y) / J_z) * rate_roll * rate_yaw + q * S * L_ * m_z_alpha * attack_ang / J_z + q * S * L_ * m_z_w * rate_pitch / J_z;

    double dfgamma_dgamma = rate_yaw * sin(bank_ang) * tan(pitch_ang) + rate_pitch * cos(bank_ang) * tan(pitch_ang);
    double dfgamma_dbeta = 0;
    double dfgamma_dalpha = (-rate_yaw * cos(bank_ang) + rate_pitch * sin(bank_ang)) / pow(cos(attack_ang), 2);
    double dfgamma_dwx = 1;
    double dfgamma_dwy = -cos(bank_ang) * tan(pitch_ang);
    double dfgamma_dwz = sin(bank_ang) * tan(pitch_ang);

    double Lf_h1 = f_gamma;
    RowVectorXd Lg1_h1{ {0,0,0} };
    RowVectorXd Lg2_h1{ {1, 0, 0, 0, 0, 0} };
    double L2f_h1 = dfgamma_dgamma * f_gamma + dfgamma_dbeta * f_beta + dfgamma_dalpha * f_alpha + dfgamma_dwx * f_wx + dfgamma_dwy * f_wy + dfgamma_dwz * f_wz;

    RowVectorXd Lg1_Lf_h1(3);
    Lg1_Lf_h1(0) = (1 / J_x) * dfgamma_dwx;
    Lg1_Lf_h1(1) = (1 / J_y) * dfgamma_dwy;
    Lg1_Lf_h1(2) = (1 / J_z) * dfgamma_dwz;

    RowVectorXd Lg2_Lf_h1{ {0,0,0,0,0,0} };
    Lg2_Lf_h1(0) = dfgamma_dgamma;
    Lg2_Lf_h1(1) = dfgamma_dbeta;
    Lg2_Lf_h1(2) = dfgamma_dalpha;
    Lg2_Lf_h1(3) = dfgamma_dwx;
    Lg2_Lf_h1(4) = dfgamma_dwy;
    Lg2_Lf_h1(5) = dfgamma_dwz;

    //Lie Derivative Expressions of Vertical Overload Subsystem
    double dfbeta_dgamma = 0;
    double dfbeta_dbeta = (q * S * C_y * attack_ang * sin(attack_ang) * cos(slip_ang) + q * S * C_z * (cos(slip_ang) - slip_ang * sin(slip_ang))) / (m * Vm);
    double dfbeta_dalpha = rate_roll * cos(attack_ang) - rate_yaw * sin(attack_ang) + q * S * C_y * (sin(attack_ang) + attack_ang * cos(attack_ang)) * sin(slip_ang) / (m * Vm);
    double dfbeta_dwx = sin(attack_ang);
    double dfbeta_dwy = cos(attack_ang);
    double dfbeta_dwz = 0;

    double Lf_h2 = q * S * C_z * f_beta / (m * g);
    RowVectorXd Lg1_h2{ {0,0,0} };
    RowVectorXd Lg2_h2{{0,0,0,0,0,0}};
    Lg2_h2(1) = q * S * C_z / (m * g);

    double L2f_h2 = dfbeta_dgamma * f_gamma + dfbeta_dbeta * f_beta + dfbeta_dalpha * f_alpha + dfbeta_dwx * f_wx + dfbeta_dwy * f_wy + dfbeta_dwz * f_wz;
    L2f_h2 = L2f_h2 * (q * S * C_z / (m * g));

    RowVectorXd Lg1_Lf_h2(3);
    Lg1_Lf_h2(0) = (q * S * C_z / (m * g)) * (1 / J_x) * dfbeta_dwx;
    Lg1_Lf_h2(1) = (q * S * C_z / (m * g)) * (1 / J_y) * dfbeta_dwy;
    Lg1_Lf_h2(2) = (q * S * C_z / (m * g)) * (1 / J_z) * dfbeta_dwz;

    RowVectorXd Lg2_Lf_h2(6);
    Lg2_Lf_h2(0) = (q * S * C_z / (m * g)) * dfbeta_dgamma;
    Lg2_Lf_h2(1) = (q * S * C_z / (m * g)) * dfbeta_dbeta;
    Lg2_Lf_h2(2) = (q * S * C_z / (m * g)) * dfbeta_dalpha;
    Lg2_Lf_h2(3) = (q * S * C_z / (m * g)) * dfbeta_dwx;
    Lg2_Lf_h2(4) = (q * S * C_z / (m * g)) * dfbeta_dwy;
    Lg2_Lf_h2(5) = (q * S * C_z / (m * g)) * dfbeta_dwz;

    //Lie Derivative Expressions of Vertical Overload Subsystem
    double dfalpha_dgamma = 0;
    double dfalpha_dbeta = rate_yaw * sin(attack_ang) / pow(cos(slip_ang), 2) - rate_roll * cos(attack_ang) / pow(cos(slip_ang), 2) - q * S * C_y * attack_ang * cos(attack_ang) * sin(slip_ang) / (m * Vm * pow(cos(slip_ang), 2));
    double dfalpha_dalpha = rate_yaw * tan(slip_ang) * cos(attack_ang) + rate_roll * tan(slip_ang) * sin(attack_ang) - q * S * C_y * (cos(attack_ang) - attack_ang * sin(attack_ang)) / (m * Vm * cos(slip_ang));
    double dfalpha_dwx = -tan(slip_ang) * cos(attack_ang);
    double dfalpha_dwy = tan(slip_ang) * sin(attack_ang);
    double dfalpha_dwz = 1;

    double Lf_h3 = q * S * C_y * f_alpha / (m * g);
    RowVectorXd Lg1_h3{ {0,0,0} };
    RowVectorXd Lg2_h3{ {0,0,0,0,0,0} };
    Lg2_h3(2) = q * S * C_y / (m * g);

    double L2f_h3 = dfalpha_dgamma * f_gamma + dfalpha_dbeta * f_beta + dfalpha_dalpha * f_alpha + dfalpha_dwx * f_wx + dfalpha_dwy * f_wy + dfalpha_dwz * f_wz;
    L2f_h3 = L2f_h3 * (q * S * C_y / (m * g));

    RowVectorXd Lg1_Lf_h3(3);
    Lg1_Lf_h3(0) = (q * S * C_y / (m * g)) * (1 / J_x) * dfalpha_dwx;
    Lg1_Lf_h3(1) = (q * S * C_y / (m * g)) * (1 / J_y) * dfalpha_dwy;
    Lg1_Lf_h3(2) = (q * S * C_y / (m * g)) * (1 / J_z) * dfalpha_dwz;

    RowVectorXd Lg2_Lf_h3(6);
    Lg2_Lf_h3(0) = (q * S * C_y / (m * g)) * dfalpha_dgamma;
    Lg2_Lf_h3(1) = (q * S * C_y / (m * g)) * dfalpha_dbeta;
    Lg2_Lf_h3(2) = (q * S * C_y / (m * g)) * dfalpha_dalpha;
    Lg2_Lf_h3(3) = (q * S * C_y / (m * g)) * dfalpha_dwx;
    Lg2_Lf_h3(4) = (q * S * C_y / (m * g)) * dfalpha_dwy;
    Lg2_Lf_h3(5) = (q * S * C_y / (m * g)) * dfalpha_dwz;

    // Calculations of Control Law
    Matrix3d Gx;
    Gx.row(0) = Lg1_Lf_h1;
    Gx.row(1) = Lg1_Lf_h2;
    Gx.row(2) = Lg1_Lf_h3;
    
    Vector3d Fx;
    Fx << L2f_h1, L2f_h2, L2f_h3;
    
    MatrixXd K = MatrixXd::Zero(3, 6);
    K(0, 0) = k0_1;
    K(0, 1) = k1_1;
    K(1, 2) = k0_2;
    K(1, 3) = k1_2;
    K(2, 4) = k0_3;
    K(2, 5) = k1_3;
    
    MatrixXd phi_x = MatrixXd::Zero(3, 6);
    phi_x.row(0) << -k1_1 * Lg2_h1 - Lg2_Lf_h1;
    phi_x.row(1) << -k1_2 * Lg2_h2 - Lg2_Lf_h2;
    phi_x.row(2) << -k1_3 * Lg2_h3 - Lg2_Lf_h3;

    double dref_1, dref_2, dref_3;                          //First derivative of the reference inputs
    double ddref_1, ddref_2, ddref_3;                       //Second derivative of the reference inputs

    if (t_meas > t_last)
    {
        dref_1 = Derivative(Ref_Last(0), Ref(0), t_last, t_meas);
        dref_2 = Derivative(Ref_Last(1), Ref(1), t_last, t_meas);
        dref_3 = Derivative(Ref_Last(2), Ref(2), t_last, t_meas);
        ddref_1 = Derivative(dRef_Last(0), dref_1, t_last, t_meas);
        ddref_2 = Derivative(dRef_Last(1), dref_2, t_last, t_meas);
        ddref_3 = Derivative(dRef_Last(2), dref_3, t_last, t_meas);
        
        Ref_Last = Ref;
        dRef_Last << dref_1, dref_2, dref_3;
        t_last = t_meas;
    }
    else
    {
        dref_1 = dRef_Last(0);
        dref_2 = dRef_Last(1);
        dref_3 = dRef_Last(2);
        ddref_1 = ddRef_Last(0);
        ddref_2 = ddRef_Last(1);
        ddref_3 = ddRef_Last(2);
    }

    
    VectorXd Y_(6), Yr_(6), Yr_1(3);
    Y_ << y(0), Lf_h1, y(1), Lf_h2, y(2), Lf_h3;
    Yr_ << Ref(0), dref_1, Ref(1), dref_2, Ref(2), dref_3;
    Yr_1 << ddref_1, ddref_2, ddref_3;
    
    //Control Law
    v = Gx.inverse() * (-Fx - K * (Y_ - Yr_) + Yr_1 + phi_x * Obs_d);
    cout << "Inverse Gx: " << std::endl << Gx.inverse() << std::endl << std::endl;
    cout << "Fx: " << std::endl << Fx << std::endl << std::endl;
    cout << "K: " << std::endl << K << std::endl << std::endl;
    cout << "phi_x" << std::endl << phi_x << std::endl << std::endl;
    cout << "Control Input: " << v.transpose() << std::endl << std::endl;
}

void MPC::Interceptor::DCA()
{
    MatrixXd Wc_0(5, 5), Wc_1(5, 5), Wc_2(5, 5);
    Wc_0 = MatrixXd::Zero(5, 5);
    Wc_0(0, 0) = 10;
    Wc_0(1, 1) = 2;
    Wc_0(2, 2) = 2;
    Wc_0(3, 3) = 5;
    Wc_0(4, 4) = 5;
    Wc_1 = Wc_0;
    Wc_2 = Wc_0;

    Matrix3d Wc_3;
    Wc_3 = Matrix3d::Zero();
    Wc_3(0, 0) = 800;
    Wc_3(1, 1) = 800;
    Wc_3(2, 2) = 800;

    MatrixXd Rc_0(5, 5), Rc_1(5, 5);
    Rc_0 = MatrixXd::Zero(5, 5);
    Rc_1 = MatrixXd::Zero(5, 5);
    int i;

    double u_upper{ 0 }, u_lower{ 0 }, u_max_Ts{ 0 }, u_min_Ts{ 0 };
    for (i = 0; i < 5; i++)
    {
        if (i >= 3)
        {
            u_max_Ts = u_T(i) + u_rate_max(i) * Ts;
            u_upper = Min(u_max(i), u_max_Ts);
            u_min_Ts = u_T(i) + u_rate_min(i) * Ts;
            u_lower = Max(u_min(i), u_min_Ts);
            Rc_0(i, i) = Scale(u(i), u_max(i), u_min(i));
        }
        else
        {
            u_max_Ts = u_T(i) + deg_to_rad(u_rate_max(i)) * Ts;
            u_upper = Min(deg_to_rad(u_max(i)), u_max_Ts);
            u_min_Ts = u_T(i) + deg_to_rad(u_rate_min(i)) * Ts;
            u_lower = Max(deg_to_rad(u_min(i)), u_min_Ts);
            Rc_0(i, i) = Scale(u(i), u_upper, u_lower);
        }
        Rc_1(i, i) = Scale(0, 0.001, -0.001);
    }
    cout << "Rc_0: " << std::endl << Rc_0 << std::endl << std::endl;
    cout << "Rc_1: " << std::endl << Rc_1 << std::endl << std::endl;

    MatrixXd Wc_Temp(5, 5), Gc_Temp(5, 5), I_Temp(5, 5);
    MatrixXd Ec(5, 5), Fc(5, 5), Hc(5, 3), Wc(5, 5), Gc(5, 3);
    
    Wc_Temp = Rc_0.pow(2) * Wc_0.pow(2) + Rc_1.pow(2) * Wc_1.pow(2) + 2 * Wc_2.pow(2) + pow(lambda_c, 2) * Bc.transpose() * Wc_3.pow(2) * Bc;
    Wc = Wc_Temp.sqrt();

    Gc_Temp = Bc * Wc.pow(-2) * Bc.transpose();
    Gc = Wc.pow(-2) * Bc.transpose() * Gc_Temp.inverse();

    I_Temp = (MatrixXd::Identity(5, 5) - Gc * Bc) * Wc.pow(-2);
    Hc = pow(lambda_c, 2) * I_Temp * Bc.transpose() * Wc_3.pow(2) + Gc;
    Fc = I_Temp * Wc_2.pow(2);
    Ec = I_Temp * Rc_1.pow(2) * Wc_1.pow(2);

    //DCA Expression
    u = Ec * u_cs + Fc * u_T + Fc * u_2T + Hc * v;
    cout << "Interceptor Input: " << u.transpose() << std::endl << std::endl;
    cout << "Ec: " << std::endl << Ec << std::endl << std::endl;
    cout << "Fc: " << std::endl << Fc << std::endl << std::endl;
    cout << "Hc: " << std::endl << Hc << std::endl << std::endl;

    u_2T = u_T;
    u_T = u;
}


void MPC::Interceptor::Check_Param() {
}
void MPC::Interceptor::Track(double t_meas) {

}

void MPC::Interceptor::First_Run(objectstate Target_State) {
}
void MPC::Interceptor::Fine_Calculations(objectstate Last_Interceptor_State, objectstate Target_State) {
}

void MPC::Interceptor::Update_Target_State(objectstate Target_State)
{
};

bool MPC::Interceptor::Get_State(objectstate& xo_state) {
	return true;
}
void MPC::Interceptor::Get_Record() {

}

bool MPC::Interceptor::Is_Started() {
	return true;
}
bool MPC::Interceptor::Is_Hit() {
	return true;
}
bool MPC::Interceptor::Is_Missed(int last_phase) {
	return true;
}

bool MPC::Interceptor::Is_Terminated(bool& xo_interceptor, bool& xo_target, int& xo_target_id)
{
	return true;
}



