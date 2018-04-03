P.gravity = 9.8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params for Aersonade UAV
%physical parameters of airframe
P.mass = 13.5;
P.Jx   = 0.8244;
P.Jy   = 1.135;
P.Jz   = 1.759;
P.Jxz  = .1204;
% aerodynamic coefficients
P.S_wing        = 0.55;
P.b             = 2.8956;
P.c             = 0.18994;
P.S_prop        = 0.2027;
P.rho           = 1.2682;
P.k_motor       = 80;
P.k_T_P         = 0;
P.k_Omega       = 0;
P.e             = 0.9;

P.C_L_0         = 0.28;
P.C_L_alpha     = 3.45;
P.C_L_q         = 0.0;
P.C_L_delta_e   = -0.36;
P.C_D_0         = 0.03;
P.C_D_alpha     = 0.30;
P.C_D_p         = 0.0437;
P.C_D_q         = 0.0;
P.C_D_delta_e   = 0.0;
P.C_m_0         = -0.02338;
P.C_m_alpha     = -0.38;
P.C_m_q         = -3.6;
P.C_m_delta_e   = -0.5;
P.C_Y_0         = 0.0;
P.C_Y_beta      = -0.98;
P.C_Y_p         = 0.0;
P.C_Y_r         = 0.0;
P.C_Y_delta_a   = 0.0;
P.C_Y_delta_r   = -0.17;
P.C_ell_0       = 0.0;
P.C_ell_beta    = -0.12;
P.C_ell_p       = -0.26;
P.C_ell_r       = 0.14;
P.C_ell_delta_a = 0.08;
P.C_ell_delta_r = 0.105;
P.C_n_0         = 0.0;
P.C_n_beta      = 0.25;
P.C_n_p         = 0.022;
P.C_n_r         = -0.35;
P.C_n_delta_a   = 0.06;
P.C_n_delta_r   = -0.032;
P.C_prop        = 1.0;
P.M             = 50;
P.epsilon       = 0.1592;
P.alpha0        = 0.4712;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.gamma = P.Jx*P.Jz - (P.Jxz)^2;
P.gamma1 = (P.Jxz*(P.Jx-P.Jy+P.Jz))/P.gamma;
P.gamma2 = (P.Jz*(P.Jz-P.Jy)+P.Jxz^2)/P.gamma;
P.gamma3 = P.Jz/P.gamma;
P.gamma4 = P.Jxz/P.gamma;
P.gamma5 = (P.Jz-P.Jx)/P.Jy;
P.gamma6 = P.Jxz/P.Jy;
P.gamma7 = ((P.Jx-P.Jy)*P.Jx + P.Jxz^2)/P.gamma;
P.gamma8 = P.Jx/P.gamma;

%%%%%%%%

P.C_p_0 = P.gamma3*P.C_ell_0 + P.gamma4*P.C_n_0;
P.C_p_beta = P.gamma3*P.C_ell_beta + P.gamma4*P.C_n_beta;
P.C_p_p = P.gamma3*P.C_ell_p + P.gamma4*P.C_n_p;
P.C_p_r = P.gamma3*P.C_ell_r + P.gamma4*P.C_n_p;
P.C_p_delta_a = P.gamma3*P.C_ell_delta_a + P.gamma4*P.C_n_delta_a;
P.C_p_delta_r = P.gamma3*P.C_ell_delta_r + P.gamma4*P.C_n_delta_r;

P.C_r_0 = P.gamma4*P.C_ell_0 + P.gamma8*P.C_n_0;
P.C_r_beta = P.gamma4*P.C_ell_beta + P.gamma8*P.C_n_beta;
P.C_r_p = P.gamma4*P.C_ell_p + P.gamma8*P.C_n_p;
P.C_r_r = P.gamma4*P.C_ell_r + P.gamma8*P.C_n_r;
P.C_r_delta_a = P.gamma4*P.C_ell_delta_a + P.gamma8*P.C_n_delta_a;
P.C_r_delta_r = P.gamma4*P.C_ell_delta_r + P.gamma8*P.C_n_delta_r;

P.AR = (P.b^2)/P.S_wing;

%%%%%%%

% wind parameters
P.wind_n = 3;%0;%
P.wind_e = -3;%0;%
P.wind_d = 0;
P.L_u = 200;
P.L_v = 200;
P.L_w = 50;
P.sigma_u = 1.06;
P.sigma_v = 1.06;
P.sigma_w = .7;


% compute trim conditions using 'mavsim_chap5_trim.slx'
% initial airspeed
P.Va0 = 35;
gamma = 0*pi/180;  % desired flight path angle (radians)
R     = Inf;         % desired radius (m) - use (+) for right handed orbit,
P.lambda = 100;
% autopilot sample rate
P.Ts = 0.01;

% first cut at initial conditions
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -100;  % initial Down position (negative altitude)
P.u0     = P.Va0; % initial velocity along body x-axis
P.v0     = 0;  % initial velocity along body y-axis
P.w0     = 0;  % initial velocity along body z-axis
P.phi0   = 0;  % initial roll angle
P.theta0 = 0;  % initial pitch angle
P.psi0   = 0;  % initial yaw angle
P.p0     = 0;  % initial body frame roll rate
P.q0     = 0;  % initial body frame pitch rate
P.r0     = 0;  % initial body frame yaw rate

                    %                          (-) for left handed orbit

% run trim commands
[x_trim, u_trim]=compute_trim('mavsim_trim',P.Va0,gamma,R);
P.u_trim = u_trim;
P.x_trim = x_trim;

% set initial conditions to trim conditions
% initial conditions
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -100;  % initial Down position (negative altitude)
P.u0     = x_trim(4);  % initial velocity along body x-axis
P.v0     = x_trim(5);  % initial velocity along body y-axis
P.w0     = x_trim(6);  % initial velocity along body z-axis
P.phi0   = x_trim(7);  % initial roll angle
P.theta0 = x_trim(8);  % initial pitch angle
P.psi0   = x_trim(9);  % initial yaw angle
P.p0     = x_trim(10);  % initial body frame roll rate
P.q0     = x_trim(11);  % initial body frame pitch rate
P.r0     = x_trim(12);  % initial body frame yaw rate

% compute different transfer functions
[T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r,T]...
    = compute_tf_model(x_trim,u_trim,P);

% linearize the equations of motion around trim conditions
[A_lon, B_lon, A_lat, B_lat] = compute_ss_model('mavsim_trim',x_trim,u_trim);
out = eig(A_lon);

P.tau = 0.05;

%%%% ROLL (phi) %%%%

P.delta_a_max = 50*pi/180;
P.e_phi_max = 20*pi/180;
P.zeta_phi = 1.107;
P.wn_phi = sqrt(abs(T.a_phi2)*P.delta_a_max/P.e_phi_max);

P.kp_phi = P.delta_a_max/P.e_phi_max*sign(T.a_phi2);
P.ki_phi = 0.00001;
P.kd_phi = (2*P.wn_phi*P.zeta_phi-T.a_phi1)/T.a_phi2;

%%%% COURSE (chi) %%%%

% P.delta_a_max = 20*pi/180;
% P.e_phi_max = 30*pi/180;
P.zeta_chi = .607;
P.W_chi = 8.0;
P.wn_chi = P.wn_phi/P.W_chi;
Vg = P.Va0;

P.kp_chi = 2*P.zeta_chi*P.wn_chi*Vg/P.gravity;
P.ki_chi = P.wn_chi^2*Vg/P.gravity;
P.kd_chi = 0.005;

%%%% Pitch (theta) outputs delta_e %%%%
P.theta_max = 45.0*pi/180;
P.delta_e_max = 20*pi/180;
P.e_theta_max = 10*pi/180;
P.zeta_theta = 1.107;

P.kp_theta = P.delta_e_max/P.e_theta_max*sign(T.a_theta3);

P.wn_theta = sqrt(T.a_theta2+P.kp_theta*T.a_theta3);

P.ki_theta = 0.0001;
P.kd_theta = (2*P.wn_theta*P.zeta_theta-T.a_theta1)/T.a_theta3;
P.K_thetaDC = P.kp_theta*T.a_theta3/(T.a_theta2+P.kp_theta*T.a_theta3);

%%%% airspeed_with_throttle_hold (V) %%%%

P.delta_t_max = 1;
P.V_max = 50;
P.e_V_max = 1;
P.zeta_V = .707;
P.W_V = 2.0;
P.wn_V = P.wn_theta/P.W_V;

P.kp_V = (2*P.wn_V*P.zeta_V-T.a_V1)/T.a_V2;
P.ki_V = P.wn_V^2/T.a_V2;
P.kd_V = 0.000001;

%%%% airspeed_with_pitch_hold (V) %%%%

P.zeta_V2 = .707;
P.W_V2 = 10.0;
P.wn_V2 = P.wn_theta/P.W_V2;

P.kp_V2 = -(2*P.wn_V2*P.zeta_V2-T.a_V1)/(P.K_thetaDC*P.gravity);
P.ki_V2 = -P.wn_V2^2/(P.K_thetaDC*P.gravity);
P.kd_V2 = 0.000000001;

%%%% Altitude_hold outputs pitch_c %%%%

P.zeta_h = .707;
P.W_h = 10.0;
P.wn_h = P.wn_theta/P.W_h;

P.kp_h = (2*P.wn_h*P.zeta_h)/(P.K_thetaDC*P.Va0);
P.ki_h = P.wn_h^2/(P.K_thetaDC*P.Va0);
P.kd_h = 0.000001;

%%%% ALTUTUDE ZONES %%%
P.altitude_take_off_zone = 50; %meters
P.altitude_hold_zone = 200; %meters


%%%% SENSOR DATA %%%%
P.Ts_gps = 0.1;
P.k_GPS = 1/1100;
P.press0 = 0.0;

P.bias_gyro_x = 0.0;
P.bias_gyro_y = 0.0;
P.bias_gyro_z = 0.0;
P.sigma_gyro_x = 0.13*pi/180;
P.sigma_gyro_y = 0.13*pi/180;
P.sigma_gyro_z = 0.13*pi/180;

P.sigma_accel_x = 0.0025*9.81;
P.sigma_accel_y = 0.0025*9.81;
P.sigma_accel_z = 0.0025*9.81;

P.B_diff_press = 0.02*1000;
P.sigma_diff_press = 0.002*1000;

P.B_abs_press = 0.125*1000;
P.sigma_abs_press = 0.01*1000;

P.B_GPS_n = 4.7;
P.B_GPS_e = 4.7;
P.B_GPS_h = 9.2;
P.sigma_GPS_n = 0.4;
P.sigma_GPS_e = 0.4;
P.sigma_GPS_h = 0.7;
P.sigma_GPS_Vg = 1.0;
P.sigma_GPS_course = P.sigma_GPS_Vg/P.Va0;

P.Error_Factor = 1;
%%%% LPF GAINS %%%%
Gain_macro_tune = 1;
P.y_gyro_x_LPFgain = exp(-50.0*Gain_macro_tune*P.Ts);      
P.y_gyro_y_LPFgain = exp(-50.0*Gain_macro_tune*P.Ts);      
P.y_gyro_z_LPFgain = exp(-50.0*Gain_macro_tune*P.Ts);      
P.y_accel_x_LPFgain = exp(-1E7*Gain_macro_tune*P.Ts);     
P.y_accel_y_LPFgain = exp(-1E7*Gain_macro_tune*P.Ts);     
P.y_accel_z_LPFgain = exp(-1E7*Gain_macro_tune*P.Ts);     
P.y_static_pres_LPFgain = exp(-5.0*Gain_macro_tune*P.Ts); 
P.y_diff_pres_LPFgain = exp(-1000.0*Gain_macro_tune*P.Ts);   
P.y_gps_n_LPFgain = exp(-15*Gain_macro_tune*P.Ts);       
P.y_gps_e_LPFgain = exp(-15*Gain_macro_tune*P.Ts);       
P.y_gps_h_LPFgain = exp(-15*Gain_macro_tune*P.Ts);       
P.y_gps_Vg_LPFgain = exp(-1E7*Gain_macro_tune*P.Ts);      
P.y_gps_course_LPFgain = exp(-1E7*Gain_macro_tune*P.Ts);  

%%%% Simplified Guidance Model %%%%
P.b_chidot = 20.0;
P.b_chi = 5.0;
P.b_hdot = 2.0;
P.b_h = 1.0;
P.b_Va = 4.0;
P.b_phi = 1.0;
P.gamma_max = 45*pi/180;
%%%% Orbit and Straight Line Following %%%%
P.korbit = 10.0;
P.kpath = 0.005;
P.chi_inf = 89.0*pi/180;

% need to add the following to your parameter file:
 
% chapter 11 - path manager
% number of waypoints in data structure
P.size_waypoint_array = 500;
P.phi_max = 60.0*pi/180;
P.R_min = P.Va0^2/P.gravity/tan(P.phi_max)*.5;

% create random city map
city_width      = 2000;  % the city is of size (width)x(width)
building_height = 300;   % maximum height of buildings
%building_height = 1;   % maximum height of buildings (for camera)
num_blocks      = 5;    % number of blocks in city
street_width    = .8;   % percent of block that is street.
% P.pd0           = -P;  % initial height of MAV
P.map = createWorld(city_width, building_height, num_blocks, street_width);

