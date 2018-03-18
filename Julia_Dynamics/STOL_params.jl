using LabelledArrays, StaticArrays


names1 = @SArray [:delta_e,:delta_a,:delta_r,:delta_t,:gravity,:mass,:Jx,:Jy,:Jz,:Jxz,:S_wing,:b,:c,:S_prop,:rho,:k_motor,:k_T_P,:k_Omega,:e,:C_L_0,:C_L_alpha,:C_L_q,:C_L_delta_e,:C_D_0,:C_D_alpha,:C_D_p,:C_D_q,:C_D_delta_e,:C_m_0,:C_m_alpha,:C_m_q,:C_m_delta_e,:C_Y_0,:C_Y_beta,:C_Y_p,:C_Y_r,:C_Y_delta_a,:C_Y_delta_r,:C_ell_0,:C_ell_beta,:C_ell_p,:C_ell_r,:C_ell_delta_a,:C_ell_delta_r,:C_n_0,:C_n_beta,:C_n_p,:C_n_r,:C_n_delta_a,:C_n_delta_r,:C_prop,:M,:epsilon,:alpha0,:gamma0,:gamma1,:gamma2,:gamma3,:gamma4,:gamma5,:gamma6,:gamma7,:gamma8,:C_p_0,:C_p_beta,:C_p_p,:C_p_r,:C_p_delta_a,:C_p_delta_r,:C_r_0,:C_r_beta,:C_r_p,:C_r_r,:C_r_delta_a,:C_r_delta_r,:AR,:wind_n,:wind_e,:wind_d,:L_u,:L_v,:L_w,:sigma_u,:sigma_v,:sigma_w,:Va0,:lambda,:Ts,:pn0,:pe0,:pd0,:u0,:v0,:w0,:phi0,:theta0,:psi0,:p0,:q0,:r0]
values_1 = @MArray [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

P = LMArray(names1,values_1)

# Control inputs
P[:delta_e] = 0.0
P[:delta_a] = 0.0
P[:delta_r] = 0.0
P[:delta_t] = 0.5


P[:gravity] = 9.8;

#################################
# Params for Aersonade UAV
#physical parameters of airframe
P[:mass] = 13.5;
P[:Jx]   = 0.8244;
P[:Jy]   = 1.135;
P[:Jz]   = 1.759;
P[:Jxz]  = .1204;
# aerodynamic coefficients
P[:S_wing]        = 0.55;
P[:b]             = 2.8956;
P[:c]             = 0.18994;
P[:S_prop]        = 0.2027;
P[:rho]           = 1.2682;
P[:k_motor]       = 80;
P[:k_T_P]         = 0;
P[:k_Omega]       = 0;
P[:e]             = 0.9;

P[:C_L_0]         = 0.28;
P[:C_L_alpha]     = 3.45;
P[:C_L_q]         = 0.0;
P[:C_L_delta_e]   = -0.36;
P[:C_D_0]         = 0.03;
P[:C_D_alpha]     = 0.30;
P[:C_D_p]         = 0.0437;
P[:C_D_q]         = 0.0;
P[:C_D_delta_e]   = 0.0;
P[:C_m_0]         = -0.02338;
P[:C_m_alpha]     = -0.38;
P[:C_m_q]         = -3.6;
P[:C_m_delta_e]   = -0.5;
P[:C_Y_0]         = 0.0;
P[:C_Y_beta]      = -0.98;
P[:C_Y_p]         = 0.0;
P[:C_Y_r]         = 0.0;
P[:C_Y_delta_a]   = 0.0;
P[:C_Y_delta_r]   = -0.17;
P[:C_ell_0]       = 0.0;
P[:C_ell_beta]    = -0.12;
P[:C_ell_p]       = -0.26;
P[:C_ell_r]       = 0.14;
P[:C_ell_delta_a] = 0.08;
P[:C_ell_delta_r] = 0.105;
P[:C_n_0]         = 0.0;
P[:C_n_beta]      = 0.25;
P[:C_n_p]         = 0.022;
P[:C_n_r]         = -0.35;
P[:C_n_delta_a]   = 0.06;
P[:C_n_delta_r]   = -0.032;
P[:C_prop]        = 1.0;
P[:M]             = 50;
P[:epsilon]       = 0.1592;
P[:alpha0]        = 0.4712;

##################################

P[:gamma0] = P[:Jx]*P[:Jz] - (P[:Jxz])^2;
P[:gamma1] = (P[:Jxz]*(P[:Jx]-P[:Jy]+P[:Jz]))/P[:gamma0];
P[:gamma2] = (P[:Jz]*(P[:Jz]-P[:Jy])+P[:Jxz]^2)/P[:gamma0];
P[:gamma3] = P[:Jz]/P[:gamma0];
P[:gamma4] = P[:Jxz]/P[:gamma0];
P[:gamma5] = (P[:Jz]-P[:Jx])/P[:Jy];
P[:gamma6] = P[:Jxz]/P[:Jy];
P[:gamma7] = ((P[:Jx]-P[:Jy])*P[:Jx] + P[:Jxz]^2)/P[:gamma0];
P[:gamma8] = P[:Jx]/P[:gamma0];

########

P[:C_p_0] = P[:gamma3]*P[:C_ell_0] + P[:gamma4]*P[:C_n_0];
P[:C_p_beta] = P[:gamma3]*P[:C_ell_beta] + P[:gamma4]*P[:C_n_beta];
P[:C_p_p] = P[:gamma3]*P[:C_ell_p] + P[:gamma4]*P[:C_n_p];
P[:C_p_r] = P[:gamma3]*P[:C_ell_r] + P[:gamma4]*P[:C_n_p];
P[:C_p_delta_a] = P[:gamma3]*P[:C_ell_delta_a] + P[:gamma4]*P[:C_n_delta_a];
P[:C_p_delta_r] = P[:gamma3]*P[:C_ell_delta_r] + P[:gamma4]*P[:C_n_delta_r];

P[:C_r_0] = P[:gamma4]*P[:C_ell_0] + P[:gamma8]*P[:C_n_0];
P[:C_r_beta] = P[:gamma4]*P[:C_ell_beta] + P[:gamma8]*P[:C_n_beta];
P[:C_r_p] = P[:gamma4]*P[:C_ell_p] + P[:gamma8]*P[:C_n_p];
P[:C_r_r] = P[:gamma4]*P[:C_ell_r] + P[:gamma8]*P[:C_n_r];
P[:C_r_delta_a] = P[:gamma4]*P[:C_ell_delta_a] + P[:gamma8]*P[:C_n_delta_a];
P[:C_r_delta_r] = P[:gamma4]*P[:C_ell_delta_r] + P[:gamma8]*P[:C_n_delta_r];

P[:AR] = (P[:b]^2)/P[:S_wing];

#######

# wind parameters
P[:wind_n] = 0;#3;
P[:wind_e] = 0;#2;
P[:wind_d] = 0;
P[:L_u] = 200;
P[:L_v] = 200;
P[:L_w] = 50;
P[:sigma_u] = 1.06;
P[:sigma_v] = 1.06;
P[:sigma_w] = .7;


# compute trim conditions using 'mavsim_chap5_trim.slx'
# initial airspeed
P[:Va0] = 17;
gamma_angle = 0*pi/180;  # desired flight path angle (radians)
R     = Inf;         # desired radius (m) - use (+) for right handed orbit,
P[:lambda] = 100;
# autopilot sample rate
P[:Ts] = 0.01;

# first cut at initial conditions
P[:pn0]    = 0;  # initial North position
P[:pe0]    = 0;  # initial East position
P[:pd0]    = 0;  # initial Down position (negative altitude)
P[:u0]     = P[:Va0]; # initial velocity along body x-axis
P[:v0]     = 0;  # initial velocity along body y-axis
P[:w0]     = 0;  # initial velocity along body z-axis
P[:phi0]   = 0;  # initial roll angle
P[:theta0] = 0;  # initial pitch angle
P[:psi0]   = 0;  # initial yaw angle
P[:p0]     = 0;  # initial body frame roll rate
P[:q0]     = 0;  # initial body frame pitch rate
P[:r0]     = 0;  # initial body frame yaw rate
