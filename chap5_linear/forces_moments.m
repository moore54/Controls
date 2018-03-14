% forces_moments.m
%   Computes the forces and moments acting on the airframe.
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack: angles in radians
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%

function out = forces_moments(x, delta, wind, P)

    % relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis
    w_wg    = wind(6); % gust along body z-axis

    R_roll = [...
            1, 0, 0;...
            0, cos(phi), sin(phi);...
            0, -sin(phi), cos(phi)];
    R_pitch = [...
            cos(theta), 0, -sin(theta);...
            0, 1, 0;...
            sin(theta), 0, cos(theta)];
    R_yaw = [...
            cos(psi), sin(psi), 0;...
            -sin(psi), cos(psi), 0;...
            0, 0, 1];

    R = R_roll*R_pitch*R_yaw;
      % note that R above either leaves the vector alone or rotates
      % a vector in a left handed rotation.  We want to rotate all
      % points in a right handed rotation, so we must transpose
    R = R';%'

    % compute wind data in NED $TODO: include gusts
    W_gustBODY = [u_wg;v_wg;w_wg];
%     W_gustNED = W_gustBODY\R
    w_n = w_ns;%+W_gustNED(1);
    w_e = w_es;%+W_gustNED(2);
    w_d = w_ds;%+W_gustNED(3);

    W_BODY = R*[w_n;w_e;w_d]+[W_gustBODY(1);W_gustBODY(2);W_gustBODY(3)];

    % compute air data in body frame
    ur = u-W_BODY(1);
    vr = v-W_BODY(2);
    wr = w-W_BODY(3);
    
    % compute total wind in ground frame
    W_vec = R'*W_BODY;
    w_n = W_vec(1);
    w_e = W_vec(2);
    w_d = W_vec(3);

    Va = sqrt(ur^2+vr^2+wr^2);
    alpha = atan2(wr,ur);%-P.alpha0;
    beta = asin(vr/Va);
    
    if (Va == 0)
        Va = P.Va0;
    end
    if ~isfinite(alpha)
        alpha = 0;
    end
    if ~isfinite(beta)
        beta = 0;
    end

   
   C_D_alpha = P.C_D_p + ((P.C_L_0 + P.C_L_alpha*alpha)^2)/(pi*P.e*P.AR);
   sigma_alpha = (1 + exp(-P.M*(alpha-P.alpha0)) + exp(P.M*(alpha+P.alpha0)))/((1+exp(-P.M*(alpha-P.alpha0)))*(1+exp(P.M*(alpha+P.alpha0))));
   C_L_alpha = (1-sigma_alpha)*(P.C_L_0 + P.C_L_alpha*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));
    
   C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
   C_X_q_alpha = -P.C_D_q*cos(alpha) + P.C_L_q*sin(alpha);
   C_X_delta_e_alpha = -P.C_D_delta_e*cos(alpha) + P.C_L_delta_e*sin(alpha);
   C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);
   C_Z_q_alpha = -P.C_D_q*sin(alpha) - P.C_L_q*cos(alpha);
   C_Z_delta_e_alpha = -P.C_D_delta_e*sin(alpha) - P.C_L_delta_e*cos(alpha);
   
    % compute external forces and torques on aircraft
    % Gravity force is in the X DIRECTION !!!
    grav_force_x = P.mass*P.gravity*(cos(theta)*cos(psi));
    grav_force_y = P.mass*P.gravity*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi));
    grav_force_z = P.mass*P.gravity*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
    
    Force(1) =  grav_force_x+0.5*P.rho*Va^2*P.S_wing*(C_X_alpha+C_X_q_alpha*P.c/(2*Va)*q+C_X_delta_e_alpha*delta_e)...
        +0.5*P.rho*P.S_prop*P.C_prop*((P.k_motor*delta_t)^2-Va^2);
    Force(2) = grav_force_y+0.5*P.rho*Va^2*P.S_wing*(P.C_Y_0+P.C_Y_beta*beta+P.C_Y_p*P.b/(2*Va)*p+...
        P.C_Y_r*P.b/(2*Va)*r+P.C_Y_delta_a*delta_a+P.C_Y_delta_r*delta_r);
    Force(3) = grav_force_z+0.5*P.rho*Va^2*P.S_wing*(C_Z_alpha+C_Z_q_alpha*P.c/(2*Va)*q+C_Z_delta_e_alpha*delta_e);

    Torque(1) = 0.5*P.rho*Va^2*P.S_wing*P.b*(P.C_ell_0 + P.C_ell_beta*beta + (P.C_ell_p*P.b*p/(2*Va)) + (P.C_ell_r*P.b*r/(2*Va)) + P.C_ell_delta_a*delta_a + P.C_ell_delta_r*delta_r)...
        - P.k_T_P*(P.k_Omega*delta_t)^2;
    Torque(2) = 0.5*P.rho*Va^2*P.S_wing*P.c*(P.C_m_0 + P.C_m_alpha*alpha + (P.C_m_q*P.c*q/(2*Va)) + P.C_m_delta_e*delta_e);   
    Torque(3) = 0.5*P.rho*Va^2*P.S_wing*P.b*(P.C_n_0 + P.C_n_beta*beta + (P.C_n_p*P.b*p/(2*Va)) + (P.C_n_r*P.b*r/(2*Va)) + P.C_n_delta_a*delta_a + P.C_n_delta_r*delta_r);
   
    out = [Force'; Torque'; Va; alpha; beta; w_n; w_e; w_d];
end
