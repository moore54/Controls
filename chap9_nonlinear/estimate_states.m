% estimate_states KEVIN
%   - estimate the MAV states using gyros, accels, pressure sensors, and
%   GPS.
%
% Outputs are:
%   pnhat    - estimated North position,
%   pehat    - estimated East position,
%   hhat     - estimated altitude,
%   Vahat    - estimated airspeed,
%   alphahat - estimated angle of attack
%   betahat  - estimated sideslip angle
%   phihat   - estimated roll angle,
%   thetahat - estimated pitch angel,
%   chihat   - estimated course,
%   phat     - estimated roll rate,
%   qhat     - estimated pitch rate,
%   rhat     - estimated yaw rate,
%   Vghat    - estimated ground speed,
%   wnhat    - estimate of North wind,
%   wehat    - estimate of East wind
%   psihat   - estimate of heading angle
%
%
% Modified:  3/15/2010 - RB
%            5/18/2010 - RB
%

function xhat = estimate_states(uu, P)

% rename inputs
y_gyro_x      = uu(1);
y_gyro_y      = uu(2);
y_gyro_z      = uu(3);
y_accel_x     = uu(4);
y_accel_y     = uu(5);
y_accel_z     = uu(6);
y_static_pres = uu(7);
y_diff_pres   = uu(8);
y_gps_n       = uu(9);
y_gps_e       = uu(10);
y_gps_h       = uu(11);
y_gps_Vg      = uu(12);
y_gps_course  = uu(13);
t             = uu(14);


% not estimating these states
alphahat = 0;
betahat  = 0;
bxhat    = 0;
byhat    = 0;
bzhat    = 0;

% Zero for flight
y_static_pres = y_static_pres-P.B_abs_press;


%Declare Persistents
persistent y_gyro_x_old
persistent y_gyro_y_old
persistent y_gyro_z_old
persistent y_accel_x_old
persistent y_accel_y_old
persistent y_accel_z_old
persistent y_static_pres_old
persistent y_diff_pres_old
persistent y_gps_n_old
persistent y_gps_e_old
persistent y_gps_h_old
persistent y_gps_Vg_old
persistent y_gps_course_old
persistent pnhat
persistent pehat
persistent Vghat
persistent chihat

persistent xhat_a
persistent P_a
% persistent phihat
% persistent thetahat

persistent xhat_p
persistent P_p
persistent psihat

wnhat = 0;
wehat = 0;
% chihat = 0;
% Vghat = 0;
% pnhat = 0;
% pehat = 0;


%Initiate Persistents

if t==0

%     y_gyro_x_old = 0.0;
%     y_gyro_y_old = 0.0;
%     y_gyro_z_old = 0.0;
    y_accel_x_old = y_accel_x;
    y_accel_y_old = y_accel_y;
    y_accel_z_old = y_accel_z;
    y_static_pres_old = P.rho*P.gravity*P.pd0;
    y_diff_pres_old =  0.5*P.rho*P.Va0^2;
    y_gps_n_old = 0.0;
    y_gps_e_old = 0.0;
    y_gps_h_old = 0.0;
    y_gps_Vg_old = P.Va0;
    y_gps_course_old = 0.0;

    xhat_a = [P.phi0;P.theta0];
    P_a = diag([(15*pi/180)^2,(15*pi/180)^2]);

    xhat_p = [0;0;P.Va0;0;0;0;0];
    P_p = diag([4^2,4^2,1^2,(10*pi/180)^2,(5*pi/180)^2,5^2,5^2]);
    pnhat = -999;
    pehat = -999;
    Vghat = P.Va0;
    chihat = 0;
    psihat = 0;

end

% Low Pass Filter on
LPFy_gyro_x = y_gyro_x;%LPF(y_gyro_x_old, y_gyro_x, P.y_gyro_x_LPFgain);
LPFy_gyro_y = y_gyro_y;%LPF(y_gyro_y_old, y_gyro_y, P.y_gyro_y_LPFgain);
LPFy_gyro_z = y_gyro_z;%LPF(y_gyro_z_old, y_gyro_z, P.y_gyro_z_LPFgain);
LPFy_accel_x = y_accel_x;%LPF(y_accel_x_old, y_accel_x, P.y_accel_x_LPFgain);
LPFy_accel_y = y_accel_y;%LPF(y_accel_y_old, y_accel_y, P.y_accel_y_LPFgain);
LPFy_accel_z = y_accel_z;%LPF(y_accel_z_old, y_accel_z, P.y_accel_z_LPFgain);
LPFy_static_pres = LPF(y_static_pres_old, y_static_pres, P.y_static_pres_LPFgain);
LPFy_diff_pres = LPF(y_diff_pres_old, y_diff_pres, P.y_diff_pres_LPFgain);
LPFy_gps_n = LPF(y_gps_n_old, y_gps_n, P.y_gps_n_LPFgain);
LPFy_gps_e = LPF(y_gps_e_old, y_gps_e, P.y_gps_e_LPFgain);
LPFy_gps_h = LPF(y_gps_h_old, y_gps_h, P.y_gps_h_LPFgain);
LPFy_gps_Vg = LPF(y_gps_Vg_old, y_gps_Vg, P.y_gps_Vg_LPFgain);
LPFy_gps_course = LPF(y_gps_course_old, y_gps_course, P.y_gps_course_LPFgain);

% Estimate States With Sensor Inversion
phat = LPFy_gyro_x;
qhat = LPFy_gyro_y;
rhat = LPFy_gyro_z;
pnhat = LPFy_gps_n; 
pehat = LPFy_gps_e;
hhat = LPFy_static_pres/(P.rho*P.gravity);
Vahat = sqrt(2/P.rho*LPFy_diff_pres);
%     phihat = atan(LPFy_accel_y/LPFy_accel_z);
%     thetahat = asin(LPFy_accel_x/P.gravity);
% Vghat = LPFy_gps_Vg;
chihat = LPFy_gps_course;
psihat = 0.0;%chihat;

%%%%% Roll and pitch extended Kalman filter %%%%

Q_a = diag([1E-8,1E-8]);
    
% Qs_a = diag([(P.sigma_accel_x)^2, (P.sigma_accel_y)^2, (P.sigma_accel_z)^2]);
R_a = [P.sigma_accel_x^2,P.sigma_accel_y^2,P.sigma_accel_z^2];

y_accel = [LPFy_accel_x,LPFy_accel_y,LPFy_accel_z];
y_accel_old = [y_accel_x_old,y_accel_y_old,y_accel_z_old];

%Estimate roll, pitch with EKF
N = 10;
for i = 1:N
    phihat = xhat_a(1);
    thetahat = xhat_a(2);
    f_a = [phat+qhat*sin(phihat)*tan(thetahat)+rhat*cos(phihat)*tan(thetahat);...
           qhat*cos(phihat)-rhat*sin(phihat)];

    A = [qhat*cos(phihat)*tan(thetahat)-rhat*sin(phihat)*tan(thetahat), (qhat*sin(phihat)+rhat*cos(phihat))/(cos(thetahat))^2;...
        -qhat*sin(phihat)-rhat*cos(phihat), 0];

%     G_a = [...
%         1, sin(phihat)*tan(thetahat), cos(phihat)*tan(thetahat), 0;...
%         0, cos(phihat), -sin(phihat), 0;...
%         ];

    xhat_a = xhat_a+(P.Ts/N)*f_a;
    P_a = P_a + (P.Ts/N) * (A*P_a+P_a*A'+Q_a);
end

phihat = xhat_a(1);
thetahat = xhat_a(2);

h_a = [...
    qhat*Vahat*sin(thetahat)+P.gravity*sin(thetahat);...
    rhat*Vahat*cos(thetahat)-phat*Vahat*sin(thetahat)-P.gravity*cos(thetahat)*sin(phihat);...
    -qhat*Vahat*cos(thetahat)-P.gravity*cos(thetahat)*cos(phihat);...
    ];

dhdx_a = [...
    0, qhat*Vahat*cos(thetahat)+P.gravity*cos(thetahat);...
    -P.gravity*cos(phihat)*cos(thetahat), -rhat*Vahat*sin(thetahat)-phat*Vahat*cos(thetahat)+P.gravity*sin(thetahat)*sin(phihat);...
    P.gravity*sin(phihat)*cos(thetahat), (qhat*Vahat+P.gravity*cos(phihat))*sin(thetahat);...
    ];


 Ri = P.sigma_accel_x^2;
 throwout = 0.1;
    if y_accel_x ~= y_accel_x_old
        % x measurement
        hi = h_a(1);
        Ci = dhdx_a(1,:);
        Li = P_a*Ci'*(Ri + Ci*P_a*Ci')^(-1);
        P_a = (eye(2) - Li*Ci)*P_a;
        if (y_accel_x - hi)/y_accel_x<throwout
            xhat_a = xhat_a + Li*(y_accel_x - hi);
        end
        
        % y measurement
        hi = h_a(2);
        Ci = dhdx_a(2,:);
        Li = P_a*Ci'*(Ri + Ci*P_a*Ci')^(-1);
        P_a = (eye(2) - Li*Ci)*P_a;
        if (y_accel_y - hi)/y_accel_y<throwout
            xhat_a = xhat_a + Li*(y_accel_y - hi);
        end
        
        % z measurement
        hi = h_a(3);
        Ci = dhdx_a(3,:);
        Li = P_a*Ci'*(Ri + Ci*P_a*Ci')^(-1);
        P_a = (eye(2) - Li*Ci)*P_a;
        if (y_accel_z - hi)/y_accel_z<throwout
            xhat_a = xhat_a + Li*(y_accel_z - hi);
        end
        
    end

phihat = xhat_a(1);
thetahat = xhat_a(2);


%%%%%%% EKF for position heading and wind %%%%%%%

Q_p = diag([1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8]);
R_p = [1,1,1,1,1,1,1];

%Estimate wind triangle
y_wn = Vahat*cos(psihat)+wnhat-Vghat*cos(chihat);
y_we = Vahat*sin(psihat)+wehat-Vghat*sin(chihat);
y_p = [LPFy_gps_n,LPFy_gps_e,LPFy_gps_Vg,LPFy_gps_course,y_wn,y_we];


N = 10;
for i = 1:N
    
    pnhat = xhat_p(1);
    pehat = xhat_p(2);
    Vghat = xhat_p(3); 
    chihat = xhat_p(4);
    wnhat = xhat_p(5);
    wehat = xhat_p(6);
    psihat = xhat_p(7);
    
    psidothat = qhat*sin(phihat)/cos(thetahat)+rhat*cos(phihat)/cos(thetahat);
f_p = [...
    Vghat*cos(chihat);...
    Vghat*sin(chihat);...
    Vahat/Vghat*psidothat*(-wnhat*sin(psihat)+wehat*cos(psihat));...
    P.gravity/Vghat*tan(phihat)*cos(chihat-psihat);...
    0;...
    0;...
    qhat*sin(phihat)/cos(thetahat)+rhat*cos(phihat)/cos(thetahat);...
    ];

dVgdotdpsi = -psidothat*Vahat*(wnhat*cos(psihat)+wehat*sin(psihat))/Vghat;
dchidotdVg = -P.gravity/Vghat^2*tan(phihat)*cos(chihat-psihat);
dchidotdchi = -P.gravity/Vghat*tan(phihat)*sin(chihat-psihat);
dchidotdpsi = P.gravity/Vghat*tan(phihat)*sin(chihat-psihat);

Vgdothat = ((Vahat*cos(psihat)+wnhat)*(-Vahat*psidothat*sin(psihat))+(Vahat*sin(psihat)+wehat)*(Vahat*psidothat*cos(psihat)))/Vghat;

dfdx_p = [...
    0,0,cos(chihat),-Vghat*sin(chihat),0,0,0;...
    0,0,sin(chihat),Vghat*cos(chihat),0,0,0;...
    0,0,-Vgdothat/Vghat,0,-psidothat*Vahat*sin(psihat)/Vghat,psidothat*Vahat*cos(psihat)/Vghat,dVgdotdpsi;...
    0,0,dchidotdVg,dchidotdchi,0,0,dchidotdpsi;...
    0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0;...
    ];


    xhat_p = xhat_p+(P.Ts/N)*f_p;
    A = dfdx_p;
    P_p = P_p + (P.Ts/N) * (A*P_p+P_p*A'+Q_p);
end

pnhat = xhat_p(1);
pehat = xhat_p(2);
Vghat = xhat_p(3); 
chihat = xhat_p(4);
wnhat = xhat_p(5);
wehat = xhat_p(6);
psihat = xhat_p(7);

h_p = [...
    pnhat;...
    pehat;...
    Vghat;...
    chihat;...
    Vahat*cos(psihat)+wnhat-Vghat*cos(chihat);...
    Vahat*sin(psihat)+wehat-Vghat*sin(chihat);...
    ];



dhdx_p = [...
    1,0,0,0,0,0,0;...
    0,1,0,0,0,0,0;...
    0,0,1,0,0,0,0;...
    0,0,0,1,0,0,0;...
    0,0,-cos(chihat),Vghat*sin(chihat),1,0,-Vahat*sin(phihat);...
    0,0,-sin(chihat),-Vghat*cos(chihat),0,1,Vahat*cos(phihat);...
    ];

%Update if new value from sensors
 throwout = 0.1;
    if y_accel_x ~= y_accel_x_old
        % LPFy_gps_n measurement
        Ri = P.sigma_GPS_n^2;
        hi = h_p(1);
        Ci = dhdx_p(1,:);
        Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
        P_p = (eye(7) - Li*Ci)*P_p;
        if (LPFy_gps_n - hi)/hi<throwout
            xhat_p = xhat_p + Li*(LPFy_gps_n - hi);
        end
        
        % LPFy_gps_e measurement
        Ri = P.sigma_GPS_e^2;
        hi = h_p(2);
        Ci = dhdx_p(2,:);
        Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
        P_p = (eye(7) - Li*Ci)*P_p;
        if (LPFy_gps_e - hi)/LPFy_gps_e<throwout
            xhat_p = xhat_p + Li*(LPFy_gps_e - hi);
        end
        
        % LPFy_gps_Vg measurement
        Ri = P.sigma_GPS_Vg^2;
        hi = h_p(3);
        Ci = dhdx_p(3,:);
        Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
        P_p = (eye(7) - Li*Ci)*P_p;
        if (LPFy_gps_Vg - hi)/LPFy_gps_Vg<throwout
            xhat_p = xhat_p + Li*(LPFy_gps_Vg - hi);
        end
        
        % LPFy_gps_course
        Ri = P.sigma_GPS_course^2;
        hi = h_p(4);
        Ci = dhdx_p(4,:);
        Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
        P_p = (eye(7) - Li*Ci)*P_p;
        if (LPFy_gps_course - hi)/LPFy_gps_course<throwout
            xhat_p = xhat_p + Li*(LPFy_gps_course - hi);
        end
        
        % y_wn
        Ri = P.sigma_u^2;
        hi = h_p(5);
        Ci = dhdx_p(5,:);
        Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
        P_p = (eye(7) - Li*Ci)*P_p;
        if (y_wn - hi)/y_wn<throwout
            xhat_p = xhat_p + Li*(y_wn - hi);
        end
        
        % y_we
        Ri = P.sigma_v^2;
        hi = h_p(6);
        Ci = dhdx_p(6,:);
        Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
        P_p = (eye(7) - Li*Ci)*P_p;
        if (y_we - hi)/y_we<throwout
            xhat_p = xhat_p + Li*(y_we - hi);
        end
        
%         % y_we
%         Ri = P.sigma_v^2;
%         hi = h_p(7);
%         Ci = dhdx_p(7,:);
%         Li = P_p*Ci'*(Ri + Ci*P_p*Ci')^(-1);
%         P_p = (eye(7) - Li*Ci)*P_p;
%         xhat_p = xhat_p + Li*(y_we - hi);
%         
    end


pnhat = xhat_p(1);
pehat = xhat_p(2);
Vghat = xhat_p(3); 
chihat = xhat_p(4);
wnhat = xhat_p(5);
wehat = xhat_p(6);
psihat = xhat_p(7);

% if abs(Vghat)<1E-6
%     Vghat = 1E-6;
% %     fprintf('Vg limited to 1E-6')
% end

% while (chihat>pi)
%     chihat = chihat-2*pi;
% end


%Update
y_gyro_x_old = LPFy_gyro_x;
y_gyro_y_old = LPFy_gyro_y;
y_gyro_z_old = LPFy_gyro_z;
y_accel_x_old = LPFy_accel_x;
y_accel_y_old = LPFy_accel_y;
y_accel_z_old = LPFy_accel_z;
y_static_pres_old = LPFy_static_pres;
y_diff_pres_old = LPFy_diff_pres;
y_gps_n_old = LPFy_gps_n;
y_gps_e_old = LPFy_gps_e;
y_gps_h_old = LPFy_gps_h;
y_gps_Vg_old = LPFy_gps_Vg;
y_gps_course_old = LPFy_gps_course;

%Output vars
xhat = [...
    pnhat;...
    pehat;...
    hhat;...
    Vahat;...
    alphahat;...
    betahat;...
    phihat;...
    thetahat;...
    chihat;...
    phat;...
    qhat;...
    rhat;...
    Vghat;...
    wnhat;...
    wehat;...
    psihat;...
    bxhat;...
    byhat;...
    bzhat;...
    ];
end

function out = LPF(old, current, gain)
out = gain*old+(1-gain)*current;
end


