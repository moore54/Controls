% estimate_states
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
    
    % Temporary
    chihat = 0;
    Vghat = 0;
    wnhat = 0;
    wehat = 0;
    
% Put my state estimation code here
    persistent pnhat
    persistent pehat
    persistent hhat
    persistent LPF_y_diff_pres
    persistent phat
    persistent qhat
    persistent rhat
%     persistent LPF_y_accel_x
%     persistent LPF_y_accel_y
%     persistent LPF_y_accel_z
    persistent x_hat
    persistent P_EKF
    persistent y_accel_x_old
    persistent y_accel_y_old
    persistent y_accel_z_old

    
    a = 15;
    a_h = 100;
    a_Va = 50;
    a_accel = 200;
    if t == 0
        pnhat = 0;
        pehat = 0;
        hhat = -P.pd0;
        LPF_y_diff_pres = 0.5*P.rho*P.Va0^2;
%         LPF_y_accel_x = 0;
%         LPF_y_accel_y = 0;
%         LPF_y_accel_z = -P.gravity;
        phat = 0;
        qhat = 0;
        rhat = 0;
        x_hat = [P.phi0; P.theta0];
        P_EKF = diag([(15*pi/180)^2, (15*pi/180)^2]);
        y_accel_x_old = y_accel_x;
        y_accel_y_old = y_accel_y;
        y_accel_z_old = y_accel_z;
    end
    
    pnhat = exp(-a*P.Ts)*pnhat + (1 - exp(-a*P.Ts))*y_gps_n;
    pehat = exp(-a*P.Ts)*pehat + (1 - exp(-a*P.Ts))*y_gps_e;
    hhat = exp(-a_h*P.Ts)*hhat + (1 - exp(-a_h*P.Ts))*y_static_pres/(P.rho*P.gravity);
    LPF_y_diff_pres = exp(-a_Va*P.Ts)*LPF_y_diff_pres + (1 - exp(-a_Va*P.Ts))*y_diff_pres;
    Vahat = sqrt(2/P.rho*LPF_y_diff_pres);
%     LPF_y_accel_x = exp(-a_accel*P.Ts)*LPF_y_accel_x + (1 - exp(-a_accel*P.Ts))*y_accel_x;
%     LPF_y_accel_y = exp(-a_accel*P.Ts)*LPF_y_accel_y + (1 - exp(-a_accel*P.Ts))*y_accel_y;
%     LPF_y_accel_z = exp(-a_accel*P.Ts)*LPF_y_accel_z + (1 - exp(-a_accel*P.Ts))*y_accel_z;
%     phihat = atan(LPF_y_accel_y/LPF_y_accel_z);
%     thetahat = asin(LPF_y_accel_x/P.gravity);
    psihat = 0;
%     phat = exp(-a*P.Ts)*phat + (1 - exp(-a*P.Ts))*y_gyro_x;
%     qhat = exp(-a*P.Ts)*qhat + (1 - exp(-a*P.Ts))*y_gyro_y;
%     rhat = exp(-a*P.Ts)*rhat + (1 - exp(-a*P.Ts))*y_gyro_z;
    phat = y_gyro_x;
    qhat = y_gyro_y;
    rhat = y_gyro_z;
    
    %% EKF Stuff
    % phihat & thetahat
    % Predict
    N = 10;
    for i = 1:N
        phihat = x_hat(1);
        thetahat = x_hat(2);
        cp = cos(phihat);
        sp = sin(phihat);
        tp = tan(phihat);
        ct = cos(thetahat);
        st = sin(thetahat);
        tt = tan(thetahat);
        f_xhat_u = [phat + qhat*sp*tt + rhat*cp*tt;
                    qhat*cp - rhat*sp];
        
        A = [qhat*cp*tt - rhat*sp*tt, (qhat*sp + rhat*cp)/ct^2;
            -qhat*sp - rhat*cp,                  0];
        
        Q = diag([1e-8, 1e-8]);
        x_hat = x_hat + (P.Ts/N)*f_xhat_u;
        P_EKF = P_EKF + (P.Ts/N)*(A*P_EKF + P_EKF*A' + Q);
    end
    phihat = x_hat(1);
    thetahat = x_hat(2);
    cp = cos(phihat);
        sp = sin(phihat);
        ct = cos(thetahat);
        st = sin(thetahat);
    % Measurement update
    h = [qhat*Vahat*st + P.gravity*st;
         rhat*Vahat*ct - phat*Vahat*st - P.gravity*ct*sp;
         -qhat*Vahat*ct - P.gravity*ct*cp];
           
    d_h = [0, qhat*Vahat*ct + P.gravity*ct;
           -P.gravity*cp*ct,  -rhat*Vahat*st - phat*Vahat*ct + P.gravity*sp*st;
           P.gravity*sp*ct, (qhat*Vahat + P.gravity*cp)*st];
       
    Ri = P.sigma_accel_x^2;
    if y_accel_x ~= y_accel_x_old
        % x measurement
        hi = h(1);
        Ci = d_h(1,:);
        Li = P_EKF*Ci'*(Ri + Ci*P_EKF*Ci')^(-1);
        P_EKF = (eye(2) - Li*Ci)*P_EKF;
        x_hat = x_hat + Li*(y_accel_x - hi);
        
        % y measurement
        hi = h(2);
        Ci = d_h(2,:);
        Li = P_EKF*Ci'*(Ri + Ci*P_EKF*Ci')^(-1);
        P_EKF = (eye(2) - Li*Ci)*P_EKF;
        x_hat = x_hat + Li*(y_accel_y - hi);
        
        % z measurement
        hi = h(3);
        Ci = d_h(3,:);
        Li = P_EKF*Ci'*(Ri + Ci*P_EKF*Ci')^(-1);
        P_EKF = (eye(2) - Li*Ci)*P_EKF;
        x_hat = x_hat + Li*(y_accel_z - hi);
        
    end
           
    phihat = x_hat(1);
    thetahat = x_hat(2);
    
%     %% pn, pe, Vg, chi, wn, we, psi
%     % Predict
%     N = 10;
%     for i = 1:N
%         pnhat = x_hat_2(1);
%         pehat = x_hat_2(2);
%         Vghat = x_hat_2(2);
%         chihat = x_hat_2(2);
%         wnhat = x_hat_2(2);
%         wehat = x_hat_2(2);
%         psihat = x_hat_2(2);
%         cp = cos(phihat);
%         sp = sin(phihat);
%         tp = tan(phihat);
%         ct = cos(thetahat);
%         st = sin(thetahat);
%         tt = tan(thetahat);
%         f_xhat_u = [phat + qhat*sp*tt + rhat*cp*tt;
%             qhat*cp - qhat*sp];
%         
%         A = [qhat*cp*tt - qhat*sp*tt, (qhat*sp + qhat*cp)/ct^2;
%             -qhat*sp - qhat*cp,                  0];
%         
%         Q = diag([1e-6, 1e-6]);
%         x_hat = x_hat + (P.Ts/N)*f_xhat_u;
%         P_EKF = P_EKF + (P.Ts/N)*(A*P_EKF + P_EKF*A' + Q);
%     end
%     
%     % Measurement update
%     h = [qhat*Vahat*st + P.gravity*st;
%                 rhat*Vahat*ct - phat*Vahat*st - P.gravity*ct*sp;
%                 -qhat*Vahat*ct - P.gravity*ct*cp];
%            
%     d_h = [0, qhat*Vahat*ct + P.gravity*ct;
%            -P.gravity*cp*ct,  -rhat*Vahat*st - phat*Vahat*ct + P.gravity*sp*st;
%            P.gravity*sp*ct, (qhat*Vahat + P.gravity*cp)*st];
%        
%     Ri = P.sigma_accel_x^2;
%     if y_accel_x ~= y_accel_x_old
%         % x measurement
%         hi = h(1);
%         Ci = d_h(1,:);
%         Li = P_EKF*Ci'*(Ri + Ci*P_EKF*Ci')^(-1);
%         P_EKF = (eye(2) - Li*Ci)*P_EKF;
%         x_hat = x_hat + Li*(y_accel_x - hi);
%         
%         % y measurement
%         hi = h(2);
%         Ci = d_h(2,:);
%         Li = P_EKF*Ci'*(Ri + Ci*P_EKF*Ci')^(-1);
%         P_EKF = (eye(2) - Li*Ci)*P_EKF;
%         x_hat = x_hat + Li*(y_accel_y - hi);
%         
%         % z measurement
%         hi = h(3);
%         Ci = d_h(3,:);
%         Li = P_EKF*Ci'*(Ri + Ci*P_EKF*Ci')^(-1);
%         P_EKF = (eye(2) - Li*Ci)*P_EKF;
%         x_hat = x_hat + Li*(y_accel_z - hi);
%         
%     end
%            
%     phihat = x_hat(1);
%     thetahat = x_hat(2);
    

%% Output Definitions    
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