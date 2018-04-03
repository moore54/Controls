function y = autopilot(uu,P)
%
% autopilot for mavsim
%
% Modification History:
%   2/11/2010 - RWB
%   5/14/2010 - RWB
%   9/30/2014 - RWB
%

% process inputs
NN = 0;
%    pn       = uu(1+NN);  % inertial North position
%    pe       = uu(2+NN);  % inertial East position
h        = uu(3+NN);  % altitude
Va       = uu(4+NN);  % airspeed
%    alpha    = uu(5+NN);  % angle of attack
%    beta     = uu(6+NN);  % side slip angle
phi      = uu(7+NN);  % roll angle
theta    = uu(8+NN);  % pitch angle
chi      = uu(9+NN);  % course angle
p        = uu(10+NN); % body frame roll rate
q        = uu(11+NN); % body frame pitch rate
r        = uu(12+NN); % body frame yaw rate
%    Vg       = uu(13+NN); % ground speed
%    wn       = uu(14+NN); % wind North
%    we       = uu(15+NN); % wind East
%    psi      = uu(16+NN); % heading
%    bx       = uu(17+NN); % x-gyro bias
%    by       = uu(18+NN); % y-gyro bias
%    bz       = uu(19+NN); % z-gyro bias
NN = NN+19;
Va_c     = uu(1+NN);  % commanded airspeed (m/s)
h_c      = uu(2+NN);  % commanded altitude (m)
chi_c    = uu(3+NN);  % commanded course (rad)
NN = NN+3;
t        = uu(1+NN);   % time

autopilot_version = 4;
TOL_Error = 0.10;
% autopilot_version == 1 <- used for tuning
% autopilot_version == 2 <- standard autopilot defined in book
% autopilot_version == 3 <- Total Energy Control for longitudinal AP
switch autopilot_version
    case 1
        [delta, x_command] = autopilot_tuning(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error);
    case 2
        [delta, x_command] = autopilot_uavbook(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error);
    case 3
        [delta, x_command] = autopilot_TECS(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error);
    case 4
        [delta, x_command] = autopilot_kevin(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error);
end
y = [delta; x_command];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autopilot versions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = roll_hold(y_c, y, ~, P,myflag,TOL_Error)

kp=P.kp_phi;
ki=P.ki_phi;
kd=P.kd_phi;
limit=P.delta_a_max;
Ts=P.Ts;
tau=P.tau;

persistent integrator_roll;
persistent differentiator_roll;
persistent error_d1_roll;

if myflag==1 % reset (initialize) persistent variables
    % when myflag==1
    integrator_roll = 0;
    differentiator_roll = 0;
    error_d1_roll = 0; % _d1 means delayed by one time step
end
error = y_c - y; % compute the current error
integrator_roll = integrator_roll + (Ts/2)*(error + error_d1_roll);
% update integrator
differentiator_roll = (2*tau-Ts)/(2*tau+Ts)*differentiator_roll...
    + 2/(2*tau+Ts)*(error - error_d1_roll);
% update differentiator

% the loop
% implement PID control

u = sat(...
    kp * error +...
    ki * integrator_roll +...
    kd * differentiator_roll,...
    limit...
    );

% implement integrator anti-windup
if abs(error_d1_roll-error)<TOL_Error
    u_unsat = kp*error + ki*integrator_roll + kd*differentiator_roll;
    integrator_roll = integrator_roll + Ts/ki * (u - u_unsat);
end

error_d1_roll = error; % update the error for next time through

    function out = sat(in, limit)
        if in > limit
            out = limit;
        elseif in < -limit
            out = -limit;
        else
            out = in;
        end
    end
end

function u = course_hold(y_c, y, ~, P,myflag,TOL_Error)

kp=P.kp_chi;
ki=P.ki_chi;
kd=P.kd_chi;
limit=P.delta_a_max;
Ts=P.Ts;
tau=P.tau;

persistent integrator_course;
persistent differentiator_course;
persistent error_d1_course;

if myflag==1 % reset (initialize) persistent variables
    % when myflag==1
    integrator_course = 0;
    differentiator_course = 0;
    error_d1_course = 0; % _d1 means delayed by one time step
end
error = y_c - y; % compute the current error
integrator_course = integrator_course + (Ts/2)*(error + error_d1_course);
% update integrator
differentiator_course = (2*tau-Ts)/(2*tau+Ts)*differentiator_course...
    + 2/(2*tau+Ts)*(error - error_d1_course);
% update differentiator

% the loop
% implement PID control

u = sat(...
    kp * error +...
    ki * integrator_course +...
    kd * differentiator_course,...
    limit...
    );

% implement integrator anti-windup
if abs(error_d1_course-error)<TOL_Error
    u_unsat = kp*error + ki*integrator_course + kd*differentiator_course;
    integrator_course = integrator_course + Ts/ki * (u - u_unsat);
end

error_d1_course = error; % update the error for next time through

    function out = sat(in, limit)
        if in > limit
            out = limit;
        elseif in < -limit
            out = -limit;
        else
            out = in;
        end
    end
end


function u = pitch_hold(y_c, y, ~, P,myflag,TOL_Error)

kp=P.kp_theta;
ki=P.ki_theta;
kd=P.kd_theta;
limit=P.delta_e_max;
Ts=P.Ts;
tau=P.tau;

persistent integrator_pitch;
persistent differentiator_pitch;
persistent error_d1_pitch;

if myflag==1 % reset (initialize) persistent variables
    % when myflag==1
    integrator_pitch = 0;
    differentiator_pitch = 0;
    error_d1_pitch = 0; % _d1 means delayed by one time step
end
error = y_c - y; % compute the current error
integrator_pitch = integrator_pitch + (Ts/2)*(error + error_d1_pitch);
% update integrator
differentiator_pitch = (2*tau-Ts)/(2*tau+Ts)*differentiator_pitch...
    + 2/(2*tau+Ts)*(error - error_d1_pitch);
% update differentiator

% the loop
% implement PID control

u = sat(...
    kp * error +...
    ki * integrator_pitch +...
    kd * differentiator_pitch,...
    limit...
    );

% implement integrator anti-windup
if abs(error_d1_pitch-error)<TOL_Error
    u_unsat = kp*error + ki*integrator_pitch + kd*differentiator_pitch;
    integrator_pitch = integrator_pitch + Ts/ki * (u - u_unsat);
end

error_d1_pitch = error; % update the error for next time through

    function out = sat(in, limit)
        if in > limit
            out = limit;
        elseif in < -limit
            out = -limit;
        else
            out = in;
        end
    end
end


function u = airspeed_with_throttle_hold(y_c, y, ~, P,myflag,TOL_Error)

kp=P.kp_V;
ki=P.ki_V;
kd=P.kd_V;
limit=P.delta_t_max;
Ts=P.Ts;
tau=P.tau;

persistent integrator_air_throttle;
persistent differentiator_air_throttle;
persistent error_d1_air_throttle;

if myflag==1 % reset (initialize) persistent variables
    % when myflag==1
    integrator_air_throttle = 0;
    differentiator_air_throttle = 0;
    error_d1_air_throttle = 0; % _d1 means delayed by one time step
end
error = y_c - y; % compute the current error
integrator_air_throttle = integrator_air_throttle + (Ts/2)*(error + error_d1_air_throttle);
% update integrator
differentiator_air_throttle = (2*tau-Ts)/(2*tau+Ts)*differentiator_air_throttle...
    + 2/(2*tau+Ts)*(error - error_d1_air_throttle);
% update differentiator

% the loop
% implement PID control

u = sat(...
    kp * error +...
    ki * integrator_air_throttle +...
    kd * differentiator_air_throttle,...
    limit...
    );

% implement integrator anti-windup
if abs(error_d1_air_throttle-error)<TOL_Error
    u_unsat = kp*error + ki*integrator_air_throttle + kd*differentiator_air_throttle;
    integrator_air_throttle = integrator_air_throttle + Ts/ki * (u - u_unsat);
end

error_d1_air_throttle = error; % update the error for next time through

    function out = sat(in, limit)
        if in > limit
            out = limit;
        elseif in < 0
            out = 0;
        else
            out = in;
        end
    end
end

function u = airspeed_with_pitch_hold(y_c, y, ~, P,myflag,TOL_Error)

kp=P.kp_V2;
ki=P.ki_V2;
kd=P.kd_V2;
limit=P.theta_max;
Ts=P.Ts;
tau=P.tau;

persistent integrator_air_pitch;
persistent differentiator_air_pitch;
persistent error_d1_air_pitch;

if myflag==1 % reset (initialize) persistent variables
    % when myflag==1
    integrator_air_pitch = 0;
    differentiator_air_pitch = 0;
    error_d1_air_pitch = 0; % _d1 means delayed by one time step
end
error = y_c - y; % compute the current error
integrator_air_pitch = integrator_air_pitch + (Ts/2)*(error + error_d1_air_pitch);
% update integrator
differentiator_air_pitch = (2*tau-Ts)/(2*tau+Ts)*differentiator_air_pitch...
    + 2/(2*tau+Ts)*(error - error_d1_air_pitch);
% update differentiator

% the loop
% implement PID control

u = sat(...
    kp * error +...
    ki * integrator_air_pitch +...
    kd * differentiator_air_pitch,...
    limit...
    );


% implement integrator anti-windup
if abs(error_d1_air_pitch-error)<TOL_Error
    u_unsat = kp*error + ki*integrator_air_pitch + kd*differentiator_air_pitch;
    integrator_air_pitch = integrator_air_pitch + Ts/ki * (u - u_unsat);
end

error_d1_air_pitch = error; % update the error for next time through

    function out = sat(in, limit)
        if in > limit
            out = limit;
        elseif in < -limit
            out = -limit;
        else
            out = in;
        end
    end
end

function u = altitude_hold(y_c, y, ~, P,myflag,TOL_Error)

kp=P.kp_h;
ki=P.ki_h;
kd=P.kd_h;
limit=P.theta_max;
Ts=P.Ts;
tau=P.tau;

persistent integrator_altitude_pitch;
persistent differentiator_air_pitch;
persistent error_d1_air_pitch;

if myflag==1 % reset (initialize) persistent variables
    % when myflag==1
    integrator_altitude_pitch = 0;
    differentiator_air_pitch = 0;
    error_d1_air_pitch = 0; % _d1 means delayed by one time step
   
end
error = y_c - y; % compute the current error

differentiator_air_pitch = (2*tau-Ts)/(2*tau+Ts)*differentiator_air_pitch...
    + 2/(2*tau+Ts)*(error - error_d1_air_pitch);

integrator_altitude_pitch = integrator_altitude_pitch + (Ts/2)*(error + error_d1_air_pitch);
%     % update integrator


% update differentiator

% the loop
% implement PID control


u = sat(...
    kp * error +...
    ki * integrator_altitude_pitch +...
    kd * differentiator_air_pitch,...
    limit...
    );


     % implement integrator anti-windup
     if ki~=0.0
     u_unsat = kp*error + ki*integrator_altitude_pitch + kd*differentiator_air_pitch;
     integrator_altitude_pitch = integrator_altitude_pitch + Ts/ki * (u - u_unsat);
     end

error_d1_air_pitch = error; % update the error for next time through

    function out = sat(in, limit)
        if in > limit
            out = limit;
        elseif in < -limit
            out = -limit;
        else
            out = in;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopilot_tuning
%   - used to tune each loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta, x_command] = autopilot_tuning(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error)

myflag = 0;
if t==0
    myflag = 1;
end

mode = 5;
switch mode
    case 1 % tune the roll loop
        phi_c = chi_c; % interpret chi_c to autopilot as course command
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        % use trim values for elevator and throttle while tuning the lateral autopilot
        delta_e = P.u_trim(1);
        delta_t = P.u_trim(4);
        theta_c = 0;
    case 2 % tune the course loop
        
        phi_c   = course_hold(chi_c, chi, r, P,myflag,TOL_Error);
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        % use trim values for elevator and throttle while tuning the lateral autopilot
        delta_e = P.u_trim(1);
        delta_t = P.u_trim(4);
        theta_c = 0;
    case 3 % tune the throttle to airspeed loop and pitch loop simultaneously
        theta_c = 20*pi/180 + h_c;
        chi_c = 0;
        phi_c   = course_hold(chi_c, chi, r, P,myflag,TOL_Error);
        delta_t = airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        % use trim values for elevator and throttle while tuning the lateral autopilot
    case 4 % tune the pitch to airspeed loop
        chi_c = 0;
        delta_t = P.u_trim(4);
        
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error);
        theta_c = airspeed_with_pitch_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        % use trim values for elevator and throttle while tuning the lateral autopilot
    case 5 % tune the pitch to altitude loop
        chi_c = 0;
        
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error);
        theta_c = altitude_hold(h_c, h,0, P,myflag,TOL_Error);
        delta_t = airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        % use trim values for elevator and throttle while tuning the lateral autopilot
end
%----------------------------------------------------------
% create outputs

% control outputs
delta = [delta_e; delta_a; delta_r; delta_t];
% commanded (desired) states
x_command = [...
    0;...                    % pn
    0;...                    % pe
    h_c;...                  % h
    Va_c;...                 % Va
    0;...                    % alpha
    0;...                    % beta
    phi_c;...                % phi
    %theta_c*P.K_theta_DC;... % theta
    theta_c;
    chi_c;...                % chi
    0;...                    % p
    0;...                    % q
    0;...                    % r
    ];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopilot_uavbook
%   - autopilot defined in the uavbook
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta, x_command] = autopilot_uavbook(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error)

myflag = 0;
if t==0
    myflag = 1;
end

% define persistent variable for state of altitude state machine
persistent altitude_state;
% initialize persistent variable

if h<=P.altitude_take_off_zone
    altitude_state = 1;
elseif h<=h_c-P.altitude_hold_zone
    altitude_state = 2;
elseif h>=h_c+P.altitude_hold_zone
    altitude_state = 3;
else
    altitude_state = 4;
end


% implement state machine
switch altitude_state
    case 1  % in take-off zone
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error);
        theta_c = P.theta_max;%altitude_hold(h_c, h,0, P,myflag,TOL_Error);
        Va_c = P.V_max;
        delta_t = airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        
    case 2  % climb zone
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error);
        theta_c = altitude_hold(h_c, h,0, P,myflag,TOL_Error);
        Va_c = P.V_max;
        delta_t = airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        
    case 3 % descend zone
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error);
        theta_c = altitude_hold(h_c, h,0, P,myflag,TOL_Error);
        delta_t = 0.0;%airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
        
    case 4 % altitude hold zone
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error);
        theta_c = altitude_hold(h_c, h,0, P,myflag,TOL_Error);
        delta_t = airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder
end



%----------------------------------------------------------
% create outputs

% control outputs
delta = [delta_e; delta_a; delta_r; delta_t];
% commanded (desired) states
x_command = [...
    0;...                    % pn
    0;...                    % pe
    h_c;...                  % h
    Va_c;...                 % Va
    0;...                    % alpha
    0;...                    % beta
    phi_c;...                % phi
    %theta_c*P.K_theta_DC;... % theta
    theta_c;
    chi_c;...                % chi
    0;...                    % p
    0;...                    % q
    0;...                    % r
    ];

y = [delta; x_command];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopilot_TECS
%   - longitudinal autopilot based on total energy control systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta, x_command] = autopilot_TECS(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P)

%----------------------------------------------------------
% lateral autopilot
if t==0,
    % assume no rudder, therefore set delta_r=0
    delta_r = 0;%coordinated_turn_hold(beta, 1, P);
    phi_c   = course_hold(chi_c, chi, r, 1, P);
    
else
    phi_c   = course_hold(chi_c, chi, r, 0, P);
    delta_r = 0;%coordinated_turn_hold(beta, 0, P);
end
delta_a = roll_hold(phi_c, phi, p, P);


%----------------------------------------------------------
% longitudinal autopilot based on total energy control


delta_e = 0;
delta_t = 0;


%----------------------------------------------------------
% create outputs

% control outputs
delta = [delta_e; delta_a; delta_r; delta_t];
% commanded (desired) states
x_command = [...
    0;...                    % pn
    0;...                    % pe
    h_c;...                  % h
    Va_c;...                 % Va
    0;...                    % alpha
    0;...                    % beta
    phi_c;...                % phi
    %theta_c*P.K_theta_DC;... % theta
    theta_c;
    chi_c;...                % chi
    0;...                    % p
    0;...                    % q
    0;...                    % r
    ];

y = [delta; x_command];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autopilot functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [delta, x_command] = autopilot_kevin(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P,TOL_Error)

myflag = 0;
if t<1E-4
    myflag = 1;
end

%         chi_c = 0.0
        phi_c   = course_hold(chi_c, chi, r,P,myflag,TOL_Error*10);
        theta_c = altitude_hold(h_c, h,0, P,myflag,TOL_Error);
        delta_t = airspeed_with_throttle_hold(Va_c, Va,0, P,myflag,TOL_Error);
        
        
        delta_a = roll_hold(phi_c, phi, p, P,myflag,TOL_Error);
        delta_e = pitch_hold(theta_c, theta, q, P,myflag,TOL_Error);
        delta_r = 0; % no rudder




%----------------------------------------------------------
% create outputs

% control outputs
delta = [delta_e; delta_a; delta_r; delta_t];
% commanded (desired) states
x_command = [...
    0;...                    % pn
    0;...                    % pe
    h_c;...                  % h
    Va_c;...                 % Va
    0;...                    % alpha
    0;...                    % beta
    phi_c;...                % phi
    %theta_c*P.K_theta_DC;... % theta
    theta_c;
    chi_c;...                % chi
    0;...                    % p
    0;...                    % q
    0;...                    % r
    ];

y = [delta; x_command];

end
