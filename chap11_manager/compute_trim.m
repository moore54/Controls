function [x_trim,u_trim] = compute_trim(filename, Va, gamma, R)
% Va is the desired airspeed (m/s)
% gamma is the desired flight path angle (radians)
% R is the desired radius (m) - use (+) for right handed orbit, 
%                                   (-) for left handed orbit


if (R~=0)
    psidot = Va/R;
else
    psidot = 0;
end

% add stuff here
X0 = [0; 0; 0; Va; 0; 0; 0; gamma; 0; 0; 0; 0] ;
IX0 = [] ;
U0 = [0; 0; 0; 1] ;
IU0 = [] ;
Y0 = [Va; gamma; 0] ;
IY0 = [1,3];

DX0 = [0; 0;-Va*sin(gamma); 0; 0; 0;0;0; psidot;0;0;0]; 
IDX = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12];

% compute trim conditions
[x_trim,u_trim,y_trim,dx_trim] = trim(filename,X0,U0,Y0,IX0,IU0,IY0,DX0,IDX);

% check to make sure that the linearization worked (should be small)
norm(dx_trim(3:end)-X0(3:end))

