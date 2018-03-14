function [A_lon,B_lon,A_lat,B_lat] = compute_ss_model(filename,x_trim,u_trim)
% x_trim is the trimmed state,
% u_trim is the trimmed input
  
% add stuff here  
[A,B,C,D]=linmod(filename,x_trim,u_trim);

E1 = [0,0,0,0,1,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,1,0,0;...
    0,0,0,0,0,0,0,0,0,0,0,1;...
    0,0,0,0,0,0,1,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,1,0,0,0;...
    ];

E2 = [0,1,0,0;...
    0,0,1,0];

E1_lon = [0,0,0,1,0,0,0,0,0,0,0,0;...
    0,0,0,0,0,1,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0,0,0,1,0;...
    0,0,0,0,0,0,0,1,0,0,0,0;...
    0,0,-1,0,1,0,0,0,0,0,0,0;...
    ];

E2_lon = [1,0,0,0;...
    0,0,0,1];


A_lat = E1 * A * E1';
B_lat = E1 * B * E2';
A_lon = E1_lon * A * E1_lon';
B_lon = E1_lon * B * E2_lon';

% eig_lat = eig(A_lat);
% eig_lon = eig(A_lon);