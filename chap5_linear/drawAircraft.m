
function drawAircraft(uu,V,F,patchcolors)

    % process inputs to function
    pn       = uu(1);       % inertial North position
    pe       = uu(2);       % inertial East position
    pd       = uu(3);
    u        = uu(4);
    v        = uu(5);
    w        = uu(6);
    phi      = uu(7);       % roll angle
    theta    = uu(8);       % pitch angle
    psi      = uu(9);       % yaw angle
    p        = uu(10);       % roll rate
    q        = uu(11);       % pitch rate
    r        = uu(12);       % yaw rate
    t        = uu(13);       % time

    % define persistent variables
    persistent vehicle_handle;
    persistent Vertices
    persistent Faces
    persistent facecolors

    % first time function is called, initialize plot and persistent vars
    if t==0,
        figure(1), clf
        [Vertices,Faces,facecolors] = defineVehicleBody;
        vehicle_handle = drawVehicleBody(Vertices,Faces,facecolors,...
                                               pn,pe,pd,phi,theta,psi,...
                                               [],'normal');
        title('Vehicle')
        xlabel('East')
        ylabel('North')
        zlabel('-Down')
        view(90,0)  % set the vieew angle for figure
        axis([-500,500,-200,1000,-500,500]);
%         axis equal
        hold on

    % at every other time step, redraw base and rod
    else
        drawVehicleBody(Vertices,Faces,facecolors,...
                           pn,pe,pd,phi,theta,psi,...
                           vehicle_handle);
    end
end


%=======================================================================
% drawVehicle
% return handle if 3rd argument is empty, otherwise use 3rd arg as handle
%=======================================================================
%
function handle = drawVehicleBody(V,F,patchcolors,...
                                     pn,pe,pd,phi,theta,psi,...
                                     handle,mode)
  V = rotate(V, 0.0, 0.0, pi);  % rotate vehicle to face the correct direction as drawn
  V = translate(V,1.5, 0.0, 0.0);  % translate vehicle so
  V = rotate(V, phi, theta, psi);  % rotate vehicle
  V = translate(V, pn, pe, pd);  % translate vehicle
  % transform vertices from NED to XYZ (for matlab rendering)
  R = [...
      0, 1, 0;...
      1, 0, 0;...
      0, 0, -1;...
      ];
  V = R*V;

  if isempty(handle),
  handle = patch('Vertices', V', 'Faces', F,...
                 'FaceVertexCData',patchcolors,...
                 'FaceColor','flat',...
                 'EraseMode', mode);
  else
    set(handle,'Vertices',V','Faces',F);
    drawnow
  end
end

%%%%%%%%%%%%%%%%%%%%%%%
function pts=rotate(pts,phi,theta,psi)

  % define rotation matrix (right handed)
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
  R = R';

  % rotate vertices
  pts = R*pts;

end
% end rotateVert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate vertices by pn, pe, pd
function pts = translate(pts,pn,pe,pd)

  pts = pts + repmat([pn;pe;pd],1,size(pts,2));

end

% end translate


%=======================================================================
% defineVehicleBody
%=======================================================================
function [V_plane, F_plane, colors_plane] = defineVehicleBody




    % planebird parameters
    scale = 50.0;
    nose_w = 1.0*scale;
    nose_h = 1.0*scale;
    nose_l = 1.0*scale;
    nose_scale = [1.0,3.0,3.0];
    tail_scale = [1.0,4.0,4.0];
    fuse_l = 5.0*scale;
    wing_Xmount = 1.0*scale;
    wing_Zmount = 0.0*scale;
    wing_chord = 1.0*scale;
    wing_span = 7.0*scale;
    tail_Xmount = 4.0*scale;
    tail_Zmount = 0.0*scale;
    tail_tipchord = 1.0*scale;
    tail_span = 2.0*scale;
    rudder_mount = 4.0*scale;
    rudder_chord = 1.0*scale;
    rudder_h = 1.0*scale;
    w = 0.25*scale;   % width of the center pod
    l = .4*scale;   % length of connecting rods
    lw = 0.05*scale; % width of rod
    r = 1.0*scale;   % radius of the rotor
    N = 10;     % number of points defining rotor

	p1 = [0,-nose_w/2,-nose_h/2]./nose_scale;
    p2 = [0,nose_w/2,-nose_h/2]./nose_scale;
    p3 = [0,nose_w/2,nose_h/2]./nose_scale;
    p4 = [0,-nose_w/2,nose_h/2]./nose_scale;
    p5 = [nose_l,-nose_w/2,-nose_h/2];
    p6 = [nose_l,nose_w/2,-nose_h/2];
    p7 = [nose_l,nose_w/2,nose_h/2];
    p8 = [nose_l,-nose_w/2,nose_h/2];
    p9 = [wing_Xmount,-wing_span/2,wing_Zmount];
    p10 = [wing_Xmount+wing_chord,-wing_span/2,wing_Zmount];
    p11 = [wing_Xmount+wing_chord,wing_span/2,wing_Zmount];
    p12 = [wing_Xmount,wing_span/2,wing_Zmount];
    p13 = [tail_Xmount,-tail_span/2,tail_Zmount];
    p14 = [tail_Xmount+tail_tipchord,-tail_span/2,tail_Zmount];
    p15 = [tail_Xmount+tail_tipchord,tail_span/2,tail_Zmount];
    p16 = [tail_Xmount,tail_span/2,tail_Zmount];
    p17 = [rudder_mount,0,tail_Zmount];
    p18 = [rudder_mount,0,-rudder_h];
    p19 = [rudder_mount+rudder_chord,0,-rudder_h];
    p20 = [rudder_mount+rudder_chord,0,tail_Zmount];
	p21 = [fuse_l,-nose_w/2,-nose_h/2]./tail_scale;
    p22 = [fuse_l,nose_w/2,-nose_h/2]./tail_scale;
    p23 = [fuse_l,nose_w/2,nose_h/2]./tail_scale;
    p24 = [fuse_l,-nose_w/2,nose_h/2]./tail_scale;

    % define colors for faces
    red    = [1, 0, 0];
    green  = [0, 1, 0];
    blue   = [0, 0, 1];
    rodcolor = [.1, 0, .5];
    yellow = [1,1,0];


% define vertices, faces, and colors for planebird body

  %--------- vertices and faces for center pod ------------
  % vertices of the center pod
  V_center = [...
    p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19;p20;p21;p22;p23;p24
    ];
  % define faces of center pod
  F_center = [...
    1,2,3,4;...
    1,5,6,2;...
    2,6,7,3;...
    3,7,8,4;...
    1,5,8,4;...
    5,21,22,6;...
    6,22,23,7;...
    7,23,24,8;...
    8,24,21,5;...
    9,10,11,12;...
    13,14,15,16;...
    17,18,19,20;...
        ];

  %--------- vertices and faces for connecting rods ------------
  % vertices for right rod
  V_rod_front = [...
    w/2, lw/2, 0;...
    w/2, -lw/2, 0;...
    w/2, 0, lw/2;...
    w/2, 0, -lw/2;...
    l+w/2, lw/2, 0;...
    l+w/2, -lw/2, 0;...
    l+w/2, 0, lw/2;...
    l+w/2, 0, -lw/2;...
    ];
  V_rod_right = V_rod_front*[0,1,0;1,0,0;0,0,1];
  F_rod_right = 8 + [...
        1, 2, 6, 5;... % x-y face
        3, 4, 8, 7;... % x-z face
    ];


  %--------- vertices and faces for rotors ------------
  V_rotor = [];
  for i=1:N
    V_rotor = [V_rotor; r*cos(2*pi*i/N), r*sin(2*pi*i/N),0];
  end
  for i=1:N
    V_rotor = [V_rotor; r/10*cos(2*pi*i/N), r/10*sin(2*pi*i/N),0];
  end
  F_rotor = [];
  for i=1:N-1
    F_rotor = [F_rotor; i, i+1, N+i+1, N+i];
  end
  F_rotor = [F_rotor; N, 1, N+1, 2*N];

  V_rotor1 = V_rotor + repmat([0, 0, 0],2*N,1);
  V_rotor = [V_rotor1(:,3),V_rotor1(:,2),V_rotor1(:,1)];
  F_rotor = 24 + F_rotor;


  % collect all of the vertices for the quadrotor into one matrix
  V_plane = [...
    V_center; V_rotor...
    ];
  % collect all of the faces for the quadrotor into one matrix
  F_plane = [...
    F_center; F_rotor...
    ];


  colors_plane = [...
    blue;blue;blue;blue;blue;blue;blue;blue;blue;blue;blue;blue;
%     green;... % rod left
%     green;... %
    ];
  for i=1:N
    colors_plane = [colors_plane; green];  % rotor
  end

  V_plane= V_plane';
%   F_plane = F_plane';
%   colors_plane = colors_plane';

end
