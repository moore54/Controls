% gps.m
%   Compute the output of gps sensor
%
%  Revised:
%   3/5/2010 - RB 
%   5/14/2010 - RB

function y = gps(uu, P)

    % relabel the inputs
    Va      = uu(1);
%    alpha   = uu(2);
%    beta    = uu(3);
    wn      = uu(4);
    we      = uu(5);
%    wd      = uu(6);
    pn      = uu(7);
    pe      = uu(8);
    pd      = uu(9);
%    u       = uu(10);
%    v       = uu(11);
%    w       = uu(12);
%    phi     = uu(13);
%    theta   = uu(14);
    psi     = uu(15);
%    p       = uu(16);
%    q       = uu(17);
%    r       = uu(18);
    t       = uu(19);
    
    persistent vn;
    persistent ve;
    persistent vd;
    

    if t==0
        vn = 0.0;
        ve = 0.0;
        vd = 0.0;
        

    end

    vn = exp(-P.k_GPS*P.Ts)*vn;
    ve = exp(-P.k_GPS*P.Ts)*ve;
    vd = exp(-P.k_GPS*P.Ts)*vd;
    
    % construct North, East, and altitude GPS measurements
    y_gps_n = pn+vn+P.B_GPS_n+randn*P.sigma_GPS_n;
    y_gps_e = pe+ve+P.B_GPS_e+randn*P.sigma_GPS_e; 
    y_gps_h = pd+vd+P.B_GPS_h+randn*P.sigma_GPS_h; 
    
    % construct groundspeed and course measurements

    y_gps_Vg     = sqrt((Va*cos(psi)+wn)^2+(Va*sin(psi)+we)^2)+randn*P.sigma_GPS_Vg;
    y_gps_course = atan2(Va*sin(psi)+we,Va*cos(psi)+wn)+randn*P.sigma_GPS_course;

    % construct total output
    y = [...
        y_gps_n;...
        y_gps_e;...
        y_gps_h;...
        y_gps_Vg;...
        y_gps_course;...
        ];
    
end



