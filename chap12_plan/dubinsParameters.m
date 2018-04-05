% dubinsParameters
%   - Find Dubin's parameters between two configurations'
%
% input is:
%   start_node  - [wn_s, we_s, wd_s, chi_s, 0, 0]
%   end_node    - [wn_e, wn_e, wd_e, chi_e, 0, 0]
%   R           - minimum turn radius
%
% output is:
%   dubinspath  - a matlab structure with the following fields
%       dubinspath.ps   - the start position in re^3
%       dubinspath.chis - the start course angle
%       dubinspath.pe   - the end position in re^3
%       dubinspath.chie - the end course angle
%       dubinspath.R    - turn radius
%       dubinspath.L    - length of the Dubins path
%       dubinspath.cs   - center of the start circle
%       dubinspath.lams - direction of the start circle
%       dubinspath.ce   - center of the end circle
%       dubinspath.lame - direction of the end circle
%       dubinspath.w1   - vector in re^3 defining half plane H1
%       dubinspath.q1   - unit vector in re^3 along straight line path
%       dubinspath.w2   - vector in re^3 defining position of half plane H2
%       dubinspath.w3   - vector in re^3 defining position of half plane H3
%       dubinspath.q3   - unit vector defining direction of half plane H3
%

function dubinspath = dubinsParameters(start_node, end_node, Rmin)

ell = norm(start_node(1:2)-end_node(1:2));
if 0%ell<2*Rmin
    disp('The distance between nodes must be larger than 2R.');
    dubinspath = [];
else

    ps   = start_node(1:3)';
    chis = start_node(4);
    pe   = end_node(1:3)';
    chie = end_node(4);


    crs = ps+Rmin*[cos(chis+pi/2);sin(chis+pi/2);0];
    cls = ps+Rmin*[cos(chis-pi/2);sin(chis-pi/2);0];
    cre = pe+Rmin*[cos(chie+pi/2);sin(chie+pi/2);0];
    cle = pe+Rmin*[cos(chie-pi/2);sin(chie-pi/2);0];


    % compute L1
    %     u = [1 2 0];
    %     v = [1 0 0];
    %     ThetaInDegrees = atan2(norm(cross(u,v)),dot(u,v));
    diffpspe = crs-cre;
    north = diffpspe(1);
    east = diffpspe(2);
    theta = atan2(east,north);
    L1 = norm(crs-cre)+Rmin*mod(2*pi,+mod(theta,-pi/2)-mod(chis,-pi/2))+Rmin*mod(2*pi,mod(chie,-pi/2)-mod(theta,-pi/2));
    % compute L2
    diffpspe = cle-cls;
    north = diffpspe(1);
    east = diffpspe(2);
    theta = atan2(east,north);
    ell = norm(cle-crs);
    theta2 = theta-pi/2+asin(2*Rmin/ell);
    if isreal(theta2)==0
        L2 = 9999;
    else
        L2 = sqrt((ell)^2-4*Rmin^2)+Rmin*mod(2*pi,theta2-mod(chis,-pi/2))+Rmin*mod(2*pi,mod(theta2,pi)-mod(chie,pi/2));
    end
    % compute L3
    diffpspe = cls-cre;
    north = diffpspe(1);
    east = diffpspe(2);
    theta = atan2(east,north);
    ell = norm(cre-cls);
    theta2 = acos(2*Rmin/ell);
    if isreal(theta2)==0
        L3 = 9999;
    else
        L3 = sqrt((ell)^2-4*Rmin^2)+Rmin*mod(2*pi,mod(chis,pi/2)-mod(theta,+theta2))+Rmin*mod(2*pi,mod(chie,-pi/2)-mod(theta,theta2-pi));
    end
    % compute L4
    diffpspe = cls-cle;
    north = diffpspe(1);
    east = diffpspe(2);
    theta = atan2(east,north);
    L4 = norm(cls-cle)+Rmin*mod(2*pi,mod(chis,pi/2)-mod(theta,pi/2)) + Rmin*mod(2*pi,mod(theta,pi/2)-mod(chie,pi/2));
    % L is the minimum distance
    [L,idx] = min([L1,L2,L3,L4]);
    e1 = [1; 0; 0];
    switch(idx)
        case 1
            cs = crs;
            lams = 1;
            ce = cre;
            lame = 1;
            q1 = (ce-cs)/norm(ce-cs);
            w1 = cs+Rmin*rotz(-pi/2)*q1;
            w2 = ce+Rmin*rotz(-pi/2)*q1;
        case 2
            cs = crs;
            lams = 1;
            ce = cle;
            lame = -1;
            diffpspe = ce-cs;
            north = diffpspe(1);
            east = diffpspe(2);
            theta = atan2(east,north);
            ell = norm(ce-cs);
            theta2 = theta-pi/2+asin(2*Rmin/ell);
            %             ell = norm(ce-cs);
            %             theta = angle(ce-cs);
            %             theta2 = theta-pi/2+asin(2*Rmin/ell);
            q1 = rotz(theta2+pi/2)*e1;
            w1 = cs+Rmin*rotz(theta2)*e1;
            w2 = ce+Rmin*rotz(theta2+pi)*e1;
        case 3
            cs = cls;
            lams = -1;
            ce = cre;
            lame = 1;

            diffpspe = ce-cs;
            north = diffpspe(1);
            east = diffpspe(2);
            theta = atan2(east,north);
            ell = norm(ce-cs);
            theta2 = acos(2*Rmin/ell);
            %             ell = norm(ce-cs);
            %             theta = angle(ce-cs);
            %             theta2 = acos(2*Rmin/ell);
            q1 = rotz(theta+theta2-pi/2)*e1;
            w1 = cs+Rmin*rotz(theta+theta2)*e1;
            w2 = ce+Rmin*rotz(theta+theta2-pi)*e1;
        case 4
            cs = cls;
            lams = -1;
            ce = cle;
            lame = -1;
            q1 = (ce-cs)/norm(ce-cs);
            w1 = cs+Rmin*rotz(pi/2)*q1;
            q2 = q1;
            w2 = ce+Rmin*rotz(pi/2)*q2;
    end
    w3 = pe;
    q3 = rotz(chie)*e1;


    % assign path variables
    dubinspath.ps   = ps;
    dubinspath.chis = chis;
    dubinspath.pe   = pe;
    dubinspath.chie = chie;
    dubinspath.R    = Rmin;
    dubinspath.L    = L;
    dubinspath.cs   = cs;
    dubinspath.lams = lams;
    dubinspath.ce   = ce;
    dubinspath.lame = lame;
    dubinspath.w1   = w1;
    dubinspath.q1   = q1;
    dubinspath.w2   = w2;
    dubinspath.w3   = w3;
    dubinspath.q3   = q3;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotz(theta)
%   rotation matrix about the z axis.
function R = rotz(theta)
R = [...
    cos(theta), -sin(theta), 0;...
    sin(theta), cos(theta), 0;...
    0, 0, 1;...
    ];
end