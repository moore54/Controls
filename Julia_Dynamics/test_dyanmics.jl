using DynamicalSystems, StaticArrays, PyPlot
PyPlot.close("all")
# Lorenz system
# Equations of motion:
@inline @inbounds function loop(u, p, t)
    σ = p[1]; ρ = p[2]; β = p[3]
    du1 = σ*(u[2]-u[1])
    du2 = u[1]*(ρ-u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    return SVector{3}(du1, du2, du3)
end
# Jacobian:
@inline @inbounds function loop_jac(u, p, t)
    σ, ρ, β = p
    J = @SMatrix [-σ  σ  0;
    ρ - u[3]  (-1)  (-u[1]);
    u[2]   u[1]  -β]
    return J
end

@inline @inbounds function forces_moments()
    # relabel the inputs
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
    w_ns    = wind(1); # steady wind - North
    w_es    = wind(2); # steady wind - East
    w_ds    = wind(3); # steady wind - Down
    u_wg    = wind(4); # gust along body x-axis
    v_wg    = wind(5); # gust along body y-axis
    w_wg    = wind(6); # gust along body z-axis

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
      # note that R above either leaves the vector alone or rotates
      # a vector in a left handed rotation.  We want to rotate all
      # points in a right handed rotation, so we must transpose
    R = R';#'

    # compute wind data in NED $TODO: include gusts
    W_gustBODY = [u_wg;v_wg;w_wg];
#     W_gustNED = W_gustBODY\R
    w_n = w_ns;#+W_gustNED(1);
    w_e = w_es;#+W_gustNED(2);
    w_d = w_ds;#+W_gustNED(3);

    W_BODY = R*[w_n;w_e;w_d]+[W_gustBODY(1);W_gustBODY(2);W_gustBODY(3)];

    # compute air data in body frame
    ur = u-W_BODY(1);
    vr = v-W_BODY(2);
    wr = w-W_BODY(3);

    # compute total wind in ground frame
    W_vec = R'*W_BODY;#'
    w_n = W_vec(1);
    w_e = W_vec(2);
    w_d = W_vec(3);

    Va = sqrt(ur^2+vr^2+wr^2);
    alpha = atan2(wr,ur);#-P[:alpha0];
    beta = asin(vr/Va);

    if (Va == 0)
        Va = P[:Va0];
    end
    if ~isfinite(alpha)
        alpha = 0;
    end
    if ~isfinite(beta)
        beta = 0;
    end


   C_D_alpha = P[:C_D_p] + ((P[:C_L_0] + P[:C_L_alpha]*alpha)^2)/(pi*P[:e]*P[:AR]);
   sigma_alpha = (1 + exp(-P[:M]*(alpha-P[:alpha0])) + exp(P[:M]*(alpha+P[:alpha0])))/((1+exp(-P[:M]*(alpha-P[:alpha0])))*(1+exp(P[:M]*(alpha+P[:alpha0]))));
   C_L_alpha = (1-sigma_alpha)*(P[:C_L_0] + P[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));

   C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
   C_X_q_alpha = -P[:C_D_q]*cos(alpha) + P[:C_L_q]*sin(alpha);
   C_X_delta_e_alpha = -P[:C_D_delta_e]*cos(alpha) + P[:C_L_delta_e]*sin(alpha);
   C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);
   C_Z_q_alpha = -P[:C_D_q]*sin(alpha) - P[:C_L_q]*cos(alpha);
   C_Z_delta_e_alpha = -P[:C_D_delta_e]*sin(alpha) - P[:C_L_delta_e]*cos(alpha);

    # compute external forces and torques on aircraft
    Force(1) =  -P[:mass]*P[:gravity]*sin(theta)+0.5*P[:rho]*Va^2*P[:S_wing]*(C_X_alpha+C_X_q_alpha*P[:c]/(2*Va)*q+C_X_delta_e_alpha*delta_e)...
        +0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*delta_t)^2-Va^2);
    Force(2) =  P[:mass]*P[:gravity]*cos(theta)*sin(phi)+0.5*P[:rho]*Va^2*P[:S_wing]*(P[:C_Y_0]+P[:C_Y_beta]*beta+P[:C_Y_p]*P[:b]/(2*Va)*p+...
        P[:C_Y_r]*P[:b]/(2*Va)*r+P[:C_Y_delta_a]*delta_a+P[:C_Y_delta_r]*delta_r);
    Force(3) =  P[:mass]*P[:gravity]*cos(theta)*cos(phi)+0.5*P[:rho]*Va^2*P[:S_wing]*(C_Z_alpha+C_Z_q_alpha*P[:c]/(2*Va)*q+C_Z_delta_e_alpha*delta_e);

    Torque(1) = 0.5*P[:rho]*Va^2*P[:S_wing]*P[:b]*(P[:C_ell_0] + P[:C_ell_beta]*beta + (P[:C_ell_p]*P[:b]*p/(2*Va)) + (P[:C_ell_r]*P[:b]*r/(2*Va)) + P[:C_ell_delta_a]*delta_a + P[:C_ell_delta_r]*delta_r)...
        - P[:k_T_P]*(P[:k_Omega]*delta_t)^2;
    Torque(2) = 0.5*P[:rho]*Va^2*P[:S_wing]*P[:c]*(P[:C_m_0] + P[:C_m_alpha]*alpha + (P[:C_m_q]*P[:c]*q/(2*Va)) + P[:C_m_delta_e]*delta_e);
    Torque(3) = 0.5*P[:rho]*Va^2*P[:S_wing]*P[:b]*(P[:C_n_0] + P[:C_n_beta]*beta + (P[:C_n_p]*P[:b]*p/(2*Va)) + (P[:C_n_r]*P[:b]*r/(2*Va)) + P[:C_n_delta_a]*delta_a + P[:C_n_delta_r]*delta_r);
end

@inline @inbounds function testUAV(du,u,p,t)

    Cphi = cos(phi);
    Ctheta = cos(theta);
    Cpsi = cos(psi);

    Sphi = sin(phi);
    Stheta = sin(theta);
    Spsi = sin(psi);

    Tphi = tan(phi);
    Ttheta = tan(theta);
    Tpsi = tan(psi);



    pndot = (Ctheta*Cpsi)*u+(Sphi*Stheta*Cpsi-Cphi*Spsi)*v+(Cphi*Stheta*Cpsi+Sphi*Spsi)*w;
    pedot = (Ctheta*Spsi)*u+(Sphi*Stheta*Spsi+Cphi*Cpsi)*v+(Cphi*Stheta*Spsi-Sphi*Cpsi)*w;
    pddot = (-Stheta)*u+(Sphi*Ctheta)*v+(Cphi*Ctheta)*w;
    udot = (r*v-q*w)+fx/P[:mass];
    vdot = (p*w-r*u)+fy/P[:mass];
    wdot = (q*u-p*v)+fz/P[:mass];
    phidot = (1)*p+(Sphi*Ttheta)*q+(Cphi*Ttheta)*r;
    thetadot = (0)*p+(Cphi)*q+(-Sphi)*r;
    psidot = (0)*p+(Sphi/Ctheta)*q+(Cphi/Ctheta)*r;
    pdot = (P[:gamma1]*p*q-P[:gamma2]*q*r) + (P[:gamma3]*ell+P[:gamma4]*n);
    qdot = (P[:gamma5]*p*r-P[:gamma6]*(p^2-r^2)) + (m/P[:Jy]);
    rdot = (P[:gamma7]*p*q-P[:gamma1]*q*r) + (P[:gamma4]*ell+P[:gamma8]*n);

    du[1] = pndot
    du[2] = pedot
    du[3] = pddot
    du[4] = udot
    du[5] = vdot
    du[6] = wdot
    du[7] = phidot
    du[8] = thetadot
    du[9] = psidot
    du[10] = pdot
    du[11] = qdot
    du[12] = rdot

end

ds = ContinuousDynamicalSystem(loop, rand(3), [10.0, 28.0, 8/3], loop_jac)
data = trajectory(ds, 10) #this returns a dataset

PyPlot.plot3D(data[:,1],data[:,2],data[:,3])

PyPlot.figure()
PyPlot.plot(data[:,1],data[:,2])

integ_1 = integrator(ds, rand(3))
for i=1:100
    # i = 1

    test2 = step!(integ_1,0.1,true)
    PyPlot.plot(integ_1.u[1],integ_1.u[2],"r.")
    PyPlot.pause(0.1)
end

set_parameter!(integ_1,[rand()*10,rand()*20,rand()*2])

for i=1:100
    # i = 1

    test2 = step!(integ_1,0.1,true)
    PyPlot.plot(integ_1.u[1],integ_1.u[2],"r.")
    PyPlot.pause(0.1)
end
