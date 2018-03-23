using DynamicalSystems, StaticArrays, PyPlot
PyPlot.close("all")

include("./STOL_params.jl")


function testUAV!(dinputs,inputs,par,time)

    pn    = inputs[1];
    pe    = inputs[2];
    pd    = inputs[3];
    u     = inputs[4];
    v     = inputs[5];
    w     = inputs[6];
    phi   = inputs[7];
    theta = inputs[8];
    psi   = inputs[9];
    p     = inputs[10];
    q     = inputs[11];
    r     = inputs[12];

    delta_e = par[:delta_e]
    delta_a = par[:delta_a]
    delta_r = par[:delta_r]
    delta_t = par[:delta_t]

    # compute air data in body frame
    ur = u#-W_BODY[1];
    vr = v#-W_BODY[2];
    wr = w#-W_BODY[3];

    Va = sqrt(ur^2+vr^2+wr^2);
    alpha = atan2(wr,ur);#-par[:alpha0];
    beta = asin(vr/Va);

    if (Va == 0)
        Va = par[:Va0];
    end

    if !isfinite(alpha)
        alpha = 0;
    end
    if !isfinite(beta)
        beta = 0;
    end

    #-------- FORCES --------#

   C_D_alpha = par[:C_D_p] + ((par[:C_L_0] + par[:C_L_alpha]*alpha)^2)/(pi*par[:e]*par[:AR]);
   sigma_alpha = (1 + exp(-par[:M]*(alpha-par[:alpha0])) + exp(par[:M]*(alpha+par[:alpha0])))/((1+exp(-par[:M]*(alpha-par[:alpha0])))*(1+exp(par[:M]*(alpha+par[:alpha0]))));
   C_L_alpha = (1-sigma_alpha)*(par[:C_L_0] + par[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));

   C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
   C_X_q_alpha = -par[:C_D_q]*cos(alpha) + par[:C_L_q]*sin(alpha);
   C_X_delta_e_alpha = -par[:C_D_delta_e]*cos(alpha) + par[:C_L_delta_e]*sin(alpha);
   C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);
   C_Z_q_alpha = -par[:C_D_q]*sin(alpha) - par[:C_L_q]*cos(alpha);
   C_Z_delta_e_alpha = -par[:C_D_delta_e]*sin(alpha) - par[:C_L_delta_e]*cos(alpha);

    # compute external forces and torques on aircraft
    Force1 =  -par[:mass]*par[:gravity]*sin(theta)+0.5*par[:rho]*Va^2*par[:S_wing]*(C_X_alpha+C_X_q_alpha*par[:c]/(2*Va)*q+C_X_delta_e_alpha*delta_e)+0.5*par[:rho]*par[:S_prop]*par[:C_prop]*((par[:k_motor]*delta_t)^2-Va^2);
    Force2 =  par[:mass]*par[:gravity]*cos(theta)*sin(phi)+0.5*par[:rho]*Va^2*par[:S_wing]*(par[:C_Y_0]+par[:C_Y_beta]*beta+par[:C_Y_p]*par[:b]/(2*Va)*p+par[:C_Y_r]*par[:b]/(2*Va)*r+par[:C_Y_delta_a]*delta_a+par[:C_Y_delta_r]*delta_r);
    Force3 =  par[:mass]*par[:gravity]*cos(theta)*cos(phi)+0.5*par[:rho]*Va^2*par[:S_wing]*(C_Z_alpha+C_Z_q_alpha*par[:c]/(2*Va)*q+C_Z_delta_e_alpha*delta_e);

    Torque1 = 0.5*par[:rho]*Va^2*par[:S_wing]*par[:b]*(par[:C_ell_0] + par[:C_ell_beta]*beta + (par[:C_ell_p]*par[:b]*p/(2*Va)) + (par[:C_ell_r]*par[:b]*r/(2*Va)) + par[:C_ell_delta_a]*delta_a + par[:C_ell_delta_r]*delta_r)- par[:k_T_P]*(par[:k_Omega]*delta_t)^2;
    Torque2 = 0.5*par[:rho]*Va^2*par[:S_wing]*par[:c]*(par[:C_m_0] + par[:C_m_alpha]*alpha + (par[:C_m_q]*par[:c]*q/(2*Va)) + par[:C_m_delta_e]*delta_e);
    Torque3 = 0.5*par[:rho]*Va^2*par[:S_wing]*par[:b]*(par[:C_n_0] + par[:C_n_beta]*beta + (par[:C_n_p]*par[:b]*p/(2*Va)) + (par[:C_n_r]*par[:b]*r/(2*Va)) + par[:C_n_delta_a]*delta_a + par[:C_n_delta_r]*delta_r);


    #-------- DYNAMICS --------#
    fx    = Force1
    fy    = Force2
    fz    = Force3
    ell   = Torque1
    m     = Torque2
    n     = Torque3

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
    udot = (r*v-q*w)+fx/par[:mass];
    vdot = (p*w-r*u)+fy/par[:mass];
    wdot = (q*u-p*v)+fz/par[:mass];
    phidot = (1)*p+(Sphi*Ttheta)*q+(Cphi*Ttheta)*r;
    thetadot = (0)*p+(Cphi)*q+(-Sphi)*r;
    psidot = (0)*p+(Sphi/Ctheta)*q+(Cphi/Ctheta)*r;
    pdot = (par[:gamma1]*p*q-par[:gamma2]*q*r) + (par[:gamma3]*ell+par[:gamma4]*n);
    qdot = (par[:gamma5]*p*r-par[:gamma6]*(p^2-r^2)) + (m/par[:Jy]);
    rdot = (par[:gamma7]*p*q-par[:gamma1]*q*r) + (par[:gamma4]*ell+par[:gamma8]*n);

    dinputs[1] = pndot
    dinputs[2] = pedot
    dinputs[3] = pddot
    dinputs[4] = udot
    dinputs[5] = vdot
    dinputs[6] = wdot
    dinputs[7] = phidot
    dinputs[8] = thetadot
    dinputs[9] = psidot
    dinputs[10] = pdot
    dinputs[11] = qdot
    dinputs[12] = rdot

end


x0  = [P[:pn0],P[:pe0],P[:pd0],P[:u0],P[:v0],P[:w0],P[:phi0],P[:theta0],P[:psi0],P[:p0],P[:q0],P[:r0]]

println("testUAV")
println("Compile & Run ContinuousDynamicalSystem()")
@time ds = ContinuousDynamicalSystem(testUAV!, x0, P)
println("Run ContinuousDynamicalSystem()")
@time ds = ContinuousDynamicalSystem(testUAV!, x0, P)

println("Compile & Run trajectory()")
@time data = trajectory(ds, 100; dt = 0.01) #this returns a dataset
println("Run trajectory()")
@time data = trajectory(ds, 100; dt = 0.01) #this returns a dataset
#
PyPlot.plot3D(data[:,1],data[:,2],data[:,3])
PyPlot.xlabel("x")
PyPlot.ylabel("y")
PyPlot.zlabel("z")

println("Compile & Run integrator()")
@time integ_1 = integrator(ds)
println("Run integrator()")
@time integ_1 = integrator(ds)

println("Compile & Run step!()")
@time step!(integ_1,0.5, true)
println("Run step!()")
@time step!(integ_1,0.5, true)

for i = 1:10
    integ_1.p[:delta_e] = 0.0*pi/180
    integ_1.p[:delta_a] = 0.0*pi/180
    integ_1.p[:delta_r] = 0.0*pi/180
    integ_1.p[:delta_t] = 0.25
    step!(integ_1,0.5, true)
    PyPlot.scatter3D(integ_1.u[1],integ_1.u[2],integ_1.u[3],"r.")
end
