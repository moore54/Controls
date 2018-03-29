using DynamicalSystems
# PyPlot.close("all")

include("./STOL_params.jl")


function testUAV!(dinputs,inputs,par,time)

    pn    = inputs[1];
    pd    = inputs[2];
    u     = inputs[3];
    w     = inputs[4];

    theta = par[:theta]
    delta_t = par[:delta_t]

    Va = sqrt(u^2+w^2);
    # println(Va)
    alpha = atan2(w,u);#-par[:alpha0];

    if (Va == 0)
        Va = par[:Va0];
    end

    if !isfinite(alpha)
        alpha = 0;
    end

    #-------- FORCES --------#

   C_D_alpha = par[:C_D_p] + ((par[:C_L_0] + par[:C_L_alpha]*alpha)^2)/(pi*par[:e]*par[:AR]);
   sigma_alpha = (1 + exp(-par[:M]*(alpha-par[:alpha0])) + exp(par[:M]*(alpha+par[:alpha0])))/((1+exp(-par[:M]*(alpha-par[:alpha0])))*(1+exp(par[:M]*(alpha+par[:alpha0]))));
   C_L_alpha = (1-sigma_alpha)*(par[:C_L_0] + par[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));

   C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
   C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);

    # compute external forces and torques on aircraft
    prop_thrust = 0.5*par[:rho]*par[:S_prop]*par[:C_prop]*((par[:k_motor]*delta_t)^2-Va^2)
    fx =  -par[:mass]*par[:gravity]*sin(theta)+0.5*par[:rho]*Va^2*par[:S_wing]*(C_X_alpha)+prop_thrust;
    fz =  par[:mass]*par[:gravity]*cos(theta)+0.5*par[:rho]*Va^2*par[:S_wing]*(C_Z_alpha);

    #-------- DYNAMICS --------#
    Ctheta = cos(theta);
    Stheta = sin(theta);

    pndot = (Ctheta)*u+(Stheta)*w;
    pddot = (-Stheta)*u+(Ctheta)*w;
    udot = fx/par[:mass];
    wdot = fz/par[:mass];

    dinputs[1] = pndot
    dinputs[2] = pddot
    dinputs[3] = udot
    dinputs[4] = wdot

end
