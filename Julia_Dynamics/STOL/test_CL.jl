using PyPlot
close("all")
rc("figure", figsize=(6.5, 3.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.89)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))

include("./STOL_dynamics_delta_e.jl")

alpha_array = linspace(-5,30,500).*pi/180
CL = zeros(alpha_array)
CD = zeros(alpha_array)
Fx = zeros(alpha_array)
Fy = zeros(alpha_array)
m_pitch = zeros(alpha_array)
theta = 0.0
Va = 20.0
q = 0.0
delta_e = -00.0*pi/180
phi = 0.0
psi = 0.0
for i = 1:length(alpha_array)
    alpha = alpha_array[i]
    C_D_alpha = P[:C_D_p] + ((P[:C_L_0] + P[:C_L_alpha]*alpha)^2)/(pi*P[:e]*P[:AR]);
    sigma_alpha = (1 + exp(-P[:M]*(alpha-P[:alpha0])) + exp(P[:M]*(alpha+P[:alpha0])))/((1+exp(-P[:M]*(alpha-P[:alpha0])))*(1+exp(P[:M]*(alpha+P[:alpha0]))));
    C_L_alpha = (1-sigma_alpha)*(P[:C_L_0] + P[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));
    CL[i] = C_L_alpha
    CD[i] = C_D_alpha
    C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
    C_X_q_alpha = -P[:C_D_q]*cos(alpha) + P[:C_L_q]*sin(alpha);
    C_X_delta_e_alpha = -P[:C_D_delta_e]*cos(alpha) + P[:C_L_delta_e]*sin(alpha);
    C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);
    C_Z_q_alpha = -P[:C_D_q]*sin(alpha) - P[:C_L_q]*cos(alpha);
    C_Z_delta_e_alpha = -P[:C_D_delta_e]*sin(alpha) - P[:C_L_delta_e]*cos(alpha);

    # compute external forces and torques on aircraft
    thrust = P[:thrust_weight]*P[:gravity]*P[:mass]
    # println(thrust)
    k = P[:ground_threshold]
    x0 = pi*e/k
    grav = P[:gravity]

    Force1 =  -P[:mass]*P[:gravity]*sin(theta)+0.5*P[:rho]*Va^2*P[:S_wing]*(C_X_alpha+C_X_q_alpha*P[:c]/(2*Va)*q+C_X_delta_e_alpha*delta_e)+thrust #0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*delta_t)^2-Va^2);
    # Force2 =  P[:mass]*P[:gravity]*cos(theta)*sin(phi)+0.5*P[:rho]*Va^2*P[:S_wing]*(P[:C_Y_0]+P[:C_Y_beta]*beta+P[:C_Y_p]*P[:b]/(2*Va)*p+P[:C_Y_r]*P[:b]/(2*Va)*r+P[:C_Y_delta_a]*delta_a+P[:C_Y_delta_r]*delta_r);
    Force3 =  P[:mass]*grav*cos(theta)*cos(phi)+0.5*P[:rho]*Va^2*P[:S_wing]*(C_Z_alpha+C_Z_q_alpha*P[:c]/(2*Va)*q+C_Z_delta_e_alpha*delta_e);
    Torque2 = 0.5*P[:rho]*Va^2*P[:S_wing]*P[:c]*(P[:C_m_0] + P[:C_m_alpha]*alpha + (P[:C_m_q]*P[:c]*q/(2*Va)) + P[:C_m_delta_e]*delta_e);
    Fx[i] = Force1
    Fy[i] = Force3
    m_pitch[i] = Torque2
end

LD,idx_LD = findmax(CL./CD)

figname = "L_D"
figure(figname)
plot(alpha_array*180/pi,CL./CD)
plot(alpha_array[idx_LD]*180/pi,LD,"o")
xlabel("AOA (deg)")
ylabel("L/D")
PyPlot.savefig("./figures/test_CL/$figname.png",transparent = true)


figname = "CD"
figure(figname)
plot(alpha_array*180/pi,CD)
plot(alpha_array[idx_LD]*180/pi,CD[idx_LD],"o")
xlabel("AOA (deg)")
ylabel("CD")
PyPlot.savefig("./figures/test_CL/$figname.png",transparent = true)
#
# figure("fxfy")
# plot(alpha_array*180/pi,Fx)
# plot(alpha_array*180/pi,Fy)
# xlabel("AOA (deg)")
# ylabel("Fx/Fy")
#
# figure("fxfy2")
# plot(alpha_array*180/pi,Fy./Fx)
# xlabel("AOA (deg)")
# ylabel("Fy/Fx")

# fxy,idx_fxy = findmin(sqrt.(Fx.^2+Fy.^2))

# figname = "norm_fxfy"
# figure("norm_fxfy")
# plot(alpha_array*180/pi,sqrt.(Fx.^2+Fy.^2))
# # plot(alpha_array[idx_fxy]*180/pi,sqrt.(Fx[idx_fxy].^2+Fy[idx_fxy].^2),"o")
# xlabel("AOA (deg)")
# ylabel("Norm(Fx,Fy)")
# PyPlot.savefig("./figures/test_CL/$figname.png",transparent = true)

clmax,idx_clmax = findmax(CL)

figname = "CL"
figure("CL")
plot(alpha_array*180/pi,CL)
plot(alpha_array[idx_LD]*180/pi,CL[idx_LD],"o",label = "Max L/D")
plot(alpha_array[idx_clmax]*180/pi,CL[idx_clmax],"o",label = "CLmax = $(round(CL[idx_clmax],3)) at $(round(alpha_array[idx_clmax]*180/pi,3))deg")
legend(loc = "best")
xlabel("AOA (deg)")
ylabel("CL")
PyPlot.savefig("./figures/test_CL/$figname.png",transparent = true)
