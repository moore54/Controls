fileLoc,_ = splitdir(@__FILE__)

using Snopt
using Dierckx
using Roots
using PyPlot
PyPlot.close("all")
rc("figure", figsize=(6.5, 4.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.69)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])


include("$fileLoc/STOL_dynamics_delta_e.jl")

function objcon(x)
    args = (ds,h_set,dt,TW,optimize,P)
    # println("1")
    return obj(x,args...)
end

function obj(x,ds,h_set,dt,TW,optimize,P)

    #Reinitizilze the integrator
    integ_1 = integrator(ds)
    dt = x[1]
    delta_e = x[2:21]


    errorsum = 0
    savestates = zeros(length(delta_e),7)
    Va = zeros(delta_e)
    C_minheight = zeros(delta_e)
    C_vel = zeros(delta_e)
    C_turnback = zeros(delta_e)
    CL = zeros(delta_e)
    alpha_save = zeros(delta_e)
    for i = 1:length(delta_e)
        integ_1.p[:delta_e] = delta_e[i]*pi/180
        # integ_1.p[:delta_t] = 0.5 #thrust

        if i>1
            integ_1.p[:thrust_weight] = TW #PW/Va[i-1]
        else
            integ_1.p[:thrust_weight] = TW #PW/integ_1.p[:Va0]
        end

        step!(integ_1,dt, true)
        Va[i] = sqrt(integ_1.u[3]^2+integ_1.u[4]^2)
        thrust = integ_1.p[:thrust_weight]*P[:gravity]*P[:mass]
        # prop_thrust = 0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*thrust)^2-Va[i]^2)
        errorsum += (h_set-(-integ_1.u[2]))
        savestates[i,1:6] = integ_1.u
        savestates[i,7] = thrust

        alpha = atan2(integ_1.u[4],integ_1.u[3])
        sigma_alpha = (1 + exp(-P[:M]*(alpha-P[:alpha0])) + exp(P[:M]*(alpha+P[:alpha0])))/((1+exp(-P[:M]*(alpha-P[:alpha0])))*(1+exp(P[:M]*(alpha+P[:alpha0]))));
        CL[i] = (1-sigma_alpha)*(P[:C_L_0] + P[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));
        alpha_save[i] = alpha*180/pi
        c_vel = integ_1.p[:Va0]-Va[i]

        k = P[:ground_threshold]
        x0 = pi*e/k
        C_vel[i] = c_vel/(1+e^(-k*(-integ_1.u[2]-x0))) # set to zero if on ground, but smoothly
        if i>1
            C_turnback[i] = savestates[i-1,1]-savestates[i,1]
        else
            C_turnback[i] = -savestates[i,1]
        end
        # println("2")
    end


    J = savestates[end,1]/10#errorsum/100
    C_minheight = 0.0-(-savestates[:,2])
    C_maxheight = h_set - -savestates[end,2]
    C_endslope = -savestates[end,2]-(-savestates[end-1,2])
    C = [C_maxheight/h_set;C_vel/10;C_minheight;C_turnback/10;C_endslope]
    # println("3")
    # println(sqrt(integ_1.u[3]^2+integ_1.u[4]^2))

    if optimize
        return J,C,false
    else
        return errorsum,savestates,C_maxheight,C_vel,C_minheight,CL,alpha_save
    end
end


# x0  = [P[:pn0],P[:pe0],P[:pd0],P[:u0],P[:v0],P[:w0],P[:phi0],P[:theta0],P[:psi0],P[:p0],P[:q0],P[:r0]]
# x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0]]
x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0],P[:theta0],P[:q0]]

# Set up ContinuousDynamicalSystem
ds = ContinuousDynamicalSystem(testUAV!, x0, P)
#Reinitizilze the integrator

delta_e = [-85.0, 1.16896, 52.1414, -32.6719, -2.07146, -15.4012, -9.35701, -11.8722, -10.8659, -11.2834, -11.2427, -10.1897, -12.6121, -9.07863, -13.8728, -6.67249, -19.487, 9.64414, -85.0, -85.0]
# PW_var = [5.0,50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 24.1828, 11.2694, 10.3522, 39.4638, 9.19733, 15.5779]
# delta_e = ones(20)*.5
h_set = 50.0
h_M = 10.0/13.5
momentum = P[:Va0]*P[:mass]
dt = .10
TW = 0.8



time_sim = collect(dt:dt:dt*length(delta_e))

lb_de = ones(delta_e)*-85
ub_de = ones(delta_e)*85



# x0 = [-85.0, 38.9728, -6.79026, 4.81065, 1.10381, -20.7914, -9.42534, -10.6752, -9.40374, -11.4258, -9.76186, -6.14866, -10.7634, -12.9517, -16.3083, -8.30803, -13.2473, 20.5823, -79.232, -85.0, 1.86897, 0.465348, 7.67862, 10.0, 10.0, 10.0, 10.0, 10.0, 9.91234, 9.15937, 5.05204, 0.568162, 0.335101, 2.02494, 9.59265, 10.0, 10.0, 10.0, 10.0, 10.0]
# x0 = [-85.0, 38.5752, -5.39709, 3.44675, 2.27539, -21.6652, -7.49748, -10.5739, -13.0014, -7.24735, -14.1903, -8.88464, -11.6929, -9.35434, -12.2483, -7.6421, -20.109, 16.0491, -85.0, -85.0, 1.93352, 2.3695, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0]
# x0 = [.1,-84.5883, 71.655, -83.1494, 71.6458, -39.0269, 3.09616, -15.8954, -12.7499, -4.09139, -18.0408, -8.04835, -3.37305, -26.406, 7.02331, -23.5117, -7.31014, 7.87563, -23.9456, 60.2123, -85.0, 0.381699, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
# x0 = [0.424256, -12.2756, 3.33995, -13.0701, -4.2479, -7.87497, -6.01732, -6.61042, -6.17238, -6.28871, -6.22684, -6.33089, -6.4083, -6.54489, -6.6766, -6.76479, -6.96485, -6.74469, -7.44002, -5.89143, 3.78196]
x0 = [0.204136, -24.3492, 24.3332, -13.7759, -8.60549, -14.6692, -4.74278, -10.1546, -10.2219, -10.4949, -10.3573, -10.3058, -9.95381, -3.39003, 1.13709, -1.46905, -7.42902, -0.221764, 2.94481, -3.73116, 4.11884]
# x0 = [.1;delta_e]
lb = [1E-3;lb_de]
ub = [5.0;ub_de]
# println("4")
# ----- Define Optimizer Options ----- #
options = Dict{String, Any}()
options["Derivative level"] = 0
# options["Function precision"] = 1.00E-4
# options["Difference interval"] = 1e-4
# options["Central difference interval"] = 1e-4
options["Iterations limit"] = 1e8
options["Major iterations limit"] = 1000
options["Minor iterations limit"]= 1e8
options["Major optimality tolerance"] = 1e-4 #Should be scaled so it is optimal with a solid 2 digits
options["Minor optimality  tolerance"] = 1e-6
options["Major feasibility tolerance"] = 1e-5
options["Minor feasibility tolerance"] = 1e-6
options["Minor print level"] = 5
options["Print frequency"] = 100
options["Scale option"] = 1
options["Scale tolerance"] = .95

# PW_array = [.3,.5,.7,.8,.9,.999]#linspace(.459,.46,10)
TW_array = [.7,.9,1.0]
# args = []
args2 = []
optimize = []
savestates = []
errorsum = zeros(TW_array)
mean_thrust2weight = zeros(TW_array)
total_energy = zeros(TW_array)
final_dist = zeros(TW_array)

function zero_height(x)
    return  zero_height(x,args2...)
end

function zero_height(x,alt_spl,h_set)
    h = alt_spl(x)
    return h-h_set
end
testrun = false
for i = 1:length(TW_array)
# i = 2
    TW = TW_array[i]
    # P[:mass] = Mass_array[i]
    # Mass = Mass_array[i]
    # P[:Va0] = momentum/Mass
    # println(P[:Va0])
    # h_set = h_M*Mass
    if !testrun
        optimize = true


        xopt, fopt, optinfo = snopt(objcon, x0, lb, ub, options;printfile = "$(fileLoc)/snopt-print$(TW).out", sumfile = "$(fileLoc)/snopt-summary.out") #!!!!!! ARGS are in objcon wrapper function !!!!!!!!#

        println("xopt: $(xopt)")
        println("####################################################")
        println(optinfo)
        println("####################################################")
    else
        xopt = x0
    end
    optimize = false


    errorsum[i],savestates,C_maxheight,C_vel,C_minheight,CL,alpha_save = objcon(xopt)
    mean_thrust2weight[i] = mean(savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]))
    Va = sqrt.(savestates[:,3].^2+savestates[:,4].^2)

    figname = "d_CL"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],CL,".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("CL")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    figname = "d_q"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],savestates[:,6],".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("q")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    figname = "d_theta"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],savestates[:,5]*180/pi,".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("theta (deg)")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    figname = "d_AOA"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],alpha_save,".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("AOA (deg)")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    figname = "d_Va"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],Va,".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("Va (m/s)")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)
    # PyPlot.ylim(0,32)

    figname = "pn_pd"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],-savestates[:,2],".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("-pd (m)")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    figname = "pn_pd_norm"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1]/maximum(savestates[:,1]),-savestates[:,2]/maximum(-savestates[:,2]),".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn/max pn")
    PyPlot.ylabel("-pd/max -pd")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)
    # PyPlot.xlim(0,70)
    # PyPlot.ylim(0,20)

    figname = "t_T"
    PyPlot.figure(figname)
    PyPlot.plot(time_sim,savestates[:,end]./(P[:mass]*P[:gravity]),".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Thrust/Weight")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    power = Va.*savestates[:,end]

    figname = "d_power"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],power./(P[:mass]*P[:gravity]),".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("Power/weight")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    figname = "d_delta_e"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],xopt[2:21],".-",label = "$(round.(TW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("Elevator Deflection (deg)")
    legend(loc="center left", title = "Thrust/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

    alt_spl = Dierckx.Spline1D(savestates[:,1],-savestates[:,2])
    va_spl = Dierckx.Spline1D(savestates[:,1],Va)
    T_spl = Dierckx.Spline1D(savestates[:,1],savestates[:,end])
    t_spl = Dierckx.Spline1D(savestates[:,1],time_sim)
    P_spl = Dierckx.Spline1D(time_sim,power)

    args2 = (alt_spl,h_set)
    dist_at_height = 0.0
    try
        dist_at_height = fzero(zero_height,0,maximum(savestates[:,1]))
    catch
        dist_at_height = maximum(savestates[:,1])
        if -savestates[end,2]<h_set*.999
            warn("Max height was $(savestates[end,2]), should be $h_set")
        end
    end
    time_at_height = t_spl(dist_at_height)
    total_energy[i] = integrate(P_spl,0,time_at_height)
    final_dist[i] = dist_at_height

end

figname = "errorsum"
PyPlot.figure(figname)
PyPlot.plot(TW_array,errorsum,".-",label = "$(round.(TW_array,3))")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("errorsum")
PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

figname = "total_energy"
PyPlot.figure(figname)
PyPlot.plot(TW_array,total_energy/(P[:mass]*P[:gravity]),".-",label = "$(round.(TW_array,3))")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("Energy/Weight (J/N)")
PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)

figname = "average_slope"
PyPlot.figure(figname)
PyPlot.plot(TW_array,h_set./final_dist,".-",label = "$(round.(TW_array,3))")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("Average Slope (H/D)")
PyPlot.savefig("./figures/H_50_0slope/$figname.png",transparent = true)
