fileLoc,_ = splitdir(@__FILE__)

using Snopt
using Dierckx
using Roots
using PyPlot
PyPlot.close("all")

include("$fileLoc/STOL_dynamics_delta_e.jl")

function objcon(x)
    args = (ds,h_set,dt,PW,optimize,P)
    # println("1")
    return obj(x,args...)
end

function obj(delta_e,ds,h_set,dt,PW,optimize,P)

    #Reinitizilze the integrator
    integ_1 = integrator(ds)

    errorsum = 0
    savestates = zeros(length(delta_e),7)
    Va = zeros(delta_e)
    C_minheight = zeros(delta_e)
    C_vel = zeros(delta_e)
    CL = zeros(delta_e)
    for i = 1:length(delta_e)

        integ_1.p[:delta_e] = delta_e[i]*pi/180
        integ_1.p[:delta_t] = 0.5 #thrust
        try
            integ_1.p[:thrust_weight] = PW/Va[i-1]
        catch
            integ_1.p[:thrust_weight] = PW/integ_1.p[:Va0]
        end

        step!(integ_1,dt, true)
        Va[i] = sqrt(integ_1.u[3]^2+integ_1.u[4]^2)
        thrust = integ_1.p[:thrust_weight]*P[:gravity]*P[:mass]
        # prop_thrust = 0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*thrust)^2-Va[i]^2)
        errorsum += (h_set-(-integ_1.u[2]))
        savestates[i,1:6] = integ_1.u
        savestates[i,7] = thrust

        alpha = atan2(integ_1.u[3],integ_1.u[4])
        sigma_alpha = (1 + exp(-P[:M]*(alpha-P[:alpha0])) + exp(P[:M]*(alpha+P[:alpha0])))/((1+exp(-P[:M]*(alpha-P[:alpha0])))*(1+exp(P[:M]*(alpha+P[:alpha0]))));
        CL[i] = (1-sigma_alpha)*(P[:C_L_0] + P[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha));

        c_vel = integ_1.p[:Va0]-Va[i]

        k = P[:ground_threshold]
        x0 = pi*e/k
        C_vel[i] = c_vel/(1+e^(-k*(-integ_1.u[2]-x0))) # set to zero if on ground, but smoothly
        # println("2")
    end


    J = errorsum/100
    C_minheight = 0.0-(-savestates[:,2])
    C = [C_vel/10;C_minheight]
    # println("3")
    # println(sqrt(integ_1.u[3]^2+integ_1.u[4]^2))

    if optimize
        return J,C,false
    else
        return errorsum,savestates,C_vel,C_minheight,CL
    end
end


# x0  = [P[:pn0],P[:pe0],P[:pd0],P[:u0],P[:v0],P[:w0],P[:phi0],P[:theta0],P[:psi0],P[:p0],P[:q0],P[:r0]]
# x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0]]
x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0],P[:theta0],P[:q0]]

# Set up ContinuousDynamicalSystem
ds = ContinuousDynamicalSystem(testUAV!, x0, P)
# integ_1 = integrator(ds)
# step!(integ_1,dt, true)
# println(integ_1.u)

delta_e = [-64.1588, -85.0, 84.6705, -67.8829, 40.7014, -15.1621, -43.0495, 2.08015, -42.6057, 4.12959, -55.8678, 85.0, -26.2024, -9.43524, -75.5494, 32.4107, -19.205, 17.8462, -85.0, -81.844]
h_set = 10.0
dt = .10
time_sim = collect(dt:dt:dt*length(delta_e))
x0 = delta_e
lb = ones(delta_e)*-85
ub = ones(delta_e)*85
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
options["Major optimality tolerance"] = 1e-5 #Should be scaled so it is optimal with a solid 2 digits
options["Minor optimality  tolerance"] = 1e-6
options["Major feasibility tolerance"] = 1e-5
options["Minor feasibility tolerance"] = 1e-6
options["Minor print level"] = 5
options["Print frequency"] = 100
options["Scale option"] = 1
options["Scale tolerance"] = .95

# PW_array = [.3,.5,.7,.8,.9,.999]#linspace(.459,.46,10)
PW_array = linspace(10,20,5)
PW = []
optimize = []
errorsum = zeros(PW_array)
mean_thrust2weight = zeros(PW_array)
total_energy = zeros(PW_array)
final_dist = zeros(PW_array)

function zero_height(x)
    return  zero_height(x,args...)
end

function zero_height(x,alt_spl,h_set)
    h = alt_spl(x)
    return h-h_set
end

for i = 1:length(PW_array)
    # i = 1

    optimize = true
    PW = PW_array[i]
    # println("5")
    xopt, fopt, optinfo = snopt(objcon, x0, lb, ub, options;printfile = "$(fileLoc)/snopt-print$PW.out", sumfile = "$(fileLoc)/snopt-summary.out") #!!!!!! ARGS are in objcon wrapper function !!!!!!!!#

    println("xopt: $(xopt)")
    println("####################################################")
    println(optinfo)
    println("####################################################")

    optimize = false

    delta_e = xopt
    errorsum[i],savestates,C_vel,C_minheight,CL = objcon(delta_e)
    mean_thrust2weight[i] = mean(savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]))
    Va = sqrt.(savestates[:,3].^2+savestates[:,4].^2)

    figname = "d_CL"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],CL,".-",label = "$(round.(PW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("CL")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)

    figname = "d_Va"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],Va,".-",label = "$(round.(PW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("Va (m/s)")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)
    # PyPlot.ylim(0,32)

    figname = "pn_pd"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],-savestates[:,2],".-",label = "$(round.(PW,4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("-pd (m)")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)
    # PyPlot.xlim(0,70)
    # PyPlot.ylim(0,20)

    figname = "t_T"
    PyPlot.figure(figname)
    PyPlot.plot(time_sim,savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]),".-",label = "$(round.(PW,3))")
    # PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Thrust/Weight")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)

    power = Va.*savestates[:,end]

    figname = "power"
    PyPlot.figure(figname)
    PyPlot.plot(time_sim,power,".-",label = "$(round.(PW,4))")
    # PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Power")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)

    alt_spl = Dierckx.Spline1D(savestates[:,1],-savestates[:,2])
    va_spl = Dierckx.Spline1D(savestates[:,1],Va)
    T_spl = Dierckx.Spline1D(savestates[:,1],savestates[:,end])
    t_spl = Dierckx.Spline1D(savestates[:,1],time_sim)
    P_spl = Dierckx.Spline1D(time_sim,power)

    args = (alt_spl,h_set)
    dist_at_height = fzero(zero_height,0,maximum(savestates[:,1]))
    time_at_height = t_spl(dist_at_height)
    total_energy[i] = integrate(P_spl,0,time_at_height)
    final_dist[i] = dist_at_height

end

figname = "errorsum"
PyPlot.figure(figname)
PyPlot.plot(PW_array,errorsum,".-",label = "$(round.(PW_array,3))")
PyPlot.xlabel("Power/Weight")
PyPlot.ylabel("errorsum")
PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)

figname = "total_energy"
PyPlot.figure(figname)
PyPlot.plot(PW_array,total_energy/(P[:mass]*P[:gravity]),".-",label = "$(round.(PW_array,3))")
PyPlot.xlabel("Power/Weight")
PyPlot.ylabel("Energy/Weight (J/N)")
PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)

figname = "average_slope"
PyPlot.figure(figname)
PyPlot.plot(PW_array,h_set./final_dist,".-",label = "$(round.(PW_array,3))")
PyPlot.xlabel("Power/Weight")
PyPlot.ylabel("Average Slope (H/D)")
PyPlot.savefig("./figures/dynamic_opt/$figname.png",transparent = true)
