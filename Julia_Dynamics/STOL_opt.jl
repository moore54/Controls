fileLoc,_ = splitdir(@__FILE__)

using Snopt
using PyPlot
PyPlot.close("all")

include("$fileLoc/STOL_dynamics_delta_e.jl")

function objcon(x)
    args = (ds,h_set,dt,TW,optimize,P)
    # println("1")
    return obj(x,args...)
end

function obj(delta_e,ds,h_set,dt,TW,optimize,P)

    #Reinitizilze the integrator
    integ_1 = integrator(ds)

    errorsum = 0
    savestates = zeros(length(delta_e),7)
    Va = zeros(delta_e)
    C_minheight = zeros(delta_e)
    C_vel = zeros(delta_e)
    for i = 1:length(delta_e)

        integ_1.p[:delta_e] = delta_e[i]*pi/180
        integ_1.p[:delta_t] = 0.5#thrust
        integ_1.p[:thrust_weight] = TW

        step!(integ_1,dt, true)
        Va[i] = sqrt(integ_1.u[3]^2+integ_1.u[4]^2)
        thrust = P[:thrust_weight]*P[:gravity]*P[:mass]
        # prop_thrust = 0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*thrust)^2-Va[i]^2)
        errorsum += (h_set-(-integ_1.u[2]))
        savestates[i,1:6] = integ_1.u
        savestates[i,7] = thrust

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
        return errorsum,savestates,C_vel,C_minheight
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

delta_e = [-80.0, -23.6311, -15.4772, 0.0, 0.0, -9.52617, -37.6217, -17.1981, -20.3374, -31.1988, -13.0533, -25.3677, -28.7171, -10.6457, -34.9761, -11.4172, -38.6286, 0.0, -80.0, -80.0]
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

TW_array = [.3,.5,.7,.8,.9,.999]#linspace(.459,.46,10)
TW = []
optimize = []
errorsum = zeros(TW_array)
mean_thrust2weight = zeros(TW_array)
total_energy = zeros(TW_array)
for i = 1:length(TW_array)

    optimize = true
    TW = TW_array[i]
    # println("5")
    xopt, fopt, optinfo = snopt(objcon, x0, lb, ub, options;printfile = "$(fileLoc)/snopt-print$TW.out", sumfile = "$(fileLoc)/snopt-summary.out") #!!!!!! ARGS are in objcon wrapper function !!!!!!!!#

    println("xopt: $(xopt)")
    println("####################################################")
    println(optinfo)
    println("####################################################")

    optimize = false

    delta_e = xopt
    errorsum[i],savestates,C_vel,C_minheight = objcon(delta_e)
    mean_thrust2weight[i] = mean(savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]))
    Va = sqrt.(savestates[:,3].^2+savestates[:,4].^2)

    figname = "d_Va"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],Va,".-",label = "T/weight $(round.(mean_thrust2weight[i],4))")
    PyPlot.xlabel("x-loc (m)")
    PyPlot.ylabel("Va (m/s)")
    PyPlot.legend(loc = "best")
    PyPlot.savefig("./figures/$figname.png",transparent = true)
    # PyPlot.ylim(0,32)

    figname = "pn_pd"
    PyPlot.figure(figname)
    PyPlot.plot(savestates[:,1],-savestates[:,2],".-",label = "T/weight $(round.(mean_thrust2weight[i],4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("pd (m)")
    PyPlot.legend(loc = "best")
    PyPlot.savefig("./figures/$figname.png",transparent = true)
    # PyPlot.xlim(0,70)
    # PyPlot.ylim(0,20)

    figname = "t_T"
    PyPlot.figure(figname)
    PyPlot.plot(time_sim,savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]),".-",label = "Thrust/Weight $(round.(TW,3))")
    # PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Thrust/Weight")
    PyPlot.legend(loc = "best")
    PyPlot.savefig("./figures/$figname.png",transparent = true)

    power = Va.*savestates[:,end]

    figname = "power"
    PyPlot.figure(figname)
    PyPlot.plot(time_sim,power,".-",label = "T/weight $(round.(mean_thrust2weight[i],4))")
    # PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Power")
    PyPlot.legend(loc = "best")
    PyPlot.savefig("./figures/$figname.png",transparent = true)

    total_energy[i] = sum(power.*dt)

end

figname = "errorsum"
PyPlot.figure(figname)
PyPlot.plot(mean_thrust2weight,errorsum,".-",label = "Thrust/Weight $(round.(TW,3))")
# PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("errorsum")
PyPlot.legend(loc = "best")
PyPlot.savefig("./figures/$figname.png",transparent = true)


figname = "total_energy"
PyPlot.figure(figname)
PyPlot.plot(mean_thrust2weight,total_energy,".-",label = "Thrust/Weight $(round.(TW,3))")
# PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("total_energy (J)")
PyPlot.legend(loc = "best")
PyPlot.savefig("./figures/$figname.png",transparent = true)
