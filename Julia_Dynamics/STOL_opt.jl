fileLoc,_ = splitdir(@__FILE__)

using Snopt
using PyPlot
PyPlot.close("all")

include("$fileLoc/STOL_dynamics_delta_e.jl")

function objcon(x)
    args = (ds,h_set,dt,PW,optimize,P)
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
    for i = 1:length(delta_e)

        integ_1.p[:delta_e] = delta_e[i]*pi/180
        # integ_1.p[:delta_t] = thrust
        integ_1.p[:power_weight] = PW

        step!(integ_1,dt, true)
        Va[i] = sqrt(integ_1.u[3]^2+integ_1.u[4]^2)
        thrust = par[:power_weight]*par[:gravity]*par[:mass]/Va[i]
        # prop_thrust = 0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*thrust)^2-Va[i]^2)
        errorsum += (h_set-(-integ_1.u[2]))
        savestates[i,1:6] = integ_1.u
        savestates[i,7] = thrust

        if -integ_1.u[2]>integ_1.p[:ground_threshold]
            C_vel[i] = integ_1.p[:Va0]-Va[i]
        else
            C_vel[i] = 0.0
        end

    end

    if optimize
        J = errorsum
        C_vel = integ_1.p[:Va0]-Va
        C_minheight = 0.0-(-savestates[:,2])
        C = [C_vel;C_minheight]

        # println(sqrt(integ_1.u[3]^2+integ_1.u[4]^2))
        return J,C,false
    else
        return errorsum,savestates
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

delta_e = ones(20)*0
h_set = 10.0
dt = .2
time_sim = collect(dt:dt:dt*length(delta_e))
x0 = delta_e
lb = ones(delta_e)*-180
ub = ones(delta_e)*180

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

PW_array = [.6]#linspace(.459,.46,10)
PW = []
optimize = []
errorsum = zeros(PW_array)
mean_thrust2weight = zeros(PW_array)
total_energy = zeros(PW_array)
for i = 1#:length(thrust_array)
    optimize = true
    PW = PW_array[i]
    xopt, fopt, optinfo = snopt(objcon, x0, lb, ub, options;printfile = "$(fileLoc)/snopt-print$h_set.out", sumfile = "$(fileLoc)/snopt-summary.out") #!!!!!! ARGS are in objcon wrapper function !!!!!!!!#

    println("xopt: $(xopt)")
    println("####################################################")
    println(optinfo)
    println("####################################################")

    optimize = false

    delta_e = xopt
    errorsum[i],savestates = objcon(delta_e)
    mean_thrust2weight[i] = mean(savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]))
    Va = sqrt.(savestates[:,3].^2+savestates[:,4].^2)

    PyPlot.figure("d_Va")
    PyPlot.plot(savestates[:,1],Va,".-",label = "T/weight $(round.(mean_thrust2weight[i],4))")
    PyPlot.xlabel("x-loc (m)")
    PyPlot.ylabel("Va (m/s)")
    PyPlot.legend(loc = "best")
    # PyPlot.ylim(0,32)

    PyPlot.figure("pn_pd")
    PyPlot.plot(savestates[:,1],-savestates[:,2],".-",label = "T/weight $(round.(mean_thrust2weight[i],4))")
    PyPlot.xlabel("pn (m)")
    PyPlot.ylabel("pd (m)")
    PyPlot.legend(loc = "best")
    # PyPlot.xlim(0,70)
    # PyPlot.ylim(0,20)

    PyPlot.figure("t_T")
    PyPlot.plot(time_sim,savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]),".-",label = "delta_t $(round.(thrust,3))")
    # PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Thrust/Weight")
    PyPlot.legend(loc = "best")

    power = Va.*savestates[:,end]

    PyPlot.figure("power")
    PyPlot.plot(time_sim,power,".-",label = "T/weight $(round.(mean_thrust2weight[i],4))")
    # PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
    PyPlot.xlabel("time (s)")
    PyPlot.ylabel("Power")
    PyPlot.legend(loc = "best")

    total_energy[i] = sum(power.*dt)

end

PyPlot.figure("errorsum")
PyPlot.plot(mean_thrust2weight,errorsum,".-",label = "Power/Weight $(round.(PW,3))")
# PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("errorsum")
PyPlot.legend(loc = "best")


PyPlot.figure("total_energy")
PyPlot.plot(mean_thrust2weight,total_energy,".-",label = "delta_t $(round.(thrust,3))")
# PyPlot.plot(time_sim,ones(time_sim)*P[:mass]*P[:gravity],"k-")
PyPlot.xlabel("Thrust/Weight")
PyPlot.ylabel("total_energy (J)")
PyPlot.legend(loc = "best")
