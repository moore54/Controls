fileLoc,_ = splitdir(@__FILE__)

using Snopt
using PyPlot
PyPlot.close("all")

include("$fileLoc/STOL_dynamics.jl")
function objcon(x)
    args = ()
    return obj(x,args...)
end



# x0  = [P[:pn0],P[:pe0],P[:pd0],P[:u0],P[:v0],P[:w0],P[:phi0],P[:theta0],P[:psi0],P[:p0],P[:q0],P[:r0]]
x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0]]

# Set up ContinuousDynamicalSystem
ds = ContinuousDynamicalSystem(testUAV!, x0, P)


theta = [1,2,3,4,5,6,-7,8,9,10.0]
function obj(theta,ds,h_set)

    #Reinitizilze the integrator
    integ_1 = integrator(ds)

    errorsum = 0
    savestates = zeros(length(theta),4)
    for i = 1:length(theta)
        integ_1.p[:theta] = theta[i]*pi/180
        integ_1.p[:delta_e] = -0.0*pi/180
        integ_1.p[:delta_a] = 0.0*pi/180
        integ_1.p[:delta_r] = 0.0*pi/180
        integ_1.p[:delta_t] = 0.5
        step!(integ_1,.5, true)
        println((-integ_1.u[2]))
        errorsum += h_set-(-integ_1.u[2])
        savestates[i,:] = integ_1.u
    end

    return errorsum,savestates
end


@time for i = 1:10
    theta = theta+1
    errorsum,savestates = obj(theta,ds,10.0)
    # println(errorsum)
    PyPlot.scatter3D(savestates[:,1],savestates[:,3],savestates[:,2])
    PyPlot.xlabel("pn")
    PyPlot.ylabel("u")
    PyPlot.zlabel("pd")

end
