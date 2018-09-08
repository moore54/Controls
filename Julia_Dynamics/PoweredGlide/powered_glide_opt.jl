fileLoc,_ = splitdir(@__FILE__)
print_iters = 0
using Snopt
using Dierckx
using Roots
using PyPlot
using Gradients
close("all")
rc("figure", figsize=(7.5, 4.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.69)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_colors=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

include("$fileLoc/STOL_dynamics_powered_glide.jl")

function objcon(x)
    args = (ds,h_set,dt,distance,optimize,scale,P)
    # println("1")
    return obj(x,args...)
end

function objgrad(x)
    function func(x)
        j,C,_ = objcon(x)
        return [j;C]
    end
    f, dfdx = Gradients.forwarddiff(func,x)
    lenx = length(x)
    println("x: $x")
    println("f: $(f[1])")
    println("c: $(f[2:end])")
    println("df: $(dfdx[1,:])")
    for i = 1:length(f[2:end])
        println("dc$i: $(dfdx[i+1,:])")
    end
    #function, constraints, dfdx, dcdx, fail
    return f[1],f[2:end],dfdx[1,:],dfdx[2:end,:],false
end

function obj(x,ds,h_set,dt,distance,optimize,scale,P)
    x = x./scale
    #Reinitizilze the integrator
    integ_1 = integrator(ds)
    dt = x[1]
    delta_e = x[2:41]
    TW = x[42:81]

    errorsum = 0
    savestates = zeros(length(delta_e),7)
    Va = zeros(delta_e)
    C_minheight = zeros(delta_e)
    C_vel = zeros(delta_e)
    C_turnback = zeros(delta_e)
    CL = zeros(delta_e)
    CD = zeros(delta_e)
    alpha_save = zeros(delta_e)
    for i = 1:length(delta_e)
        integ_1.p[:delta_e] = delta_e[i]
        # integ_1.p[:delta_t] = 0.5 #thrust

        integ_1.p[:thrust_weight] = TW[i] #PW/Va[i-1] #P[:thrust_weight] #


        step!(integ_1,dt, true)
        Va[i] = sqrt(integ_1.u[3]^2+integ_1.u[4]^2)
        thrust = integ_1.p[:thrust_weight]*P[:gravity]*P[:mass]
        # prop_thrust = 0.5*P[:rho]*P[:S_prop]*P[:C_prop]*((P[:k_motor]*thrust)^2-Va[i]^2)
        errorsum += (h_set-(-integ_1.u[2]))
        savestates[i,1:6] = integ_1.u
        savestates[i,7] = thrust

        alpha = atan2(integ_1.u[4],integ_1.u[3])
        sigma_alpha = (1 + exp(-P[:M]*(alpha-P[:alpha0])) + exp(P[:M]*(alpha+P[:alpha0])))/((1+exp(-P[:M]*(alpha-P[:alpha0])))*(1+exp(P[:M]*(alpha+P[:alpha0]))));
        CL[i] = (1-sigma_alpha)*(P[:C_L_0] + P[:C_L_alpha]*alpha) + sigma_alpha*(2*sign(alpha)*(sin(alpha)^2)*cos(alpha))
        CD[i] = P[:C_D_p] + ((P[:C_L_0] + P[:C_L_alpha]*alpha)^2)/(pi*P[:e]*P[:AR]);
        alpha_save[i] = alpha*180/pi
        c_vel = integ_1.p[:Va0]-Va[i] # velocity must be above stall velocity

        k = P[:ground_threshold]
        x0 = pi*e/k
        C_vel[i] = c_vel/(1+e^(-k*(-integ_1.u[2]-x0))) # set to zero if on ground, but smoothly
        if i>1
            C_turnback[i] = savestates[i-1,1]-savestates[i,1]  #Don't go straight up or down and turn backwards
        else
            C_turnback[i] = -savestates[i,1]
        end
        # println("2")
    end

    power = Va.*savestates[:,end] / 0.75 #convert to electric power
    d_height = savestates[1,2] - savestates[end,2]
    # J = sum(power*dt)/1E4 #savestates[end,1]/10#errorsum/100 # minimum energy expended
    J = (sum(power*dt)-P[:gravity]*P[:mass]*d_height)/1E4 #savestates[end,1]/10# minimize energy spent both propulsive and potential
    C_dist = distance - savestates[end,1] #distance traveled must be greater or equal to distance required
    C_minheight = 0.0-(-savestates[:,2]) # height cannot dip below 0
    C_maxheight = h_set - (maximum(-savestates[end,2]))
    C_endslope = -savestates[end,2]-(-savestates[end-1,2])
    C_power = power/(P[:mass]*P[:gravity])-P[:power_weight]
    Cvel_max = Va-integ_1.p[:Va0]*1.01
    # C = [C_vel/10;C_minheight;C_turnback/10;C_power;Cvel_max]
    C = [C_vel/10;C_dist/100;C_minheight;C_turnback/10]
    # println("3")
    # println(sqrt(integ_1.u[3]^2+integ_1.u[4]^2))

    global printIter
    printIter+=1
    # println(printIter)
    if (printIter%(length(x)*5)==0.0) || (printIter == 1) #print out every 5 major iterations
        println("
        x: $x

        J: $J  dt: $dt

        C_maxheight: $C_maxheight

        C_vel: $C_vel

        C_minheight: $C_minheight

        C_turnback: $C_turnback

        C_dist: $C_dist

        C_power: $C_power

        propulsion_energy: $(sum(power*dt))

        potential_energy : $(-P[:gravity]*P[:mass]*d_height)

        ")

        figure("path")
        clf()
        plot(savestates[:,1],-savestates[:,2],".-",label = "$printIter")
        xlabel("Distance (m)")
        ylabel("Height (m)")
        legend(loc = "best")
        pause(0.001)
    end

    if optimize
        return J,C,false
    else
        return errorsum,savestates,C_maxheight,C_vel,C_minheight,CL,CD,alpha_save
    end
end


# x0  = [P[:pn0],P[:pe0],P[:pd0],P[:u0],P[:v0],P[:w0],P[:phi0],P[:theta0],P[:psi0],P[:p0],P[:q0],P[:r0]]
# x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0]]
x0_states  = [P[:pn0],P[:pd0],P[:u0],P[:w0],P[:theta0],P[:q0]]

# Set up ContinuousDynamicalSystem
ds = ContinuousDynamicalSystem(testUAV!, x0_states, P)
#Reinitizilze the integrator

# delta_e = [-85.0, 1.16896, 52.1414, -32.6719, -2.07146, -15.4012, -9.35701, -11.8722, -10.8659, -11.2834, -11.2427, -10.1897, -12.6121, -9.07863, -13.8728, -6.67249, -19.487, 9.64414, -85.0, -85.0]
delta_e = ones(40)*-10*pi/180
TW_var = ones(delta_e)*0.8

h_set = 20.0
h_M = 10.0/13.5
momentum = P[:Va0]*P[:mass]
dt = .10
TW = 0.8


lb_de = ones(delta_e)*-85*pi/180
ub_de = ones(delta_e)*85*pi/180

lb_TW = ones(delta_e)*0.0
ub_TW = ones(delta_e)*1.0


# x0 = [1.1, -0.357509, 0.0247765, -0.0444137, 0.0545875, -0.00775476, 0.129438, -0.0316016, 0.0352566, -0.133037, 0.0291004, -0.112904, -0.353401, -0.245529, 0.00146116, 0.00579712, -0.297543, -0.0167796, -0.262352, 0.106507, -0.112886, 0.506774, 0.391332, 0.356882, 0.280428, 0.208753, 0.169694, 0.145843, 0.0826858, 0.0763246, 0.178571, 0.258273, 0.331779, 0.36597, 0.406327, 0.424674, 0.474412, 0.380176, 0.381731, 0.316303, 0.295168]
x0 = [.459487, -0.123406, -0.115946, -0.168572, -0.130133, -0.205725, -0.0787378, -0.212405, -0.0901605, -0.168617, -0.13109, -0.143469, -0.143642, -0.143034, -0.149085, -0.146726, -0.150737, -0.150169, -0.150887, -0.149252, -0.148149, -0.14311, -0.136763, -0.127008, -0.118405, -0.106209, -0.0969768, -0.087155, -0.0830964, -0.0802628, -0.0880927, -0.0979006, -0.132952, -0.154202, -0.163709, -0.168885, -0.162402, -0.14449, -0.0908417, -0.173984, -1.48353, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# x0 = ones(length(delta_e)*2+1)
# x0[1] = 0.5
x0[62:end] = x0[42:61]
x0[22:42] = x0[2:22]


scale = ones(x0)
scale[1] = 10.0
x0 = x0.*scale
# x0 = [.1;delta_e]
# lb = [1E-3;lb_de]
# ub = [5.0;ub_de]
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
options["Minor feasibility tolerance"] = 1e-5
options["Minor print level"] = 5
options["Print frequency"] = 100
options["Scale option"] = 1
options["Scale tolerance"] = .95



PW_array = [17.0,13.6,6.8,3.4]#linspace(.459,.46,10)
TW_array = reverse([.2,.4,.8,1.0])
# distance_factor = [.7,.9,1.0]
TWmax = 1.0
distance_factor = [5.0]
distance = []
# args = []
args2 = []
optimize = []
savestates = []
time_sim = []
t_save = zeros(PW_array)
errorsum = zeros(PW_array)
mean_thrust2weight = zeros(PW_array)
propulsion_total_energy = zeros(PW_array)
final_dist = zeros(PW_array)
energy_an = zeros(PW_array)

save_savestates = zeros(length(delta_e),7,length(PW_array))
alpha_num = zeros(length(delta_e),length(PW_array))

function zero_height(x)
    return  zero_height(x,args2...)
end

function zero_height(x,alt_spl,h_set)
    h = alt_spl(x)
    return h-h_set
end
testrun = false
printIter = 0 #this is a global that triggers a print every few major iterations

# for i = 1:length(PW_array)
    i = 1
    distance = 300.0#h_set*distance_factor[1]
    # P[:thrust_weight] = 0.5
    P[:power_weight] = PW_array[i]
    # TW = distance_factor[i]
    # P[:mass] = Mass_array[i]
    # Mass = Mass_array[i]
    # P[:Va0] = momentum/Mass
    # println(P[:Va0])
    # h_set = h_M*Mass
    printIter = 0 #this is a global that triggers a print every few major iterations
    if !testrun
        optimize = true
        lb = [1E-3;lb_de;lb_TW].*scale
        ub = [5.0;ub_de;ub_TW].*scale

        xopt, fopt, optinfo = snopt(objcon, x0, lb, ub, options;printfile = "$(fileLoc)/SNOPT_out/snopt-print$(i).out", sumfile = "$(fileLoc)/SNOPT_out/snopt-summary.out") #!!!!!! ARGS are in objcon wrapper function !!!!!!!!#

        println("xopt: $(xopt)")
        println("####################################################")
        println(optinfo)
        println("####################################################")
    else
        xopt = x0
    end
    optimize = false
    dt = x0[1]/scale[1]
    time_sim = collect(dt:dt:dt*length(delta_e))

    errorsum[i],savestates,C_maxheight,C_vel,C_minheight,CL,CD,alpha_save = objcon(xopt)
    mean_thrust2weight[i] = mean(savestates[:,end]./(ones(time_sim)*P[:mass]*P[:gravity]))
    Va = sqrt.(savestates[:,3].^2+savestates[:,4].^2)

    save_savestates[:,:,i] = savestates
    alpha_num[:,i] = alpha_save
    # Va_min = maximum(Va)


    #-------- PLOTS -------#
    figname = "d_CL"
    figure(figname)
    plot(savestates[:,1],CL,".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Lift Coefficient")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "d_CD"
    figure(figname)
    plot(savestates[:,1],CD,".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Drag Coefficient")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "d_q"
    figure(figname)
    plot(savestates[:,1],savestates[:,6],".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Pitch Rate")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "d_theta"
    figure(figname)
    plot(savestates[:,1],savestates[:,5]*180/pi,".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Pitch (deg)")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "d_AOA"
    figure(figname)
    plot(savestates[:,1],alpha_save,".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Angle of Attack (deg)")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "d_Va"
    figure(figname)
    plot(savestates[:,1],Va,".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Airspeed (m/s)")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)
    # ylim(0,32)

    figname = "pn_pd"
    figure(figname)
    plot(savestates[:,1],-savestates[:,2],".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    # plot(x_an,y_an,"-",label = "analytical $(round.(PW_array[i],4))")
    xlabel("Distance (m)")
    ylabel("Height (m)")
    # axis("equal")
    # xlim(minimum(savestates[:,1]),maximum(savestates[:,1]))
    # ylim(minimum(-savestates[:,2]),maximum(-savestates[:,2]))
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "pn_pd_norm"
    figure(figname)
    plot(savestates[:,1]/maximum(savestates[:,1]),-savestates[:,2]/maximum(-savestates[:,2]),".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("pn/max pn")
    ylabel("-pd/max -pd")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)
    # xlim(0,70)
    # ylim(0,20)

    figname = "t_T"
    figure(figname)
    plot(time_sim,savestates[:,end]./(P[:mass]*P[:gravity]),".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Time (s)")
    ylabel("Thrust/Weight")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    power = Va.*savestates[:,end]

    figname = "d_power"
    figure(figname)
    plot(savestates[:,1],power./(P[:mass]*P[:gravity]),".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Power/Weight")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    figname = "d_delta_e"
    figure(figname)
    plot(savestates[:,1],xopt[2:21],".-",label = "$(PW_array[i]) | $((TW_array[i]))")
    xlabel("Distance (m)")
    ylabel("Elevator Deflection (deg)")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)

    alt_spl = Dierckx.Spline1D(savestates[:,1],-savestates[:,2])
    va_spl = Dierckx.Spline1D(savestates[:,1],Va)
    T_spl = Dierckx.Spline1D(savestates[:,1],savestates[:,end])
    t_spl = Dierckx.Spline1D(savestates[:,1],time_sim)
    P_spl = Dierckx.Spline1D(time_sim,power)

    t_save[i] = time_sim[end]
    #
    args2 = (alt_spl,h_set)
    dist_at_height = 0.0
    try
        dist_at_height = fzero(zero_height,0,maximum(savestates[:,1]))
    catch
        dist_at_height = maximum(savestates[:,1])
        if -savestates[end,2]<h_set*.999
            warn("Max height was $(-savestates[end,2]), should be $h_set")
        end
    end
    time_at_height = t_spl(dist_at_height)
    propulsion_total_energy[i] = integrate(P_spl,0,time_at_height)
    final_dist[i] = dist_at_height

# end

# figname = "errorsum"
# figure(figname)
# plot(PW_array,errorsum,".-",label = "$(round.(PW_array,4))")
# xlabel("Power/Weight Constraint")
# ylabel("Error Sum")
# savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)
#
# figname = "propulsion_total_energy"
# figure(figname)
# plot(PW_array,propulsion_total_energy/(P[:mass]*P[:gravity]),".-",label = "Numerical")
# # plot(PW_array,energy_an/(P[:mass]*P[:gravity]),".-",label = "Analytical")
# # legend(loc = "best")
# xlabel("Power/Weight Constraint")
# ylabel("Energy/Weight (J/N)")
# savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)
#
# figname = "average_slope"
# figure(figname)
# plot(PW_array,h_set./final_dist,".-",label = "$(round.(PW_array,4))")
# xlabel("Power/Weight Constraint")
# ylabel("Average Slope (H/D)")
# savefig("$(fileLoc)/../figures/powered_glide/$figname.pdf",transparent = true)
