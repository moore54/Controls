fileLoc,_ = splitdir(@__FILE__)
print_iters = 0
using Snopt
using Dierckx
using Roots
using PyPlot
using Gradients
close("all")
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
    delta_e = x[2:21]
    TW = x[22:41]

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

    power = Va.*savestates[:,end]
    J = sum(power*dt)/1E4 #savestates[end,1]/10#errorsum/100
    C_dist = savestates[end,1]-distance
    C_minheight = 0.0-(-savestates[:,2])
    C_maxheight = h_set - (maximum(-savestates[end,2]))
    # C_endslope = -savestates[end,2]-(-savestates[end-1,2])
    C_power = power/(P[:mass]*P[:gravity])-P[:power_weight]
    # Cvel_max = Va-integ_1.p[:Va0]*10.0001
    C = [C_maxheight;C_vel/10;C_minheight;C_turnback/10;C_dist;C_power]
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


        ")

        # figure("path")
        # plot(savestates[:,1],-savestates[:,2],".-")
        # xlabel("pn (m)")
        # ylabel("-pd (m)")
    end

    if optimize
        return J,C,false
    else
        return errorsum,savestates,C_maxheight,C_vel,C_minheight,CL,CD,alpha_save
    end
end


# x0  = [P[:pn0],P[:pe0],P[:pd0],P[:u0],P[:v0],P[:w0],P[:phi0],P[:theta0],P[:psi0],P[:p0],P[:q0],P[:r0]]
# x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0]]
x0  = [P[:pn0],P[:pd0],P[:u0],P[:w0],P[:theta0],P[:q0]]

# Set up ContinuousDynamicalSystem
ds = ContinuousDynamicalSystem(testUAV!, x0, P)
#Reinitizilze the integrator

# delta_e = [-85.0, 1.16896, 52.1414, -32.6719, -2.07146, -15.4012, -9.35701, -11.8722, -10.8659, -11.2834, -11.2427, -10.1897, -12.6121, -9.07863, -13.8728, -6.67249, -19.487, 9.64414, -85.0, -85.0]
delta_e = ones(20)*-10*pi/180
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


x0 = [.157933, -1.09447, 0.865694, -0.470354, -0.0576175, -0.221146, -0.137922, -0.187087, -0.1693, -0.108763, -0.216798, 0.0355408, -0.126528, -0.158557, -0.0667918, 0.0105445, -0.158397, 0.135525, -0.121021, -0.155036, 0.924553, 0.572882, 0.552737, 0.537172, 0.422904, 0.349736, 0.296045, 0.281233, 0.257718, 0.249144, 0.314345, 0.358956, 0.463606, 0.5434, 0.558226, 0.56457, 0.589074, 0.563715, 0.527053, 0.476165, 0.416111]
# x0 = x0[1:21]
# x0 = [.1;delta_e;TW_var]
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
options["Major optimality tolerance"] = 1e-3 #Should be scaled so it is optimal with a solid 2 digits
options["Minor optimality  tolerance"] = 1e-6
options["Major feasibility tolerance"] = 1e-5
options["Minor feasibility tolerance"] = 1e-5
options["Minor print level"] = 5
options["Print frequency"] = 100
options["Scale option"] = 1
options["Scale tolerance"] = .95



PW_array = [20.0,10.0,5.0,2.5]#linspace(.459,.46,10)
# distance_factor = [.7,.9,1.0]
TWmax = 1.0
distance_factor = [5.0]
distance = []
# args = []
args2 = []
optimize = []
savestates = []
time_sim = []
errorsum = zeros(PW_array)
mean_thrust2weight = zeros(PW_array)
total_energy = zeros(PW_array)
final_dist = zeros(PW_array)
energy_an = zeros(PW_array)

function zero_height(x)
    return  zero_height(x,args2...)
end

function zero_height(x,alt_spl,h_set)
    h = alt_spl(x)
    return h-h_set
end
testrun = true
printIter = 0 #this is a global that triggers a print every few major iterations

for i = 1:length(PW_array)
    # i = 1
    distance = h_set*distance_factor[1]
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

        xopt, fopt, optinfo = snopt(objcon, x0, lb, ub, options;printfile = "$(fileLoc)/snopt-print$(distance).out", sumfile = "$(fileLoc)/snopt-summary.out") #!!!!!! ARGS are in objcon wrapper function !!!!!!!!#

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
    Va_min = maximum(Va)

    #-------- ANALYTICAL -------#
    T = mean(xopt[21:end])*P[:mass]*P[:gravity]
    L = maximum(CL)*0.5*P[:rho]*Va_min^2*P[:S_wing]
    D = minimum(CD)*0.5*P[:rho]*Va_min^2*P[:S_wing]
    gama = asin((T-D)/(P[:mass]*P[:gravity]))
    gamad = gama*180/pi
    alpha = maximum(alpha_save*pi/180)

    G = linspace(0,gama,100)-pi/2

    r = P[:mass]*Va_min^2/(L+T*sin(alpha)-P[:mass]*P[:gravity])

    arclength = r*gama
    t_tr = arclength/Va_min

    xtr = r*cos.(G)
    ytr = r*sin.(G)+r

    H_left = h_set-ytr[end]
    x_left = H_left/tan(gama)
    xcl = [xtr[end],xtr[end]+x_left]
    ycl = [ytr[end],ytr[end]+H_left]
    t_cl = sqrt(x_left^2+H_left^2)/Va_min
    t_an = t_tr+t_cl
    p_an = T*Va_min
    energy_an[i] = p_an*t_an
    x_an = [xtr;xcl]
    y_an = [ytr;ycl]

    #-------- PLOTS -------#
    figname = "d_CL"
    figure(figname)
    plot(savestates[:,1],CL,".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("CL")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "d_CD"
    figure(figname)
    plot(savestates[:,1],CD,".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("CD")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "d_q"
    figure(figname)
    plot(savestates[:,1],savestates[:,6],".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("q")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "d_theta"
    figure(figname)
    plot(savestates[:,1],savestates[:,5]*180/pi,".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("theta (deg)")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "d_AOA"
    figure(figname)
    plot(savestates[:,1],alpha_save,".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("AOA (deg)")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "d_Va"
    figure(figname)
    plot(savestates[:,1],Va,".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("Va (m/s)")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)
    # ylim(0,32)

    figname = "pn_pd"
    figure(figname)
    plot(savestates[:,1],-savestates[:,2],".-",label = "dynamic $(round.(PW_array[i],4))")
    plot(x_an,y_an,"-",label = "analytical $(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("-pd (m)")
    axis("equal")
    # xlim(minimum(savestates[:,1]),maximum(savestates[:,1]))
    # ylim(minimum(-savestates[:,2]),maximum(-savestates[:,2]))
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "pn_pd_norm"
    figure(figname)
    plot(savestates[:,1]/maximum(savestates[:,1]),-savestates[:,2]/maximum(-savestates[:,2]),".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn/max pn")
    ylabel("-pd/max -pd")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)
    # xlim(0,70)
    # ylim(0,20)

    figname = "t_T"
    figure(figname)
    plot(time_sim,savestates[:,end]./(P[:mass]*P[:gravity]),".-",label = "$(round.(PW_array[i],4))")
    xlabel("time (s)")
    ylabel("Inverse Slope")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    power = Va.*savestates[:,end]

    figname = "d_power"
    figure(figname)
    plot(savestates[:,1],power./(P[:mass]*P[:gravity]),".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("Power/weight")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    figname = "d_delta_e"
    figure(figname)
    plot(savestates[:,1],xopt[2:21],".-",label = "$(round.(PW_array[i],4))")
    xlabel("pn (m)")
    ylabel("Elevator Deflection (deg)")
    legend(loc="center left", title = "Power/Weight",bbox_to_anchor=(1, 0.5))
    savefig("./figures/free_analytical/$figname.pdf",transparent = true)

    alt_spl = Dierckx.Spline1D(savestates[:,1],-savestates[:,2])
    va_spl = Dierckx.Spline1D(savestates[:,1],Va)
    T_spl = Dierckx.Spline1D(savestates[:,1],savestates[:,end])
    t_spl = Dierckx.Spline1D(savestates[:,1],time_sim)
    P_spl = Dierckx.Spline1D(time_sim,power)
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
    total_energy[i] = integrate(P_spl,0,time_at_height)
    final_dist[i] = dist_at_height

end

figname = "errorsum"
figure(figname)
plot(PW_array,errorsum,".-",label = "$(round.(PW_array,4))")
xlabel("Thrust/Weight")
ylabel("errorsum")
savefig("./figures/free_analytical/$figname.pdf",transparent = true)

figname = "total_energy"
figure(figname)
plot(PW_array,total_energy/(P[:mass]*P[:gravity]),".-",label = "Numerical")
plot(PW_array,energy_an/(P[:mass]*P[:gravity]),".-",label = "Analytical")
legend(loc = "best")
xlabel("Thrust/Weight")
ylabel("Energy/Weight (J/N)")
savefig("./figures/free_analytical/$figname.pdf",transparent = true)

figname = "average_slope"
figure(figname)
plot(PW_array,h_set./final_dist,".-",label = "$(round.(PW_array,4))")
xlabel("Thrust/Weight")
ylabel("Average Slope (H/D)")
savefig("./figures/free_analytical/$figname.pdf",transparent = true)
