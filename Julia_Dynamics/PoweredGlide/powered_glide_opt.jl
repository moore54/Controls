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
    delta_e = x[2:81]
    TW = x[82:161]

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
    C_minheight = 996.0-(-savestates[:,2]) # height cannot dip below 0
    C_maxheight = (-savestates[:,2]) - 1000.0 # height cannot dip below 0
    C_endheight = 1000.0 - (maximum(-savestates[end,2]))
    C_endslope = -savestates[end,2]-(-savestates[end-1,2])
    C_power = power/(P[:mass]*P[:gravity])-P[:power_weight]
    Cvel_max = Va-integ_1.p[:Va0]*1.01
    # C = [C_vel/10;C_minheight;C_turnback/10;C_power;Cvel_max]
    C = [C_vel/10;C_dist/100;C_minheight;C_maxheight;C_turnback/10]
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
delta_e = ones(80)*-10*pi/180
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


# No lower bound optimal
x0_no_LB = [0.466575, -0.116498, -0.126046, -0.150804, -0.16734, -0.148954, -0.148027, -0.142953, -0.141099, -0.14214, -0.139366, -0.141895, -0.141678, -0.144924, -0.145326, -0.146451, -0.146114, -0.14771, -0.149292, -0.145477, -0.150261, -0.144572, -0.149451, -0.145335, -0.147801, -0.143293, -0.146998, -0.144187, -0.144631, -0.145323, -0.144011, -0.148223, -0.147968, -0.148988, -0.151256, -0.159734, -0.155608, -0.15543, -0.175694, -0.0972616, -0.0963656, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# lower alt bound
x0_LB = [0.460027, -0.128612, -0.114631, -0.157402, -0.159845, -0.14763, -0.14319, -0.144587, -0.13773, -0.136452, -0.135793, -0.131173, -0.132128, -0.125613, -0.128362, -0.122043, -0.125165, -0.122173, -0.124378, -0.126942, -0.128556, -0.137935, -0.135962, -0.143039, -0.142999, -0.145928, -0.144939, -0.145919, -0.145027, -0.145399, -0.145198, -0.146337, -0.147997, -0.150511, -0.153326, -0.155256, -0.157887, -0.156925, -0.160871, -0.148076, 0.136099, 0.13054, 0.0146938, 0.0, 0.00494878, 0.0363091, 0.0692158, 0.0930316, 0.105056, 0.117837, 0.130434, 0.140369, 0.148327, 0.150682, 0.150128, 0.141926, 0.131029, 0.113666, 0.0942773, 0.0718675, 0.0483467, 0.0259255, 0.00358225, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#lower alt bound 80 control inputs
x0_LB = [0.459479, -0.120404, -0.123309, -0.148089, -0.161699, -0.14119, -0.140543, -0.142273, -0.133867, -0.134716, -0.132838, -0.128931, -0.127049, -0.122731, -0.119943, -0.115811, -0.113371, -0.109881, -0.108048, -0.106856, -0.10565, -0.108182, -0.109633, -0.115414, -0.122439, -0.134465, -0.134785, -0.142876, -0.146605, -0.147371, -0.147843, -0.148985, -0.146496, -0.148712, -0.143347, -0.150549, -0.138101, -0.153896, -0.134097, -0.155608, -0.133157, -0.154598, -0.136057, -0.15018, -0.141222, -0.145691, -0.144623, -0.14391, -0.145207, -0.144232, -0.144662, -0.144827, -0.144377, -0.144872, -0.144489, -0.144892, -0.144513, -0.144783, -0.144786, -0.144761, -0.144559, -0.144502, -0.144409, -0.144054, -0.14361, -0.143347, -0.142995, -0.142501, -0.142236, -0.142439, -0.143171, -0.144271, -0.14645, -0.14921, -0.15254, -0.15441, -0.157479, -0.156136, -0.161032, -0.145347, 0.13426, 0.139224, 0.027914, 0.0, 0.0223178, 0.0591317, 0.0979931, 0.126368, 0.145521, 0.170509, 0.197784, 0.224291, 0.249733, 0.270634, 0.287098, 0.296847, 0.300167, 0.295199, 0.282699, 0.261904, 0.233569, 0.198994, 0.158951, 0.117268, 0.074614, 0.0332183, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#Lower and upper bound test
[0.451721, -0.131648, -0.111727, -0.156306, -0.163182, -0.149912, -0.144753, -0.136418, -0.137234, -0.13515, -0.130278, -0.137436, -0.128262, -0.136958, -0.12834, -0.133284, -0.128909, -0.12969, -0.12835, -0.127717, -0.127249, -0.127148, -0.127052, -0.127636, -0.128435, -0.129718, -0.131356, -0.131862, -0.131822, -0.133869, -0.133202, -0.135803, -0.134872, -0.138276, -0.136568, -0.140584, -0.137487, -0.141869, -0.137727, -0.141146, -0.137318, -0.13863, -0.135501, -0.134961, -0.133402, -0.134187, -0.132219, -0.132625, -0.131521, -0.131062, -0.130116, -0.129383, -0.128586, -0.127964, -0.127482, -0.127222, -0.127245, -0.127552, -0.128308, -0.129285, -0.131177, -0.132272, -0.132002, -0.133054, -0.133617, -0.134349, -0.135398, -0.13669, -0.138177, -0.139775, -0.141535, -0.143405, -0.145773, -0.148379, -0.152251, -0.154274, -0.156782, -0.155776, -0.161005, -0.144199, 0.136713, 0.133334, 0.01252, 0.0, 0.0, 0.0183671, 0.0467506, 0.068962, 0.0809905, 0.081327, 0.0803564, 0.0836426, 0.0854448, 0.0931007, 0.0960254, 0.102654, 0.105408, 0.108271, 0.108874, 0.107645, 0.104679, 0.0998864, 0.0934635, 0.0859446, 0.0775287, 0.0686438, 0.0597044, 0.0512252, 0.044182, 0.0381267, 0.031336, 0.0254486, 0.0198287, 0.0166705, 0.0143941, 0.0154481, 0.0170621, 0.0224357, 0.0272544, 0.0355113, 0.0425057, 0.0514397, 0.0591237, 0.0671818, 0.0738792, 0.0795071, 0.0838444, 0.0892768, 0.0945164, 0.0994971, 0.103459, 0.106214, 0.107417, 0.106998, 0.104801, 0.100899, 0.0955537, 0.0889659, 0.0814688, 0.0733855, 0.0650635, 0.0567154, 0.0492763, 0.0431623, 0.0365602, 0.0292178, 0.0218779, 0.0153268, 0.0101033, 0.00659079, 0.00497307, 0.00498537, 0.00608715, 0.00718903, 0.00705149, 0.00356669, 0.0, 0.0, 0.0, 0.0, 0.0]
x0 = x0_LB




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
    distance = 600.0#h_set*distance_factor[1]
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
    plot(savestates[:,1],xopt[2:41],".-",label = "$(PW_array[i]) | $((TW_array[i]))")
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
