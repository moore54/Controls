using PyPlot
fileLoc,_ = splitdir(@__FILE__)
PyPlot.close("all")
rc("figure", figsize=(6.5, 4.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false,)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.69)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])

colors=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

# include("$fileLoc/STOL_dynamics_delta_e.jl")

m = 5.0 #kg
g = 9.81
L = m*g*2
D = L/30

alpha = 15.0*pi/180
Va = 17.0
H_c = 20.0

TW_array = reverse([.2,.4,.8,1.0])
for i = 1:length(TW_array)
    T = m*g*TW_array[i]


    gama = asin((T*cos(alpha)-D)/(m*g))
    gamad = gama*180/pi

    G = linspace(0,gama,100)-pi/2
    Fy = L*cos(gama) - D*sin(gama) + T*sin(alpha+gama) - m*g
    Fx = -D*cos(gama) - L*sin(gama) + T*cos(alpha+gama)
    Fl = L - m*g*cos(gama) + T*sin(alpha)
    r = m*Va^2/Fl

    arclength = r*gama
    t_tr = arclength/Va



    xtr = r*cos.(G)
    ytr = r*sin.(G)+r

    H_left = H_c-ytr[end]
    x_left = H_left/tan(gama)
    xcl = [xtr[end],xtr[end]+x_left]
    ycl = [ytr[end],H_c]

    xrol = [0,x_left/3.0]
    yrol = [0,0]

    rc("figure", figsize=(7.5, 1.6))
    rc("figure.subplot", left=0.1, bottom=0.27, top=0.97, right=.7)
    figure("position")
    plot(xrol,yrol,"-",color = colors[i])#,label = "roll")
    plot(xtr+xrol[end],ytr,"--",color = colors[i])#,label = "transition")
    plot(xcl+xrol[end],ycl,"-",color = colors[i],label = "$(round.(TW_array[i]*Va,1)) | $((TW_array[i]))")#,label = "climb")
    legend(loc="center left", title = "Power/Weight | Thrust/Weight",bbox_to_anchor=(1, 0.5))
    xlabel("Distance (m)")
    ylabel("Height (m)")
    # axis("equal")
end
figname = "position"
figure("position")
savefig("./figures/free_analytical/$figname.pdf",transparent = true)
