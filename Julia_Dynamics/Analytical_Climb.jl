using PyPlot
PyPlot.close("all")

rc("figure", figsize=(6.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.69)
# rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
# TW = linspace(.5,.999,100)
# height = 20
# dist = height * sqrt.(1.0./TW.^2.0-1)
#
# PyPlot.figure("1")
# PyPlot.plot(TW,dist)
# PyPlot.xlabel("thrust/weight")
# PyPlot.ylabel("distance traveled (m) for $(height)m elevation")

height = 20.0 #m
mass = 13.5 #kg
gravity = 9.81
weight = mass*gravity

CLa = 6.87
CL0 = 0.0
rho = 1.18
Sref = .55
AR = 15
e = 1.0

N_PW = 30
N_Va = 30
Dp = 1.0
PW_array = linspace(4.0,20.0,N_PW)
Va_array = linspace(5,30,N_Va)
dist_save = zeros(N_PW*N_Va)
Va_save = zeros(dist_save)
TW_save = zeros(dist_save)
PW_save = zeros(dist_save)
alpha_save = zeros(dist_save)


drag = []
CL  =[]
TD = []
k = 1
dist = []
dist_star = 100

for j = 1:N_Va
   Va = Va_array[j]

   for i = 1:N_PW
      PW = PW_array[i]
      thrust = PW*weight/Va #N
      if thrust>weight
         thrust = 0.99*weight
      end
      err = 1E20

      while err>1E-6
         drag = (mass*gravity)^2 / ((1+(height/dist_star)^2)*0.5*rho*Va^2*Sref*pi*e*AR) + Dp

         dist = height*sqrt((mass*gravity/(thrust-drag))^2.0-1.0)
         err=abs(dist_star-dist)
         dist_star = copy(dist)



      end

      flight_angle = asin((thrust-drag)/(mass*gravity))
      lift = sqrt((mass*gravity)^2-(thrust-drag)^2)
      CL = lift/(0.5*rho*Va^2*Sref)
      alpha = (CL-CL0)/CLa
      # println(alpha_array[i]*180/pi)


      if alpha*180/pi>15 || thrust/weight > 0.989
         dist_save[k] = -500000.0
         Va_save[k] = Va
         TW_save[k] = thrust/weight#PW
         PW_save[k] = PW
         alpha_save[k] = -10.0
      else
         dist_save[k] = dist
         Va_save[k] = Va
         TW_save[k] = thrust/weight#PW
         PW_save[k] = PW
         alpha_save[k] = alpha
      end


      # PyPlot.semilogy(alpha_array[i]*180/pi,dist_save[i],"o",label = "thrust/weight: $(round(PW,3))")
      k+=1
   end

# println(alpha_array)


end
PyPlot.figure("Distance_Thrust_Weight")
PyPlot.scatter3D(Va_save,TW_save,(dist_save))
PyPlot.xlabel("Airspeed (m/s)")
PyPlot.ylabel("Thrust / Weight")
PyPlot.zlabel("Climb Distance (m)")
PyPlot.zlim([0,maximum(dist_save)])


PyPlot.figure("AOA_Thrust_Weight")
PyPlot.scatter3D(Va_save,TW_save,alpha_save*180/pi)
PyPlot.xlabel("Airspeed (m/s)")
PyPlot.ylabel("Thrust / Weight")
PyPlot.zlabel("AOA (deg)")
PyPlot.zlim([0,maximum(alpha_save*180/pi)])

PyPlot.figure("Distance_Power_Weight")
PyPlot.scatter3D(Va_save,PW_save,(dist_save))
PyPlot.xlabel("Airspeed (m/s)")
PyPlot.ylabel("Power / Weight")
PyPlot.zlabel("Climb Distance (m)")
PyPlot.zlim([0,maximum(dist_save)])


PyPlot.figure("AOA_Power_Weight")
PyPlot.scatter3D(Va_save,PW_save,alpha_save*180/pi)
PyPlot.xlabel("Airspeed (m/s)")
PyPlot.ylabel("Power / Weight")
PyPlot.zlabel("AOA (deg)")
PyPlot.zlim([0,maximum(alpha_save*180/pi)])
