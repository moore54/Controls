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
e_wing = 1.0

N_PW = 60
# N_CL = 10
Dp = 1.0
CDp = 0.01
TW_array = linspace(.2,.9999,N_PW)
CL_array = [.4,.6,.8,1.0,1.2,1.5,2.0,5.0]
N_CL = length(CL_array)
dist_save = zeros(N_PW*N_CL)
Va_save = zeros(dist_save)
TW_save = zeros(dist_save)
PW_save = zeros(dist_save)
CL_save = zeros(dist_save)



drag = []
CL  =[]
TD = []
thrust = []
Va = []
k = 1
dist = []
dist_star = 100

for j = 1:N_CL
   CL = CL_array[j]
   dist_plot = [] # zeros(N_PW)
   TW_plot = [] # zeros(PW_array)
   D_perc_plot = [] # zeros(PW_array)
   flight_angle_plot = [] # zeros(PW_array)
   Va_plot = [] # zeros(PW_array)
   Power_plot = [] # zeros(PW_array)
   alpha_plot = [] # zeros(PW_array)
   PW_plot = []
   energy_plot = []
   time_plot = []
   for i = 1:N_PW
      # PW = PW_array[i]
      TW = TW_array[i]

      err = 1E20
      while err>1E-6

         Va = sqrt(2*mass*gravity/(CL*rho*Sref*sqrt(1+(height/dist_star)^2)))

         thrust = TW*weight#/Va #N
         if thrust/weight> 1-1E-5
            break
            thrust = (1-1E-5)*weight
            break
         end


         CD = (CL)^2 / (pi*e_wing*AR) + CDp
         drag = CD*0.5*rho*Va^2*Sref

         dist = height*sqrt((mass*gravity/(thrust-drag))^2.0-1.0)
         err=abs(dist_star-dist)
         dist_star = copy(dist)
      end




      if thrust/weight > 1-1E-5
         break
         # dist_save[k] = 0.0
         # dist_plot[i] = 0.0
         # Va_plot[i] = 0.0# Va
         # Power_plot[i] = 0.0# thrust*Va
         # TW_save[k] = 0.0# thrust/weight#PW
         # TW_plot[i] = 0.0# thrust/weight#PW
         # PW_save[k] = 0.0# PW
         # alpha_plot[i] = 0.0
         # CL_save[k] = 0.0
         # D_perc_plot[i] = 0
         # flight_angle_plot[i] = 0.0#flight_angle
         break
      else
         flight_angle = asin((thrust-drag)/(mass*gravity))
         lift = sqrt((mass*gravity)^2-(thrust-drag)^2)
         CL = lift/(0.5*rho*Va^2*Sref)
         alpha = (CL-CL0)/CLa
         flight_dist = sqrt(height^2+dist^2)
         time = flight_dist/Va

         dist_save[k] = dist
         push!(dist_plot, dist)#[i] = dist
         push!(Va_plot, Va)#[i] = Va
         push!(Power_plot, thrust*Va)#[i] = thrust*Va
         TW_save[k] = thrust/weight
         push!(TW_plot, thrust/weight)#[i] = thrust/weight
         PW_save[k] = thrust*Va/weight
         push!(alpha_plot, alpha)#[i] = alpha
         CL_save[k] = CL
         push!(D_perc_plot, 1-(thrust-drag)/thrust)#[i] = 1-(thrust-drag)/thrust
         push!(flight_angle_plot, flight_angle)#[i] = flight_angle
         push!(PW_plot,thrust*Va/weight)
         push!(energy_plot,time*thrust*Va)
         push!(time_plot,time)
      end

      k+=1
   end
# PyPlot.figure("PW")
# PyPlot.plot(PW_plot,dist_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Power / Weight (w/N)")
# PyPlot.ylabel("Climb Distance (m)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
#

PyPlot.figure("TW")
PyPlot.plot(TW_plot,height./dist_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Thrust / Weight (N/N)")
PyPlot.ylabel("Slope")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

# PyPlot.figure("Net_TW")
# PyPlot.plot((TW_plot.*weight-D_perc_plot.*TW_plot.*weight)/weight,height./dist_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Net Thrust / Weight (N/N)")
# PyPlot.ylabel("Slope")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

PyPlot.figure("D_perc")
PyPlot.plot(TW_plot,D_perc_plot*100,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Thrust / Weight")
PyPlot.ylabel("Drag Percent of Thrust (%)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

# PyPlot.figure("Net_D_perc")
# PyPlot.plot((TW_plot.*weight-D_perc_plot.*TW_plot.*weight)/weight,D_perc_plot*100,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Net Thrust / Weight")
# PyPlot.ylabel("Drag Percent of Thrust (%)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

PyPlot.figure("flight_angle")
PyPlot.plot(TW_plot,flight_angle_plot*180/pi,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Thrust / Weight")
PyPlot.ylabel("Flight Path Angle (deg)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))


PyPlot.figure("flight_angle_power")
PyPlot.plot(PW_plot,flight_angle_plot*180/pi,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Power / Weight")
PyPlot.ylabel("Flight Path Angle (deg)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))


# PyPlot.figure("net_flight_angle")
# PyPlot.plot((TW_plot.*weight-D_perc_plot.*TW_plot.*weight)/weight,flight_angle_plot*180/pi,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Net Thrust / Weight")
# PyPlot.ylabel("Flight Path Angle (deg)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))


# PyPlot.figure("Va")
# PyPlot.plot(TW_plot,Va_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Thrust / Weight")
# PyPlot.ylabel("Flight Speed Va (m/s)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

PyPlot.figure("Va_dist")
PyPlot.plot(Va_plot,height./dist_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Flight Speed Va (m/s)")
PyPlot.ylabel("Slope")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))



PyPlot.figure("Va_power")
PyPlot.plot(PW_plot,Va_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Power / Weight")
PyPlot.ylabel("Flight Speed Va (m/s)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))


PyPlot.figure("Va_TW")
PyPlot.plot(TW_plot,Va_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Thrust / Weight")
PyPlot.ylabel("Flight Speed Va (m/s)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

PyPlot.figure("Power_thrust")
PyPlot.plot(TW_plot,PW_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Thrust / Weight")
PyPlot.ylabel("Power / Weight")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("Power_power")
# PyPlot.plot(PW_plot,Power_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Power / Weight")
# PyPlot.ylabel("Power (Watts)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
PyPlot.figure("Power_dist")
PyPlot.plot(PW_plot,height./dist_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Power / Weight")
PyPlot.ylabel("Slope")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))


# PyPlot.figure("AOA")
# PyPlot.plot(TW_plot,alpha_plot*180/pi,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Thrust / Weight")
# PyPlot.ylabel("AOA (deg)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

# PyPlot.figure("AOA_Va")
# PyPlot.plot(Va_plot,alpha_plot*180/pi,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Va (m/s)")
# PyPlot.ylabel("AOA (deg)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))



PyPlot.figure("Energy_Climb")
PyPlot.plot(height./dist_plot,energy_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Slope")
PyPlot.ylabel("Energy (J)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

PyPlot.figure("Energy_TW")
PyPlot.plot(TW_plot,energy_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Thrust / Weight")
PyPlot.ylabel("Energy (J)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("Energy_PW")
# PyPlot.plot(PW_plot,energy_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Power / Weight")
# PyPlot.ylabel("Energy (J)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("Time_TW")
# PyPlot.plot(TW_plot,time_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Thrust / Weight")
# PyPlot.ylabel("Time (s)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("Time_Va")
# PyPlot.plot(Va_plot,time_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Va (m/s)")
# PyPlot.ylabel("Time (s)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("Time_Power")
# PyPlot.plot(Power_plot,time_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Power (w)")
# PyPlot.ylabel("Time (s)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("Time_Dist")
# PyPlot.plot(time_plot,dist_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Time (s)")
# PyPlot.ylabel("Distance (m)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))
#
# PyPlot.figure("TW_Vc")
# PyPlot.plot(TW_plot,height./time_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Thrust / Weight")
# PyPlot.ylabel("Climb Rate (m/s)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

# PyPlot.figure("Power_Vc")
# PyPlot.plot(PW_plot,height./time_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.xlabel("Power/Weight")
# PyPlot.ylabel("Climb Rate (m/s)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

PyPlot.figure("Net Power_Vc")
PyPlot.plot((PW_plot.*weight-D_perc_plot.*TW_plot.*weight.*Va_plot)/weight,height./time_plot,".-",label = "$(round(CL_array[j],3))")
PyPlot.xlabel("Net Power / Weight ")
PyPlot.ylabel("Climb Rate (m/s)")
legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

# PyPlot.figure("Energy3D")
# PyPlot.plot3D(TW_plot,energy_plot,dist_plot,".-",label = "$(round(CL_array[j],3))")
# PyPlot.zlabel("Climb Distance (m)")
# PyPlot.xlabel("Thrust / Weight")
# PyPlot.ylabel("Energy (J)")
# legend(loc="center left", title = "CL",bbox_to_anchor=(1, 0.5))

end
