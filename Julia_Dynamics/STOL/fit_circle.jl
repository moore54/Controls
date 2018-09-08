using PyPlot
PyPlot.close("all")
rc("figure", figsize=(6.5, 4.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.69)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])


pn = [5.20127, 10.402, 15.5576, 20.7141, 25.8523, 30.9586, 36.0118, 40.9922, 45.8837, 50.6867, 55.4229, 60.1316, 64.8596, 69.6499, 74.5318, 79.5136, 84.5813, 89.7062, 94.8557, 100.0]
pd = -[-0.00654977, -0.0453209, -0.177847, -0.409054, -0.887638, -1.63808, -2.68702, -4.03961, -5.68284, -7.56668, -9.60997, -11.7127, -13.7682, -15.6692, -17.3158, -18.6245, -19.5423, -20.0584, -20.2146, -20.0]
slope = (pd[2:end]-pd[1:end-1])./(pn[2:end]-pn[1:end-1])
aoc = atan.(slope)*180/pi
normal_angle = 90.0+aoc

PyPlot.figure()
PyPlot.plot(pn,pd,".-")
# axis("equal")
H = 100.0
for i = 1:length(normal_angle)
    normal_line_x = pn[i]+H*cos(normal_angle[i]*pi/180)
    normal_line_y = pd[i]+H*sin(normal_angle[i]*pi/180)

    PyPlot.plot([pn[i],normal_line_x],[pd[i],normal_line_y])
end

# circx = 25.0
# circy = -47.0
# circ_r = 46.0
#
# x = linspace(-circ_r+circx,circ_r+circx,100)
# y = sqrt.(circ_r.^2-(x-circx).^2)+circy
#
# figname = "fit_circle"
# PyPlot.figure()
# PyPlot.plot(pn,pd)
# PyPlot.plot(x,-y)
# # axis("equal")
