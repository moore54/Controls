using PyPlot
using Dierckx
using QuadGK
close("all")

D = 1.0
# L = 10.0
m = 5.0 #kg
g = 9.81
LD = logspace(0,1.5,10)
T = 2.0
alpha = 10.0*pi/180

gama = linspace(-90*pi/180,90*pi/180,100)
gama_dot = []
gama_ddot = []

PyPlot.figure("PFA_dot")

for i = 1:length(LD)
    # D = D*1.1
    L = LD[i]*D
    gama_dot = atan2.((L*cos.(gama)-D*sin.(gama)+T*sin.(alpha+gama)-m*g),(-D*cos.(gama)-L*sin.(gama)+T*cos.(alpha+gama)))
    gama_spl = Dierckx.Spline1D(gama,gama_dot)
    gama_ddot = derivative(gama_spl,gama)

    plot(gama*180/pi,gama_dot,label = "L/D $(round.(L/D,0))")
    # plot(gama*180/pi,gama_ddot,label = "gama_ddot  L/D: $(round.(L/D,0))")

end
# plot(gama[end]*180/pi,gama_ddot[end],label = "gama_ddot  L/D: all")
xlabel("Flight Path Angle (deg)")
ylabel("Rate (rad/s or rad/s^2)")
legend(loc = "best")

D = 1.0
L = 20.0
m = 5.0 #kg
gama = 0.0
T = m*g*.9
alpha = 5.0*pi/180



function vy_i(gama)
    vy_i(m,L,D,gama,alpha,T)
end
function vy_i(m,L,D,gama,alpha,T)
    return (L*cos(gama)-D*sin(gama)+T*cos(alpha+gama)-m*g)/m
end

function vx_i(gama)
    vx_i(m,L,D,gama,alpha,T)
end
function vx_i(m,L,D,gama,alpha,T)
    return (-D*cos(gama)-L*sin(gama)+T*sin(alpha+gama))/m
end

# gama_dot0 = atan((-D+T*cos(alpha+gama))/(L+T*sin(alpha+gama)_)
dt = 0.01
t_total = 1.0
n_iters = t_total/dt
n_iters = round(Int, n_iters)
gama_save = zeros(n_iters)
t_save = zeros(n_iters)
Power = 25
x = 0.0
y = 0.0
Vx = 17.0
Vy = 0.0
t = 0.0
figure("path")
for i = 1:n_iters
    Va = sqrt(Vx^2+Vy^2)
    # T = Power/Va
    gama_dot = atan2.((L*cos.(gama)-D*sin.(gama)+T*sin.(alpha+gama)-m*g),(-D*cos.(gama)-L*sin.(gama)+T*cos.(alpha+gama)))
    gama = gama+gama_dot*dt
    Vx = Vx + vx_i(gama)*dt
    Vy = Vy + vy_i(gama)*dt
    x = x + Vx*dt
    y = y + Vy*dt
    t = t+dt
    figure("path")
    plot(x,y,"k.")

    # figure("gama")
    # plot(t,gama,"k.")
    gama_save[i] = gama
    t_save[i] = t

    figure("Va")
    plot(t,Va,"k.")

    # pause(0.0001)
end
# dx = QuadGK.quadgk(dx_i,0.0,-gama)
# dy = QuadGK.quadgk(dy_i,0.0,-gama)
figure("path")
axis("equal")
# plot(-dx[1]*t,-dy[1]*t,"r.",label = "Single Integration")
legend(loc = "best")
xlabel("x (m)")
ylabel("y (m)")


gama_spl = Dierckx.Spline1D(t_save,gama_save)
gama_ddot = derivative(gama_spl,t_save)

gamad_spl = Dierckx.Spline1D(t_save,gama_ddot)
gama_dddot = derivative(gamad_spl,t_save)

gamadd_spl = Dierckx.Spline1D(t_save,gama_dddot)
gama_ddddot = derivative(gamadd_spl,t_save)

figure("gama")
plot(t_save,gama_save,"r.")
plot(t_save,gama_ddot,"b.")
# plot(t_save,gama_dddot,"g.")
# plot(t_save,gama_ddddot,"y.")

#
# gama = linspace(0,90*pi/180,100)
# dy = zeros(gama)
# for i = 1:length(gama)
#     dy[i] = dy_i(gama[i])
# end
# figure("testy")
# xlabel("gama (deg)")
# ylabel("Vy_dot (m/s^2)")
# plot(gama*180/pi,dy)
#
# gama = linspace(0,90*pi/180,100)
# dx = zeros(gama)
# for i = 1:length(gama)
#     dx[i] = dx_i(gama[i])
# end
# figure("testx")
# xlabel("gama (deg)")
# ylabel("Vx_dot (m/s^2)")
# plot(gama*180/pi,dx)
