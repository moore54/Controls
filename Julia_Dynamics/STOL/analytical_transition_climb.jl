using PyPlot
using Dierckx
using QuadGK
close("all")

D = 1.0
L = 30.0


gama = 45*pi/180
gamaddot = -1.0
t_total = (sqrt(4*gama*gamaddot+atan(-D/L)^2)-atan(-D/L))/(2*gamaddot)
m = 1 #kg
gama = gamaddot*t_total^2+atan(-D/L)*t_total

function dx_i(gama)
    dx_i(m,L,D,gama)
end
function dx_i(m,L,D,gama)
    return 1/m*(L*cos(gama)-D*sin(gama))
end

function dy_i(gama)
    dy_i(m,L,D,gama)
end
function dy_i(m,L,D,gama)
    return 1/m*(-D*cos(gama)-L*sin(gama))
end


dt = 0.01
n_iters = t_total/dt
round(Int, n_iters)
x = 0.0
y = 0.0
t = 0.0
figure("path")
for i = 1:n_iters
    gama = gamaddot*t^2+atan(-D/L)*t
    x = x + dx_i(-gama)*dt
    y = y + dy_i(-gama)*dt
    t = t+dt
    figure("path")
    axis("equal")
    plot(x,y,"k.")

    figure("gama")
    plot(t,gama*180/pi,"k.")

    # pause(0.0001)
end
dx = QuadGK.quadgk(dx_i,0.0,gama)
dy = QuadGK.quadgk(dy_i,0.0,gama)
figure("path")
axis("equal")
plot(dx[1]*t_total,dy[1]*t_total,"r.",label = "Single Integration")
legend(loc = "best")
xlabel("x (m)")
ylabel("y (m)")



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
