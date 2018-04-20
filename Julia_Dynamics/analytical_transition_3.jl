using PyPlot



m = 5.0 #kg
g = 9.81
L = m*g*2
D = L/30
T = m*g*.9
alpha = 15.0*pi/180
Va = 17.0
H_c = 20.0

gama = asin((T-D)/(m*g))
gamad = gama*180/pi

G = linspace(0,gama,100)-pi/2

r = m*Va^2/(L+T*sin(alpha)-m*g)

arclength = r*gama
t_tr = arclength/Va

xtr = r*cos.(G)
ytr = r*sin.(G)+r

H_left = H_c-ytr[end]
x_left = H_left/tan(gama)
xcl = [xtr[end],xtr[end]+x_left]
ycl = [ytr[end],H_c]

figure("position")
plot(xtr,ytr,label = "transition")
plot(xcl,ycl,label = "climb")
xlabel("x")
ylabel("y")
axis("equal")
