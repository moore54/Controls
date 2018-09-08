using PyPlot
close("all")
x = linspace(-.05,.05,100)

gravity = 9.81
k_array = [500,1000,2000]
PyPlot.figure("test")
for i = 1:length(k_array)
k = k_array[i]
x0 = pi*e/k
lfun = gravity./(1+e.^(-k*(x-x0)))

PyPlot.plot(x,lfun,label = "$k")
end
PyPlot.legend(loc = "best")
# PyPlot.ylim(-.0010,.01)
PyPlot.xlim(-.05,.05)
