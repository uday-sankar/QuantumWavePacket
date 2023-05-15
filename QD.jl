using QuantumOptics
using Plots
plotlyjs()
##
m = 1.0
w = 1
xmin = -5.0
xmax = 5.0
Npoints = 100
##
b_position = PositionBasis(xmin,xmax,Npoints)
# Hamiltonian in real space basis
p = momentum(b_position) # Dense operator
x = position(b_position) # Sparse operator
H = p^2/2m + 1/2*m*w^2*dense(x^2)
##
b_momentum = MomentumBasis(b_position);
# Hamiltonian
p = momentum(b_momentum) # Sparse operator
x = position(b_momentum) # Dense operator
H = dense(p^2)/2m + 1/2*m*w^2*x^2
##
# Transforms a state multiplied from the right side from real space
# to momentum space.
T_px = transform(b_momentum, b_position)
##
T_xp = dagger(T_px)

x = position(b_position)
p = momentum(b_momentum)

H_kin = LazyProduct(T_xp, p^2/2m, T_px)
V = 1/2*w*x^2
H = LazySum(H_kin, V)
##
x0 = 1.5
p0 = 0
sigma0 = 2
psi_0 = gaussianstate(b_position, x0, p0, sigma0) #+ gaussianstate(b_position, x0, p0, sigma0)
E0 = abs(dagger(psi_0)*H*psi_0)
T = [0:0.025:10;]
tout, psi_t = timeevolution.schroedinger(T, psi_0, H)
##
x_points = samplepoints(b_position)
##
plot()
anim = @animate for i in 1:length(T)
    n = 40*abs.(psi_t[i].data).^2 .+ E0
    plot(x_points,0.5*w*x_points.^2)
    plot!(x_points,n)
end
##
gif(anim,"wavepacket.gif",fps=20)
##
x_points = samplepoints(b_position)
g=0
plot(x_points,w*x_points.^2)
for i in 1:length(T)
    n = 10*abs.(psi_t[i].data).^2 .+ 2
    plot!(x_points,n)
end
#plot!(x_points,w*x_points.^2)
savefig("gaussian.png")
##
plot(x_points,10*abs.(psi_t[end].data).^2)
##
