using QuantumOptics
using Plots
plotlyjs()
##
m = 1.0
w = 1
xmin = -20.0
xmax = 20.0
Npoints = 500
dx = (xmax-xmin)/Npoints
##
b_position = PositionBasis(xmin,xmax,Npoints)
# Hamiltonian in real space basis
p = momentum(b_position) # Dense operator
x = position(b_position) # Sparse operator
##
b_momentum = MomentumBasis(b_position)
##
V0 = 2
d = 8
n=3
function V_barrier(x)
    if x > -d/2 && x < d/2
        return -(abs(x)/2)^n + (d/4)^n
    elseif x < -15 || x >15
        return 20
    else
        return 0
    end
end
V = potentialoperator(b_position, V_barrier)
T_px = transform(b_momentum, b_position)
x_points = samplepoints(b_position)
plot(x_points,V_barrier)
##
# Transforms a state multiplied from the right side from real space
# to momentum space.
T_xp = dagger(T_px)

x = position(b_position)
p = momentum(b_momentum)

H_kin = LazyProduct(T_xp, p^2/2m, T_px)
H = LazySum(H_kin, V)
##
x0 = -10
p0 = 4.2
sigma0 = 1
psi_0 = gaussianstate(b_position, x0, p0, sigma0) #+ gaussianstate(b_position, -x0, p0, sigma0)
E0 = abs(dagger(psi_0)*H*psi_0)
T = LinRange(0,50,250)
x_points = samplepoints(b_position)
plot(x_points,(40*abs2.(psi_0.data).+E0))
plot!(x_points,V_barrier)
##
tout, psi_t = timeevolution.schroedinger(T, psi_0, H)
##
x_points = samplepoints(b_position)
plot()
anim=@animate for i in 1:250
    j=i
    n = (40*abs2.(psi_t[j].data)) .+ E0
    avg_x = real(dagger(psi_t[j])*x*psi_t[j])
    plot(x_points,V_barrier)
    plot!(x_points,n)
    scatter!([avg_x],[E0])
end
##
gif(anim,"wavepacket_SmoothBarrier_E=$(round(E0;digits=1)).gif",fps=20)
