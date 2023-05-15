using QuantumOptics
using Plots
plotlyjs()
##
m = 1.0
w = 1
xmin = -18.0
xmax = 18.0
Npoints = 500
dx = (xmax-xmin)/Npoints
##
b_position = PositionBasis(xmin,xmax,Npoints)
# Hamiltonian in real space basis
p = momentum(b_position) # Dense operator
x = position(b_position) # Sparse operator
##
b_momentum = MomentumBasis(b_position);
##
V0 = 2
d = 8
V_DobWell(x)= (x^4-200*x^2)*1e-3+10
V = potentialoperator(b_position, V_DobWell)
T_px = transform(b_momentum, b_position)
plot(LinRange(xmin,xmax,Npoints),V_DobWell)
##
# Transforms a state multiplied from the right side from real space
# to momentum space.
T_xp = dagger(T_px)

x = position(b_position)
p = momentum(b_momentum)
H_kin = LazyProduct(T_xp, p^2/(2m), T_px)
const ϵ = LazySum(H_kin,V)
function H(t,psi)
    H_kin = LazyProduct(T_xp, p^2/(2m), T_px)
    return  LazySum(H_kin,V)
end
##
x0 = -10#-1
p0 = 4.4#1#4.5
sigma0 = 1
psi_0 = gaussianstate(b_position, x0, p0, sigma0) #+ gaussianstate(b_position, -x0, p0, sigma0)
E0 = abs(dagger(psi_0)*H(0,psi_0)*psi_0)
T = LinRange(0,50,100)
x_points = samplepoints(b_position)
plot(x_points,20*abs2.(psi_0.data).+E0)
plot!(x_points,V_DobWell)
##
tout, psi_t = timeevolution.schroedinger_dynamic(T, psi_0, H)
##
x_points = samplepoints(b_position)
plot()
anim=@animate for i in 1:100
    j=i#*2
    n = 20*(abs.(psi_t[j].data).^2)./dx .+ E0
    avg_x = real(dagger(psi_t[j])*x*psi_t[j])
    plot(x_points,V_DobWell)
    plot!(x_points,n)
    scatter!([avg_x],[E0])
end
##
gif(anim,"wavepacket_DW_E=$(round(E0;digits=1)).gif",fps=20)
