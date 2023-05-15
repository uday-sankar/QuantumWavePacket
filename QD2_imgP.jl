using QuantumOptics
using Plots
plotlyjs()
##
m = 1.0
w = 1
xmin = -20.0
xmax = 20.0
Npoints = 200
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
    elseif x < -15
        return 20
    elseif x > 15
        return -im*(x-15)*10/5
    else
        return 0
    end
end
V = potentialoperator(b_position, V_barrier)
T_px = transform(b_momentum, b_position)
x_points = samplepoints(b_position)
x_ax = -20:0.1:15
plot(x_ax,V_barrier)
##
# Transforms a state multiplied from the right side from real space
# to momentum space.
T_xp = dagger(T_px)

x = position(b_position)
p = momentum(b_momentum)

H_kin = LazyProduct(T_xp, p^2/2m, T_px)
H = LazySum(H_kin, V)
##
x0 =-12
p0 = 3.8
sigma0 = 1
psi_0 = gaussianstate(b_position, x0, p0, sigma0) #+ gaussianstate(b_position, -x0, p0, sigma0)
E0 = abs(dagger(psi_0)*H*psi_0)
T = LinRange(0,20,120)
x_points = samplepoints(b_position)
plot(x_points,(15*abs2.(psi_0.data).+E0))
plot!(x_ax,V_barrier)
##
tout, psi_t = timeevolution.schroedinger(T, psi_0, H)
##
Cros_prob = zeros(length(T))
Et = ones(length(T))
for i in 1:length(T)
    ψn = psi_t[i]
    Et[i] = abs(dagger(ψn)*H*ψn)
    tot_pro = sum(abs2.(ψn.data))
    Cros_prob[i] = 1 - tot_pro
end
Cros_prob
##
plot(T,Cros_prob)
plot(T,Et)
##
x_points = samplepoints(b_position)
plot()
anim=@animate for i in 1:length(T)
    j=i
    ψn = psi_t[j]
    En = real(dagger(ψn)*H*ψn)
    n = (30*abs2.(ψn.data)) .+ En
    avg_x = real(dagger(psi_t[j])*x*psi_t[j])*sum(abs2.(ψn.data)) + 15*(1-sum(abs2.(ψn.data)))
    plot(x_ax,V_barrier,dpi=650)
    plot!(x_points,n,dpi=650)
    scatter!([avg_x],[E0],dpi=650)
end
##
gif(anim,"Animation/SmoothBarrier_imgWell=$(round(E0;digits=1)).mp4",fps=20)
##
using PlotlyJS
config = PlotConfig(
    toImageButtonOptions=attr(
        format="png", # one of png, svg, jpeg, webp
        filename="test",
        height=800,
        width=1200,
        scale=1 # Multiply title/legend/axis/canvas sizes by this factor
    ).fields
)
##
pl=plot(x_ax,V_barrier,)
png(pl, "test.png")