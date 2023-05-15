using QuantumOptics
using Plots 
plotlyjs()
##
m = 1.0
w = 1
xmin = -20.0
xmax = 20.0
Npoints = 400
dx = (xmax-xmin)/Npoints
##Setting up the basis
b_position = PositionBasis(xmin,xmax,Npoints)
b_momentum = MomentumBasis(b_position)
# Hamiltonian in real space basis
x = position(b_position) # Sparse operator
p = momentum(b_momentum) # Dense operator
T_px = transform(b_momentum, b_position)
T_xp = dagger(T_px)
##
V0 = 2
d = 10
n=3
function V_barrier(x)
    if x > -d/2 && x < d/2
        return 10#-(abs(x)/2)^n + (d/4)^n
    elseif x < -15
        return 20
    elseif x > 15
        return -im*(x-15)*10/5
    else
        return 0
    end
end
V = potentialoperator(b_position, V_barrier)
x_points = samplepoints(b_position)
x_ax = -20:0.1:15
plot(x_ax,V_barrier)
##
# Transforms a state multiplied from the right side from real space
# to momentum space.
H_kin = LazyProduct(T_xp, p^2/2m, T_px)
H = LazySum(H_kin, V)
##
x0 =-12
p0 = 3.8
sigma0 = 1
k=+4.5
#psi_0 = gaussianstate(b_position, x0, p0, sigma0) #+ gaussianstate(b_position, -x0, p0, sigma0)
#psi_0 = normalize!(Ket(b_position,[(exp(-(x+10.0)^2) + 0*im)*exp(-k*im*x) for x in LinRange(xmin,xmax,Npoints)]))
psi_0_dat = []
for x in LinRange(xmin,xmax,Npoints)
    if x > -15 && x < -5
        append!(psi_0_dat,(sin(2pi*(x-15)/20) + 0*im)*exp(-k*im*x))
    else
        append!(psi_0_dat,0*im)
    end
end
psi_0_dat = Vector{ComplexF64}(psi_0_dat)
psi_0 = normalize!(Ket(b_position,psi_0_dat))
E0 = abs(dagger(psi_0)*H*psi_0)
##
T = LinRange(0,20,120)
x_points = samplepoints(b_position)
plot(x_points,(15*abs2.(psi_0.data).+E0))
plot!(x_ax,V_barrier)
##
tout, psi_t = timeevolution.schroedinger(T, psi_0, H)
##
x_points = samplepoints(b_position)
plot()
anim=@animate for i in 1:length(T)
    j=i
    ψn = psi_t[j]
    En = real(dagger(ψn)*H*ψn)
    n = (30*abs2.(ψn.data)) .+ En
    avg_x = real(dagger(psi_t[j])*x*psi_t[j])*sum(abs2.(ψn.data)) + 15*(1-sum(abs2.(ψn.data)))
    plot(x_ax,V_barrier)
    plot!(x_points,n)
    scatter!([avg_x],[En])
end
##
gif(anim,"Animation/StepBarrier_userdefined_function=$(round(E0;digits=1)).gif",fps=20)
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