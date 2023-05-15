using QuantumOptics
using Plots
plotlyjs()
##
m = 1.0
w = 1
xmin = -10.0
xmax = 10.0
Npoints = 400
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
function V_box(x)
    if abs(x) > 5
        return 100
    else
        return 0
    end
end
V = potentialoperator(b_position, V_box)
T_px = transform(b_momentum, b_position)
plot(LinRange(xmin,xmax,Npoints),V_box)
##
# Transforms a state multiplied from the right side from real space
# to momentum space.
T_xp = dagger(T_px)
x = position(b_position)
#p = momentum(b_momentum)
p = momentum(b_position)
#H_kin = LazyProduct(T_xp, p^2/(2m), T_px)
#H = LazySum(H_kin,V)
H = p^2/(2m) + dense(V)
##
x0 = -3#-1
p0 = 0#.6283185#1#4.5
sigma0 = 1
ψ0_data = [(sin((x+5)*pi/10)*(abs(x)<5) +0*im) for x in LinRange(xmin,xmax,Npoints)]
psi_0 =  normalize!((Ket(b_position,ψ0_data)))
#gaussianstate(b_position, x0, p0, sigma0) #+ gaussianstate(b_position, -x0, p0, sigma0)
E0 = abs(dagger(psi_0)*H*psi_0)
T = LinRange(0,60,400)
x_points = samplepoints(b_position)
plot(x_points,5*real(psi_0.data) .+E0)
plot!(x_points,V_box)
##
tout, psi_t = timeevolution.schroedinger(T, psi_0, H)
##
plot()
x_points = samplepoints(b_position)
anim=@animate for i in 1:200
    j=i*2
    n = 20*(abs.(psi_t[j].data).^2)./dx .+ E0
    avg_x = real(dagger(psi_t[j])*x*psi_t[j])
    plot(x_points,V_box,label=["V" "x"])
    plot!(x_points,n,label=["Ψ" "x" ])
    scatter!([avg_x],[E0],label=["<x>" "V" ],ylims=(-1,20))
end
##
gif(anim,"wavepacket_Box_E_n.gif",fps=20)
##
p_points = samplepoints(b_momentum)
anim=@animate for i in 1:400
    j=i#*2
    Psi_p = T_px*psi_t[j]
    n = 20*(abs2.(Psi_p.data))
    avg_p = real(dagger(Psi_p)*p*Psi_p)
    #plot(x_points,V_box,label=["V" "x"])
    plot!(p_points,n,label=["Ψp" "p" ])
    scatter!([avg_p],[E0],label=["<p>" "V" ])
end
##
gif(anim,"wavepacket_Box_E_p.gif",fps=20)
