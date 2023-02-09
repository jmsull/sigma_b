#jsull
using DelimitedFiles
using Interpolations

# Load P(k) and (log) interpolate 
ret = open("./planck15_x0.dat","r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
class_pk = reduce(hcat,ret)
itpclass = LinearInterpolation(log10.(class_pk[1,:]),log10.(class_pk[2,:]),
        extrapolation_bc=Line())
P(k) = 10^(itpclass(log10(k)))

# The integrand for doing log integral
function logI(lk1,lk2,lk3,L)
    ω = L/2
    k1,k2,k3 = 10^lk1,10^lk2,10^lk3
    k = sqrt(k1^2 + k2^2 + k3^2)
    W =  sin(ω*k1) * sin(ω*k2) * sin(ω*k3) 
    return 1.0/(2π)^3 * abs(W)^2 / (k1*k2*k3*ω^6) * P(k) 
end

# The window function integral on its own
function logWsq(lk1,lk2,lk3,L)
    ω = L/2
    k1,k2,k3 = 10^lk1,10^lk2,10^lk3
    k = sqrt(k1^2 + k2^2 + k3^2)
    W =  sin(ω*k1) * sin(ω*k2) * sin(ω*k3) 
    return 1.0/(2π)^3 * abs(W)^2 / (k1*k2*k3*ω^6)   
end


# Set up inntegral
L = 625.0 #Box side length in Mpc/h
Vbox = L^3
kmin=5e-6
kmax=0.5/(√3)
N=128
logkx=range(log10(kmin),log10(kmax),length=N)
kx = 10 .^(logkx)

Integrand = [logI(i,j,k,L) for i=logkx,j=logkx,k=logkx];
Window = [logWsq(i,j,k,L) for i=logkx,j=logkx,k=logkx];

#very simple integration routine, surely will be slow 
function trapz_3d(x,y,z,f)
    nx,ny,nz = length(x),length(y),length(z)
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1] #this is dlog10kx
    I=0.0
    for i in 1:nx
        wx = dx/(1.0 + (i==1) + (i==nx) )
        for j in 1:ny
            wy = dy/(1.0 + (j==1) + (j==ny))
            for k in 1:nz
                wz = dz/(1.0 + (k==1) + (k==nz))
                I += f[i,j,k] *wx*wy*wz
            end
        end
    end
    return I
end

println("Is the window close to 1? ", 8trapz_3d(logkx,logkx,logkx,Window)* log(10)^3 *Vbox - 1.0 < 0.03)

σbsq = 8trapz_3d(logkx,logkx,logkx,Integrand)* log(10)^3 #factor of 8 for each octant, log factor for logarithmic integration measure

println("Power spectrum variance in a cubic box of side length $(L) Mpc/h is $(σbsq).")

