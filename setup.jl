# setup file, nice visualization can be found here: - https://visualpde.com/nonlinear-physics/cahn-hilliard

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
# Pkg.add("LaTeXStrings")
# Pkg.add("Random")
#Pkg.add("FFTW")
# Pkg.add("Plots")

using FFTW
using Plots
using LaTeXStrings
using Random

figdir = joinpath(@__DIR__, "Plots") #plotting path

# nicer plots 
 function plot_kwargs()
        fontsizes = (
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
            legendfontsize = 14,
        )
        (;
            #gridlinewidth = 2,            
            linewidth = 3,
            markersize = 6,
            markerstrokewidth = 0,
            fontsizes...,
            size = (600, 500),
        )
    end
#######################################################################
# 1D implementation
# implementation was tested on the examples in 4.2 in https://link.springer.com/article/10.1007/s10915-021-01471-6

function cahn_hilliard_1d(;N = 300 , L = 6.0 , dt = 0.001, T_end = 360.0, M = 1.0, kappa = 0.0225, phi0 = nothing )
    N = N
    L = L        
    dx = L / N

    dt = dt
    T_end = T_end
    nsteps = Int(round(T_end / dt))

    M = M
    kappa = kappa

    # grid
    x = [j * dx for j in 0:N-1]

    # sample frequencies for a DFT of length n
    k = fftfreq(N, N/L) * 2π  

    # derivatives in fourier space
    k2 = k.^2
    k4 = k.^4

    # denominator of update formula
    denom = 1 .+ dt .* M .* kappa .* k4

    # initial condition
    if phi0 === nothing
         u0 = cos.(2 .* x) .+ 0.01 .* exp.(cos.(x .+ 0.1))
    else
        u0 = phi0.(x)  
    end

    # initial condition in fourrier space

    phi_hat = fft(u0)

    # solving as in the script

    for _ in 1:nsteps

        # position space
        phi = real.(ifft(phi_hat))

        # nonlinear term in position space
        nonlinear = phi.^3 #.- phi

        # nonlinear term in fourier space
        nonlinear_hat = fft(nonlinear)

        # update in fourier space
        phi_hat = (phi_hat .- dt .* M .* k2 .* nonlinear_hat) ./ denom
    end

    # transform back in position space
    u = real.(ifft(phi_hat))
   return (; u0, u, x)

end

 result = cahn_hilliard_1d()
 x = result.x
 phi0 = result.u0
 u = result.u

# plot
p_1D = plot(x, phi0, label=L"IC:  \cos(2x)+0.01 \exp(\cos(x + 0.1))", xlabel = L"x", ylabel = L"\phi", title = "Cahn-Hilliard"; plot_kwargs()...);
plot!(p_1D, x, u, label=L"\phi(x,360)" ; plot_kwargs()...);
savefig(p_1D, joinpath(figdir, "default_1D.pdf"))

#######################################################################
#  2D implementation
# implementation was tested on the examples in 3.2 in https://www.mdpi.com/2227-7390/8/8/1385
function cahn_hilliard_2d(;
    Nx = 1000,
    Ny = 100,
    Lx = 1.0,
    Ly = 1.0,
    dt = 0.00001,
    T = 0.001,
    M = 1.0,
    epsilon = 0.0025,
    phi0 = nothing
)
    T = T
    dx = Lx / Nx
    dy = Ly / Ny
    Nt = Int(round(T / dt))
    kappa = epsilon^2

    x = range(-0.5Lx + dx, 0.5Lx, length=Nx)
    y = range(-0.5Ly + dy, 0.5Ly, length=Ny)

    # wavenumbers
    kx = fftfreq(Nx, Nx/Lx) * 2π
    ky = fftfreq(Ny, Ny/Ly) * 2π
    K2 = kx.^2 .+ (ky.^2)'
    K4 = K2.^2
    denom = 1 .+ dt .* M .* kappa .* K4

    # initial condition
    if phi0 === nothing
        Random.seed!(1234)
        u0 =  0.05 .* (2 .* rand(Nx, Ny) .- 1)
    else
        u0 = [phi0(xi, yi) for xi in x, yi in y]
    end

    # time evolution
    u_hat = fft(u0)

    for _ in 1:Nt
        u = real.(ifft(u_hat))
        nonlinear = u.^3 .- u
        nonlinear_hat = fft(nonlinear)
        u_hat = (u_hat .- dt .* M .* K2 .* nonlinear_hat) ./ denom
    end

    u = real.(ifft(u_hat))

    return (;u0, u , x, y, T)
end

# plotting
data = cahn_hilliard_2d()
T = data.T
x = data.x
y = data.y
u0 = data.u0
u = data.u
# plot
p1 = heatmap(x, y, u0',
    aspect_ratio=1,
    title="Initial condition",
    xlabel="x", ylabel="y"
);
p2 = heatmap(x, y, u',
    aspect_ratio=1,
    title="Cahn-Hilliard, T = $T",
    xlabel="x", ylabel="y"
);

p3 = plot(p1, p2, layout=(1,2), size=(1000,450));

savefig(p3, joinpath(figdir, "cahn_hilliard_2d.pdf"))

#######################################################################
# some animation 
# implementation was tested on the examples in 3.2 in https://www.mdpi.com/2227-7390/8/8/1385

function cahn_hilliard_2d_anim(;
    Nx = 1000,
    Ny = 1000,
    Lx = 1.0,
    Ly = 1.0,
    dt = 0.00001,
    T = 0.1,
    M = 1.0,
    epsilon = 0.0025,
    nframes = 100,
    fps = 20,
    filename = joinpath(figdir, "cahn_hilliard_2d.gif"),
    phi0 = nothing
)
    
    dx = Lx / Nx
    dy = Ly / Ny
    Nt = Int(round(T / dt))
    kappa = epsilon^2

    x = range(-0.5Lx + dx, 0.5Lx, length=Nx)
    y = range(-0.5Ly + dy, 0.5Ly, length=Ny)

    kx = fftfreq(Nx, Nx/Lx) * 2π
    ky = fftfreq(Ny, Ny/Ly) * 2π
    K2 = kx.^2 .+ (ky.^2)'
    K4 = K2.^2
    denom = 1 .+ dt .* M .* kappa .* K4

    if phi0 === nothing
        Random.seed!(1234)
        u0 = -0.45 .+ 0.05 .* (2 .* rand(Nx, Ny) .- 1)
    else
        u0 = [phi0(xi, yi) for xi in x, yi in y]
    end

    u_hat = fft(u0)
    anim = Animation()
    save_every = max(1, Nt ÷ nframes)

    for step in 1:Nt
        # update
        u = real.(ifft(u_hat))
        nonlinear_hat = fft(u.^3 .- u)
        u_hat = (u_hat .- dt .* M .* K2 .* nonlinear_hat) ./ denom

        # save frame
        if step % save_every == 0 || step == 1
            u_plot = real.(ifft(u_hat))
            heatmap(x, y, u_plot',
                aspect_ratio = 1,
                xlabel = "x", ylabel = "y",
                title = "t = $(round(step * dt, digits=5))",
                clims = (-1, 1),
                size = (600, 550)
            )
            frame(anim)
        end
    end

    gif(anim, filename, fps=fps)
end

cahn_hilliard_2d_anim()