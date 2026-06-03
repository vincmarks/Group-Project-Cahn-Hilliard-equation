# nice visualization can be found here: - https://visualpde.com/nonlinear-physics/cahn-hilliard

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

figdir = joinpath(@__DIR__, "Plots") # plotting path

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
# 1D

function swift_hohenberg_1d(;
    N = 300,
    L = 6.0,
    q0 = 1.0,
    epsilon = 0.25,
    dt = 0.001,
    T_end = 20.0,
    M = 1.0,
    phi0 = nothing,
)
    dx = L / N
    nsteps = Int(round(T_end / dt))

    # grid
    x = [j * dx for j in 0:N-1]

    # wavenumbers for a periodic domain of length L
    k = fftfreq(N, N / L) .* 2π
    k2 = k.^2

    denom = 1 .+ dt .* M .* k2 .* ((q0^2 .- k2).^2 .- epsilon)

    # initial condition
    if phi0 === nothing
        Random.seed!(1234)

        # Small perturbation
        u0 = 0.05 .* rand(N) 
    else
        u0 = phi0.(x)
    end

    phi_hat = fft(u0)

    for _ in 1:nsteps

        phi = real.(ifft(phi_hat))

        nonlinear_hat = fft(phi.^3)

        phi_hat = (phi_hat .- dt .* M .* k2 .* nonlinear_hat) ./ denom
    end

    u = real.(ifft(phi_hat))
    return (; u0, u, x, epsilon, q0, T_end)
end

result = swift_hohenberg_1d()
x = result.x
phi0 = result.u0
u = result.u
T_end = result.T_end
p_1D = plot(
    x,
    phi0,
    label = L"IC: 0.05  \mathrm{rand}(x)",
    xlabel = L"x",
    ylabel = L"\phi",
    title = "Swift-Hohenberg 1D";
    plot_kwargs()...,
)
plot!(p_1D, x, u, label = L"\phi(x,t = 20)" ; plot_kwargs()...)
savefig(p_1D, joinpath(figdir, "1d_default.pdf"))

#######################################################################
# 2D

function swift_hohenberg_2d(;
    Nx = 400,
    Ny = 400,
    Lx = 60.0,
    Ly = 60.0,
    dt = 0.1,
    T = 50.0,
    M = 1.0,
    epsilon = 0.25,
    q0 = 1.0,
    phi0 = nothing,
)
    dx = Lx / Nx
    dy = Ly / Ny
    Nt = Int(round(T / dt))

    x = range(-0.5Lx + dx, 0.5Lx, length = Nx)
    y = range(-0.5Ly + dy, 0.5Ly, length = Ny)

    # wavenumbers
    kx = fftfreq(Nx, Nx / Lx) .* 2π
    ky = fftfreq(Ny, Ny / Ly) .* 2π

    K2 = kx.^2 .+ (ky.^2)'

    # denominator of update formula
    denom = 1 .+ dt .* M .* K2 .* ((q0^2 .- K2).^2 .- epsilon)

    if phi0 === nothing
        Random.seed!(1234)

        # Small perturbation
        u0 = 0.05 .* rand(Nx, Ny) 
    else
        u0 = [phi0(xi, yi) for xi in x, yi in y]
    end

    u_hat = fft(u0)

    for _ in 1:Nt

        u = real.(ifft(u_hat))

        nonlinear_hat = fft(u.^3)

        u_hat = (u_hat .- dt .* M .* K2 .* nonlinear_hat) ./ denom
    end

    u = real.(ifft(u_hat))
    return (; u0, u, x, y, T, epsilon, q0)
end

data = swift_hohenberg_2d()
T = data.T
x = data.x
y = data.y
u0 = data.u0
u = data.u

p1 = heatmap(
    x,
    y,
    u0',
    aspect_ratio = 1,
    title = "Initial condition",
    xlabel = "x",
    ylabel = "y",
    c = :viridis,
);

p2 = heatmap(
    x,
    y,
    u',
    aspect_ratio = 1,
    title = "Final state, t = $T",
    xlabel = "x",
    ylabel = "y",
    c = :viridis,
);

p3 = plot(p1, p2, layout = (1, 2), size = (1200, 550), dpi = 300)
savefig(p3, joinpath(figdir, "2d_default.png"))

#######################################################################
# 2D animation

function swift_hohenberg_2d_anim(;
    Nx = 400 ,
    Ny = 400 ,
    Lx = 60.0 ,
    Ly = 60.0 ,
    dt = 0.1,
    T = 50.0,
    M = 1.0,
    epsilon = 0.25,
    q0 = 1.0,
    nframes = 500,
    phi0 = nothing,
)
    dx = Lx / Nx
    dy = Ly / Ny
    Nt = Int(round(T / dt))

    x = range(-0.5Lx + dx, 0.5Lx, length = Nx)
    y = range(-0.5Ly + dy, 0.5Ly, length = Ny)

    kx = fftfreq(Nx, Nx / Lx) .* 2π
    ky = fftfreq(Ny, Ny / Ly) .* 2π
    K2 = kx.^2 .+ (ky.^2)'

    denom = 1 .+ dt .* M .* K2 .* ((q0^2 .- K2).^2 .- epsilon)

    if phi0 === nothing
        Random.seed!(1234)
        u0 = 0.05 .* rand(Nx, Ny)
    else
        u0 = [phi0(xi, yi) for xi in x, yi in y]
    end

    u_hat = fft(u0)

    anim = Animation()
    save_every = max(1, Nt ÷ nframes)

    # include the initial condition
    plt = heatmap(
        x, y, u0',
        aspect_ratio = 1,
        xlabel = "x", ylabel = "y",
        title = "PFC: t = 0.0",
        size = (800, 700),
        c = :viridis,
    )
    frame(anim, plt)

    for step in 1:Nt
        u = real.(ifft(u_hat))
        nonlinear_hat = fft(u.^3)
        u_hat = (u_hat .- dt .* M .* K2 .* nonlinear_hat) ./ denom

        if step % save_every == 0 
            u_plot = real.(ifft(u_hat))
            plt = heatmap(
                x, y, u_plot',
                aspect_ratio = 1,
                xlabel = "x", ylabel = "y",
                title = "Swift-Hohenberg: t = $(round(step * dt, digits = 5))",
                size = (800, 700),
                c = :viridis,
            )
            frame(anim, plt)
        end
    end

    return anim
end

anim = swift_hohenberg_2d_anim()
fps = 20
filename = joinpath(figdir, "default_2d.gif")
gif(anim, filename, fps = fps)