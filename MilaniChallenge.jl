push!(LOAD_PATH, pwd())

import Pkg
Pkg.activate("SALAMANDER")
Pkg.instantiate()

using DiffEqFlux
using Plots
#using OrdinaryDiffEq
using DifferentialEquations
using LinearAlgebra
using DiffEqSensitivity

#using ModelingToolkit, DataDrivenDiffEq
using BenchmarkTools
#using StaticArrays
using Random 
using Optim
using Flux
using TaylorIntegration

rng = MersenneTwister(1234)

#* Set Constants
#* Gravitational Parameters
GME = 3.986004407799724E5  #km^3s^-2
GMS = 1.32712440018E11     #km^3s^-2
GMM = 4.9028E3             #km^3s^-2

#* Earth's Gravity Field Parameters
RE = 6378.1363             #km
C20 = -4.84165371736E-4
C22 = 2.43914352398E-6
S22 = -1.40016683654E-6

#* Initial Angles
θG = (π/180) * 280.4606
νE = (π/180) * 4.178074622024230E-3
νS = (π/180) * 1.1407410259335311E-5
νMa = (π/180) * 1.512151961904581E-4
νMp = (π/180) * 1.2893925235125941E-6
νMs = (π/180) * 6.128913003523574E-7

#* Other Constants
aS = 1.49619E8             #km
ϵ = (π/180) * 23.4392911
ρS0 = (π/180) * 357.5256
ΩSωS = (π/180) * 282.94
PSRP = 4.56E-6

#* Initial Guess at AOM
AOM_guess = 1.0
### "True" AOM
AOM_true = .9658888019

#* Parameter Vector
p_ = Float64[GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP]
p_true = Float64[AOM_true]
p_guess = Float64[AOM_guess]

function Kepler(u, p, t, p_true)
    @inbounds begin
        X, Y, Z, Ẋ, Ẏ, Ż = u
        GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP = p_true
        AOM = p[1]

        R = √(X^2 + Y^2 + Z^2)

        return [-(GME*X)/(R^3); -(GME*Y)/(R^3); -(GME*Z)/(R^3)]
    end
end

function J2(u, p, t, p_true)
    @inbounds begin
        X, Y, Z, Ẋ, Ẏ, Ż = u
        GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP = p_true
        AOM = p[1]

        R = √(X^2 + Y^2 + Z^2)

        acc = zeros(Float64, 3)

        return [(GME*RE^2*√(5)*C20*X)/(2*R) * (3/(R^4) - (15*Z^2)/(R^6));
                (GME*RE^2*√(5)*C20*Y)/(2*R) * (3/(R^4) - (15*Z^2)/(R^6));
                (GME*RE^2*√(5)*C20*Z)/(2*R) * (9/(R^4) - (15*Z^2)/(R^6))]
    end
end

function C22S22(u, p, t, p_true)
    @inbounds begin
        X, Y, Z, Ẋ, Ẏ, Ż = u
        GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP = p_true
        AOM = p[1]

        R = √(X^2 + Y^2 + Z^2)
        
        x = X*cos(θG + νE*t) + Y*sin(θG + νE*t)
        y = -X*sin(θG + νE*t) + Y*cos(θG + νE*t)
        z = Z

        r = √(x^2 + y^2 + z^2)

        fC22x = (5*GME*RE^2*√(15)*C22*x*(y^2 - x^2))/(2*r^7) + (GME*RE^2*√(15)*C22*x)/(r^5)
        fC22y = (5*GME*RE^2*√(15)*C22*y*(y^2 - x^2))/(2*r^7) - (GME*RE^2*√(15)*C22*y)/(r^5)
        fC22z = (5*GME*RE^2*√(15)*C22*z*(y^2 - x^2))/(2*r^7)

        fS22x = -(5*GME*RE^2*√(15)*S22*x^2*y)/(r^7) + (GME*RE^2*√(15)*S22*y)/(r^5)
        fS22y = -(5*GME*RE^2*√(15)*S22*x*y^2)/(r^7) + (GME*RE^2*√(15)*S22*x)/(r^5)
        fS22z = -(5*GME*RE^2*√(15)*S22*x*y*z)/(r^7) 

        fC22X = fC22x*cos(θG + νE*t) - fC22y*sin(θG + νE*t)
        fC22Y = fC22x*sin(θG + νE*t) + fC22y*cos(θG + νE*t)
        fC22Z = fC22z

        fS22X = fS22x*cos(θG + νE*t) - fS22y*sin(θG + νE*t)
        fS22Y = fS22x*sin(θG + νE*t) + fS22y*cos(θG + νE*t)
        fS22Z = fS22z

        return [fC22X + fS22X;
                fC22Y + fS22Y;
                fC22Z + fS22Z]
    end
end

function SolarTide(u, p, t, p_true)
    @inbounds begin
        X, Y, Z, Ẋ, Ẏ, Ż = u
        GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP = p_true
        AOM = p[1]

        lS = ρS0 + νS*t
        rS = (149.619 - 2.499*cos(lS) - .021*cos(2*lS)) * 1E6
        λS = ΩSωS + lS +(π/180) * ((6892/3600)*sin(lS) + (72/3600)*sin(2*lS))

        XS = rS*cos(λS)
        YS = rS*sin(λS) * cos(ϵ)
        ZS = rS*sin(λS) * sin(ϵ)

        RS = √(XS^2 + YS^2 + ZS^2)
        RSsc = √((X-XS)^2 + (Y-YS)^2 + (Z-ZS)^2)

        acc = zeros(Float64, 3)

        return [-GMS * ((X-XS)/(RSsc^3) + XS/(RS^3));
                -GMS * ((Y-YS)/(RSsc^3) + YS/(RS^3));
                -GMS * ((Z-ZS)/(RSsc^3) + ZS/(RS^3))]
    end
end

function LunarTide(u, p, t, p_true)
    @inbounds begin
        X, Y, Z, Ẋ, Ẏ, Ż = u
        GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP = p_true
        AOM = p[1]

        ϕM = νS*t
        ϕMa =νMa * t
        ϕMp = νMp * t
        ϕMs = νMs * t

        L0 = ϕMp + ϕMa + (π/180)*218.31617
        lM = ϕMa + (π/180)*134.96292
        l′M = ϕM + (π/180)*357.5256
        FM = ϕMp + ϕMa + ϕMs + (π/180)*93.27283
        DM = ϕMp + ϕMa - ϕMs + (π/180)*297.85027

        rM = 385000 - 20905*cos(lM) - 3699*cos(2*DM - lM) -
             2956*cos(2*DM) - 570*cos(2*lM) +
             246*cos(2*lM - 2*DM) - 205*cos(l′M - 2*DM) -
             171*cos(lM + 2*DM) -
             152*cos(lM + l′M - 2*DM)
        
        λM = L0 + (π/180)*((22640/3600)*sin(lM) + (769/3600)*sin(2*lM)) -
             (4856/3600)*sin(lM - 2*DM) + (2370/3600)*sin(2*DM) -
             (668/3600)*sin(l′M) - (412/3600)*sin(2*FM) -
             (212/3600)*sin(2*lM - 2*DM) - (206/3600)*sin(lM + l′M - 2*DM) +
             (192/3600)*sin(lM + 2*DM) - (165/3600)*sin(l′M - 2*DM) +
             (148/3600)*sin(lM - l′M) - (125/3600)*sin(DM) -
             (110/3600)*sin(lM + l′M) - (55/3600)*sin(2*FM - 2*DM)

        βM = (π/180)*((18520/3600)*sin(FM + λM - L0 + (π/180)*((412/3600)*sin(2*FM) + (541/3600)*sin(l′M)))) -
             (526/3600)*sin(FM - 2*DM) +
             (44/3600)*sin(lM + FM - 2*DM) - (31/3600)*sin(-lM + FM -2*DM) -
             (25/3600)*sin(-2*lM + FM) - (23/3600)*sin(l′M + FM - 2*DM) +
             (21/3600)*sin(-lM + FM) + (11/3600)*sin(-l′M + FM - 2*DM)


        XM = rM*cos(λM)*cos(βM)
        YM = rM*sin(λM)*cos(βM)*cos(ϵ) - rM*sin(λM)*sin(ϵ)
        ZM = rM*sin(λM)*cos(βM)*sin(ϵ) + rM*sin(λM)*cos(ϵ)

        RM = √(XM^2 + YM^2 + ZM^2)
        RMsc = √((X-XM)^2 + (Y-YM)^2 + (Z-ZM)^2)

        return [-GMM * ((X-XM)/(RMsc^3) + XM/(RM^3));
                -GMM * ((Y-YM)/(RMsc^3) + YM/(RM^3));
                -GMM * ((Z-ZM)/(RMsc^3) + ZM/(RM^3))]
    end
end

function SRP(u, p, t, p_true)
    @inbounds begin
        X, Y, Z, Ẋ, Ẏ, Ż = u
        GME, GMS, GMM, RE, C20, C22, S22, θG, νE, νS, νMa, νMp, νMs, aS, ϵ, ρS0, ΩSωS, PSRP = p_true
        AOM = p[1]
+
        ZS = rS*sin(λS) * sin(ϵ)

        RS = √(XS^2 + YS^2 + ZS^2)
        RSsc = √((X-XS)^2 + (Y-YS)^2 + (Z-ZS)^2)

        return [AOM * (PSRP * aS^2 * (X-XS))/(RSsc^3);
                AOM * (PSRP * aS^2 * (Y-YS))/(RSsc^3);
                AOM * (PSRP * aS^2 * (Z-ZS))/(RSsc^3)]
    end
end

function EOM!(du, u, p, t, p_true)
    @inbounds begin
        du[1:3] .= u[4:6]
        du[4:6] .= Kepler(u, p, t, p_true) .+ J2(u, p, t, p_true) + C22S22(u, p, t, p_true) +
                   LunarTide(u, p, t, p_true) + SolarTide(u, p, t, p_true) + SRP(u, p, t, p_true)
    end                
end

@taylorize EOM_dyn!(du, u, p, t) = EOM!(du, u, p, t, p_)

#* Initial Condition
u0 = Float64[6800.; 0.0; 500.0; 0.0; √(GME/√(6800.0^2 + 520.0^2)); 0.0]

#* Generate Dataset over 30 Years
datasize = 200
T = 3600*24
#T = 3600*24
tspan = (0.0f0, T)
tsteps = range(tspan[1], tspan[2], length=datasize)

#* Setup ODE Problem
EOM_ode = ODEProblem(EOM_dyn!, u0, tspan, p_true)
# TODO: SWITCH TO TAYLOR INTEGRATOR
EOM_sol = solve(EOM_ode, TaylorMethod(20), p=p_true, abstol=1E-15)
ode_data = Array(EOM_sol)
#* Add Some Measurement Noise
noise = 1e-3 * randn(rng, Float64, size(ode_data))    
noise[4:6, :] = 1e-3 * noise[4:6, :]
noisy_data = ode_data + noise


#* Function to Predict Trajectory from Parameter Guess
function predict_univ(θ)
    return Array(solve(EOM_ode, Tsit5(), p=θ, abstol=1E-6, reltol=1E-6, saveat=EOM_sol.t, sensealg=ReverseDiffAdjoint()))
end

test = predict_univ(p_true)

#* Loss Function of the UDE
function loss_univ(θ)
    pred = predict_univ(θ)
    if size(pred)[2] == size(ode_data)[2]
        loss = sum(abs, ode_data .- pred)# + 1e-5 * sum(abs2, θ[2:end])
    else
        loss = Inf
    end
    return loss, pred
end

callback = function(p, l, pred; doplot=false)
    display(l)
    display(p)
    if doplot
        # Plot current prediction against data
        plt = Plots.scatter(noisy_data[1, :], noisy_data[2, :], label = "data")
        Plots.scatter!(plt, pred[1,:], pred[2,:], label = "prediction")

        plt2 = Plots.scatter(noisy_data[3, :], noisy_data[4, :], label = "data")
        Plots.scatter!(plt2, pred[3,:], pred[4,:], label = "prediction")

        plt3 = Plots.scatter(noisy_data[5, :], noisy_data[6, :], label = "data")
        Plots.scatter!(plt3, pred[5,:], pred[6,:], label = "prediction")
    
        display(plot(plt, plt2, plt3))
    end
    return false
end

result_univ = DiffEqFlux.sciml_train(loss_univ, p_guess,
                                     ADAM(), cb = callback,
                                     maxiters = 100)

result_univ2 = DiffEqFlux.sciml_train(loss_univ, result_univ.minimizer,
                                     LBFGS(), cb = callback,
                                     maxiters = 50)

println("AOM Error:", result_univ.minimizer - p_true)