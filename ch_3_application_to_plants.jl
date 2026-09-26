



# ===============================================================================
# 1. Make a graph with temperature (°C) as the ordinate (y-axis) and
# radiation intensity (W m-2) as the abscissa (x-axis) and plot the blackbody
# radiation law, R = σ(Tₑ + 273.15)^4 The temperature scale should go from 0° to 60°C and the radiation scale
#from 300 to 700 W/m^2.

using Plots
using Unitful

σ = 5.67e-8u"W/m^2/K^4"  # Stefan-Boltzmann constant (W m^-2 K^-4)
T_e = [0:1:60;]u"°C" # No units because the next line 
R = σ * (uconvert.(u"K", T_e)).^4  # Blackbody radiation law (W/m^2)

# this is the same as applying the radiation law I just converted to Kelvin inside 
# instead of adding 273.15u"K" manually.

plot(R, T_e, xlabel = "Radiation Intensity", 
            ylabel = "Temperature", 
            title="Blackbody Radiation Law", 
            xlims = (300, 700),
            ylims = (0, 60))


# ===============================================================================
# 2. The rate at which energy C is transferred by convection to or from a
# flat plate (such as a broad leaf) and the surrounding air is given by the
# equation below:

C = h_c * (T_l - T_a) = k1 * (V^0.5 / D^0.5) * (T_l - T_a)

# where he is the convection coefficient, given by

h_c = k1 * (V^0.5 / D^0.5)  # Convection coefficient (W/m^2/K)

# Compute a table of convection coefficients (in watts per square meter per
# degree Centigrade) for the following values of V and D: 
k1 = 0.664 # for smooth plates with Re < 5*10^5

velocities = [0.1, 1, 5] # m/s
diameters = [0.01, 0.05, 0.1, 0.3] # m

h_c_table = k1 * (velocities.^0.5' / diameters.^0.5)u"W/m^2/K" 

# The amount of energy going into convection can be estimated by multiplying
# these values by a temperature differential between the leaf and air. Assume
# that this difference is 5°C and determine the amount of energy consumed
# by convection if D = 0.05 m and V = 1.0 m/s.

C = h_c_table[2, 2] * 5u"K"  # Select the 0.05m, 1 m/s entry and apply 5°C difference
C_watts = uconvert(u"W/m^2", C)







# ===============================================================================
# 3. The rate at which water vapor diffuses from a leaf is given by ....
# too much to write read p54

Ta = 30u"°C"
da =  0.0304u"kg/m^3" # the saturation density of water vapor in the air at 30°C 
dl =  0.03038u"kg/m^3" # the saturation density of water vapor in the leaf at 30°C
hs = [0.1, 0.5, 0.9]
# water vapour density of the air = h * da  at 30°C 
rlra = 200u"s/m"

water_vapour_density_air = hs .* da
water_vapour_density_leaf = dl

E = (water_vapour_density_leaf .- water_vapour_density_air) / (rlra) # rate of water vapour diffusion from a leaf 
E = uconvert.(u"kg/m^2/s", E)

λ = 2.430e6u"J/kg"

LE = λ * E  # latent heat flux due to water vapor diffusion
LE = uconvert.(u"W/m^2", LE)

# This gives the answer for the effect of humidity just at 30C, but the question actually 
# asks for a much more complicated thing which needs the functions created below 
# to calculate water vapour densities at different temperatures... thanks Gates. 

# Start by defining constants and relative humidities for the calculations
Rv = 461.5u"J/(kg*K)"          # gas constant for water vapour
Tair = (30.0 + 273.15)u"K"     # air temperature
rtotal = 200.0u"s/m"           # leaf + aerodynamic resistance

relative_humidities = [0.1, 0.5, 0.9]

# Saturated water-vapour density
function saturation_density(T)
    # Convert Kelvin temperature to a numerical Celsius value
    Tc = ustrip(u"K", T) - 273.15

    # Saturation vapour pressure (Tetens equation)
    es = 0.61078 *
         exp(17.2694 * Tc / (Tc + 237.3)) *
         u"kPa"

    # Ideal-gas relationship: density = pressure / (Rv T)
    ρsat = es / (Rv * uconvert(u"K", T))

    return uconvert(u"g/m^3", ρsat)
end

# ------------------------------------------------------------
# Latent heat of vaporisation
# ------------------------------------------------------------

function latent_heat(T)
    Tc = ustrip(u"K", T) - 273.15

    # Temperature-dependent latent heat
    λ = (2.501e6 - 2361.0 * Tc) * u"J/kg"

    return λ
end

# ------------------------------------------------------------
# Evaporation and latent energy flux
# ------------------------------------------------------------ 

function evaporation_rate(Tleaf, RH)
    leaf_vapour = saturation_density(Tleaf)
    air_vapour = RH * saturation_density(Tair)

    E = (leaf_vapour - air_vapour) / rtotal

    return uconvert(u"kg/(m^2*s)", E)
end

function latent_energy_flux(Tleaf, RH)
    E = evaporation_rate(Tleaf, RH)
    λ = latent_heat(Tleaf)

    return uconvert(u"W/m^2", λ * E)
end



# Now plotting the graph of leaf temperatures versus latent energy flux 
# i.e. water vapour diffusion from the leaf

using Plots
using Unitful

leaf_temps_C = collect(20:1:40)u"°C"
leaf_temps_K = uconvert.(u"K", leaf_temps_C)

p = plot(
    xlabel = "Latent Energy Flux λE (W/m^2)",
    ylabel = "Leaf Temperature (°C)",
    title = "Leaf Temperature vs Latent Energy Flux at Tair = 30°C",
    legend = :bottomright,
    linewidth = 2
)


for RH in relative_humidities
    LE = latent_energy_flux.(
        leaf_temps_K, 
        Ref(RH)
    )

    LE_values = ustrip.(u"W/m^2", LE)

    plot!(
        p,
        LE_values, 
        leaf_temps_C,
        label = "$(Int(RH * 100))%"
    )
end

vline!(p, [0.0], colour = :black, linestyle = :dash, label = "Zero Flux")

display(p)








# ===============================================================================
# 4. Plot leaf temperature on the y axis and the Radiation + Convection + Evaporative 
# on the x axis for some specified conditions 
# 
# Okay so here I flipped the X and Y because it was super hard to understand the way Gates
# asks for it so keep that in mind when it comes time to plot
# ===============================================================================
using Plots
using Unitful

D = 0.05u"m"
rtotal = 200.0u"s/m"
V = 1.0u"m/s"
Tair = 30u"°C" #just for the latent_energy_flux function which takes Tair
h = 0.5
ϵ = 0.96 # Emissivity of the leaf
σ = 5.67e-8u"W/m^2/K^4" # Stefan Boltzmann constant 
k1 = 9.14u"J/m^2/s^0.5/K" 

T_ls = collect(20:0.5:45) .* u"°C"



R = uconvert.(u"W/m^2", ϵ * σ .* uconvert.(u"K", T_ls) .^ 4)  # W/m^2
C = uconvert.(u"W/m^2", k1 * sqrt(V / D) .* (uconvert.(u"K", T_ls) .- uconvert(u"K", Tair)))  # W/m^2
E = latent_energy_flux.(uconvert.(u"K", T_ls), h) # Defined in previous question 

RCE = R .+ C .+ E

# Plotting total energy loss (Radiation + Convection + Evaporation) at different leaf temperatures
p = plot(
    T_ls,
    RCE,
    xlabel = "Leaf temp",
    ylabel = "Total Energy Loss",
    title = "Leaf Temperature vs Total Energy Loss \nTair = 30°C, rh = 50%, V = 0.1m/s, D = 0.05m",
    legend = :bottomright,
    linewidth = 2,
)

plot!(p, T_ls, R, label = "Radiation")
plot!(p, T_ls, C, label = "Convection")
plot!(p, T_ls, E, label = "Evaporation")

display(p)


# Using this graph, determine the leaf temperature when the amount of
# radiation absorbed by the leaf is Qa = 700 W m-2
# Remembering here that the energy budget is R + C + E = Qa
# So if we're getting 700 W/m^2 in then RCE must equal that so select 
# for the leaf temp where that condition is met

target_Qa = 700.0u"W/m^2"
i = argmin(abs.(RCE .- target_Qa))
leaf_temp_at_target_Qa = T_ls[i] 
println("Leaf temperature at target when radiation absorbed is Qa = 700 W/m^2: $leaf_temp_at_target_Qa")


# What percentage of the total radiation absorbed is reradiated (R), 
# lost to convection (C), and goes into the evaporation of water (AB)? Presumable at Qa = 700 W/m^2

prop_R = R[i] / RCE[i]
prop_C = C[i] / RCE[i]
prop_E = E[i] / RCE[i]

println("Proportion of total energy loss at Qa = 700 W/m^2:
Radiation: $(round(prop_R * 100, digits=1))%
Convection: $(round(prop_C * 100, digits=1))%
Evaporation: $(round(prop_E * 100, digits=1))%")


# What is the rate of water loss in kilograms per square meter per second?
# λ is the latent heat of vaporization of water
λ = 2.45e6u"J/kg" # J/kg
water_loss_rate = E[i] / latent_heat(u"30°C")
water_loss_rate = uconvert(u"kg/m^2/s", water_loss_rate)
println("Rate of water loss at Qa = 700 W/m^2: $water_loss_rate")



# 5. For an amount of absorbed radiation of 500 W m-2 , determine the
# leaf temperature. Again, what percentage of the total radiation absorbed
# is reradiated (R), lost to convection (C) and goes into the evaporation of
# water (AB)?

target_Qa = 500.0u"W/m^2"
i = argmin(abs.(RCE .- target_Qa))
leaf_temp_at_target_Qa = T_ls[i] 
println("Leaf temperature at target when radiation absorbed is Qa = 500 W/m^2: $leaf_temp_at_target_Qa")
prop_R = R[i] / RCE[i]
prop_C = C[i] / RCE[i]
prop_E = E[i] / RCE[i]
println("Proportion of total energy loss at Qa = 500 W/m^2:
Radiation: $(round(prop_R * 100, digits=1))%
Convection: $(round(prop_C * 100, digits=1))%
Evaporation: $(round(prop_E * 100, digits=1))%")





# 6. For a relative humidity of 0.1 (10%) and all other quantities the
# same as in Problem 4, determine the leaf temperature when Qa = 700 W/m^2

# Keeping all the same objects and methods from earlier just changing h here
h = 0.1 # relative humidity of 10%

# Only need to recalculate the evaporation term E based on the new relative humidity h
E = latent_energy_flux.(uconvert.(u"K", T_ls), h) # Defined in previous question 

RCE = R .+ C .+ E

# Plotting total energy loss (Radiation + Convection + Evaporation) at different leaf temperatures
p = plot(
    T_ls,
    RCE,
    xlabel = "Leaf temp",
    ylabel = "Total Energy Loss",
    title = "Leaf Temperature vs Total Energy Loss \nTair = 30°C, rh = 10%, V = 0.1m/s, D = 0.05m",
    legend = :bottomright,
    linewidth = 2,
)

plot!(p, T_ls, R, label = "Radiation")
plot!(p, T_ls, C, label = "Convection")
plot!(p, T_ls, E, label = "Evaporation")

display(p)

# Now can use the minarg to find the leaf temperature corresponding to 
# the minimum temp where Qa = 700 W/m^2
target_Qa = 700.0u"W/m^2"
i = argmin(abs.(RCE .- target_Qa))
leaf_temp_at_target_Qa = T_ls[i] 
println("Leaf temperature at target when radiation absorbed is Qa = 700 W/m^2: $leaf_temp_at_target_Qa")


# How are the percentages of R, C, and λE changed?
prop_R = R[i] / RCE[i]
prop_C = C[i] / RCE[i]
prop_E = E[i] / RCE[i]

println("Proportion of total energy loss at Qa = 700 W/m^2:
Radiation: $(round(prop_R * 100, digits = 1))%
Convection: $(round(prop_C * 100, digits = 1))%
Evaporation: $(round(prop_E * 100, digits = 1))%")


# What is the rate of water loss in kilograms per square meter per
# second?
water_loss_rate = E[i] / latent_heat(30u"°C")
water_loss_rate = uconvert(u"kg/m^2/s", water_loss_rate)
println("Water loss rate at Qa = 700 W/m^2: $water_loss_rate")




# 7. If W = m D, substitution into Eq. (3.5) gives
#
#    r_a = k_2 * m^0.2 * (D / V)^(1/2)
#
# Determine values of the boundary layer resistance in seconds per meter
# for the following values of m and V/D:
#
#                         m
# V/D      0.1 s m^-1     1.0 s m^-1     10.0 s m^-1
# -------------------------------------------------
#   1
#  10
# 100


# Gates is tricky here because m is actually dimensionless and he gives 
# V/D instead of D/V so you gotta divide by sqrt(V/D) in the equation. 

m = [0.1, 1.0, 10.0] # values of m should actually be dimensionless
V_D = [1, 10, 100]u"s^-1" # values of V/D should be in s^-1
k2 = 200 * u"s^0.5/m"

ra = round.(typeof(1u"s/m"), (k2 .* (m' .^0.2) ./ sqrt.(V_D)))
ra = uconvert.(u"s/m", ra)

using PrettyTables
pretty_table(ra, 
    column_labels = ["m = 0.1", "m = 1", "m = 10"],
    row_labels = ["V/D = 1", "V/D = 10", "V/D = 100"],
    title = "Boundary Layer Resistance (s/m)"
    )

# So you can see here the boundary layer resistance correctly grows with m, 
# and decreases with increasing V/D i.e. wind speed. 
# To understand this properly you have to know that 
# increasing m does not necessarily mean a shorter leaf. It only describes the leaf’s proportions as fatter.






# 8. Determine the leaf temperature and rate of water loss for a leaf with
# dimensions W and D = 0.05 m and internal resistance re = 100 s m-1 if air
# temperature is 30°C, wind speed is 0.1 m/s, relative humidity is 0.5
# (50%), and amount of radiation absorbed is 700 W m-2 

W = 0.05u"m"
D = 0.05u"m"
re = 100u"s/m"
T_air = 30u"°C"
V = 0.1u"m/s"
RH = 0.5
Qa = 700u"W/m^2"
k2 = 200u"s^0.5/m"
m = W/D # remember this aspect ratio is dimensionless 


# Start by calculating the boundary layer resistance based on the leaf dimensions and wind speed.
ra = round(typeof(1u"s/m"), (k2 * (m ^ 0.2) / sqrt.(V / D)))

# Get the total resistance by adding the internal resistance to the boundary layer resistance.
rtotal = round(typeof(1u"s/m"), (ra + re))
rtotal = 400u"s/m"
println("Total resistance: $rtotal
Internal resistance: $re
Boundary layer resistance: $ra")

# Now you can use the total resistance to calculate the leaf temperature and rate of water loss.
# Staring with the C and R pathways can still use T_ls

T_ls = collect(20:1:50)u"°C"

R = uconvert.(u"W/m^2", ϵ * σ .* uconvert.(u"K", T_ls) .^ 4)  # W/m^2
C = uconvert.(u"W/m^2", k1 * sqrt(V / D) .* (uconvert.(u"K", T_ls) .- uconvert(u"K", T_air)))  # W/m^2

# Noting that rtotal is used in the calculation of E
E = latent_energy_flux.(uconvert.(u"K", T_ls), h) # Defined in previous question 

# Total pathway flux (R + C + E)
RCE = R .+ C .+ E

# Now plot the pathways including total as before 
using Plots
p1 = plot(T_ls, RCE, 
    xlabel = "Leaf Temperature (°C)",
    ylabel = "Total Pathway Flux (W/m^2)",
    title = "Total Pathway Flux vs Leaf Temperature\n Tair = $T_air, V = $V, RH = $(RH * 100)%, Qa = $Qa",
    legend = :bottomright,
    label = "Total Pathway Flux"
)

plot!(p1, T_ls, R, label = "Radiative Flux")
plot!(p1, T_ls, C, label = "Convective Flux")
plot!(p1, T_ls, E, label = "Evaporative Heat Flux")

display(p1)


# Now determine the leaf temperature at Qa = 700 W/m^2
target_Qa = 700.0u"W/m^2"
i = argmin(abs.(RCE .- target_Qa))
leaf_temp = T_ls[i]
println("Leaf temperature at target when radiation absorbed is Qa = 700 W/m^2: $leaf_temp")


# Now determine the rate of water loss at this leaf temperature
E_i = E[i] / latent_heat(T_air)
water_loss = uconvert(u"kg/m^2/s", E_i)

println("Rate of water loss at leaf temperature $leaf_temp: $water_loss")







# 9. In the first example in the chapter, a plant leaf had the following
# properties: P_MLT = 0.05 mmole m-2 s-1, R = 200 s m-1, K = 10 mmole
# m-3 , and K_L = 100 W m-2. The environmental conditions were L = 400
# W m-2 , C_a = 12.5 mmole m-3 , and temperature such that G(T) = 1.0.
# Calculation showed that P = 0.0185 mmole m-2 S-l.

# Start by defining the specified variables 
P_MLT = 0.05u"mmol/m^2/s"
R = 200u"s/m"
K = 10u"mmol/m^3"
K_L = 100u"W/m^2"
L = 400u"W/m^2"
C_a = 12.5u"mmol/m^3"
G_T = 1.0


# What is the photosynthetic rate if K = 5 mmole m-3 ? (All other parameters
# and conditions are as specified initially.)

K = 5u"mmol/m^3"

# use equation 3.13 here
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 

# Then use Pm in the calculation for P
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at K = 5 mmole m-3: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# close enough 

# If K = 20 mmole m-3?
K = 20u"mmol/m^3"
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at K = 20 mmole m-3: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")



# If PMLT = 0.025 mmole m-2 S-l and all other parameters and conditions
# are as specified initially?
P_MLT = 0.025u"mmol/m^2/s"
R = 200u"s/m"
K = 10u"mmol/m^3"
K_L = 100u"W/m^2"
L = 400u"W/m^2"
C_a = 12.5u"mmol/m^3"
G_T = 1.0

Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at P_MLT = 0.025 mmole m-2 s-1: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")


# If PMLT = 0.10 mmole m-2 S-l?
P_MLT = 0.10u"mmol/m^2/s"

Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at P_MLT = 0.10 mmole m-2 s-1: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")



# If R = 400 s m-1 and all other parameters and conditions are as specified
# initially? R bring the resistance to CO2 diffusion.
P_MLT = 0.05u"mmol/m^2/s"
R = 400u"s/m"
K = 10u"mmol/m^3"
K_L = 100u"W/m^2"
L = 400u"W/m^2"
C_a = 12.5u"mmol/m^3"
G_T = 1.0

Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at R = 400 s m-1: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")



# 10. From the calculations made in the chapter and those in the above
# problem, rank the following plant parameters with regard to their effect
# upon photosynthesis, from the most to least important: R, K, and PMLT •
#
# Effect on photosynthesis:
# PMLT: −45% to +66%
# K:    approximately ±30%
# R:    −17% to +9%
#
# Therefore, the order of importance is:
# PMLT > K > R




# 11. Test the photosynthetic sensitivity of the plant to changes in environmental
# conditions by letting L = 200 and then 800 W m-2 , with all
# other conditions constant. Do the same for Ca = 6.25 and 25 mmole m-3 •
# Test the sensitivity to G(T) by letting it equal 0.5. From these calculations,
# rank the response of the plant to light, carbon dioxide concentration,
# and temperature in order of greatest to least sensitivity.

P_MLT = 0.05u"mmol/m^2/s"
R = 200u"s/m"
K = 10u"mmol/m^3"
K_L = 100u"W/m^2"
L = 400u"W/m^2"
C_a = 12.5u"mmol/m^3"
G_T = 1.0

Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at L = 200 W/m^2: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at L = 400 W/m^2: 0.0187 mmol m^-2 s^-1

L = 200u"W/m^2"
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at L = 200 W/m^2: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at L = 200 W/m^2: 0.0161 mmol m^-2 s^-1

L = 800u"W/m^2"
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at L = 800 W/m^2: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at L = 800 W/m^2: 0.0203 mmol m^-2 s^-1


# Do the same for Ca = 6.25 and 25 mmole m-3
P_MLT = 0.05u"mmol/m^2/s"
R = 200u"s/m"
K = 10u"mmol/m^3"
K_L = 100u"W/m^2"
L = 400u"W/m^2"
C_a = 12.5u"mmol/m^3"
G_T = 1.0

Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at C_a = 12.5 mmol/m^3: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at C_a = 12.5 mmol/m^3: 0.0187 mmol m^-2 s^-1

C_a = 6.25u"mmol/m^3"
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at C_a = 6.25 mmol/m^3: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at C_a = 6.25 mmol/m^3: 0.0114 mmol m^-2 s^-1

C_a = 25u"mmol/m^3"
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at C_a = 25 mmol/m^3: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at C_a = 25 mmol/m^3: 0.0265 mmol m^-2 s^-1


# Now for GT
# Do the same for Ca = 6.25 and 25 mmole m-3
P_MLT = 0.05u"mmol/m^2/s"
R = 200u"s/m"
K = 10u"mmol/m^3"
K_L = 100u"W/m^2"
L = 400u"W/m^2"
C_a = 12.5u"mmol/m^3"
G_T = 1.0

Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at G_T = 1.0: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at G_T = 1.0: 0.0187 mmol m^-2 s^-1


G_T = 0.5
Pm = (P_MLT * G_T) / (1 + (K_L / L)) 
P = ((C_a + K + R*Pm) - sqrt((C_a + K + R*Pm)^2 - 4*C_a*R*Pm)) / (2*R)
println("Photosynthetic rate at G_T = 0.5: $(round(typeof(0.0001u"mmol/m^2/s"), P, digits=4))")
# Photosynthetic rate at G_T = 0.5: 0.0102 mmol m^-2 s^-1


# Performance from varying the different environmental parameters:
# L (Light intensity)
# Photosynthetic rate at L = 200 W/m^2: 0.0161 mmol m^-2 s^-1 (reduction = 13.9%)
# Photosynthetic rate at L = 400 W/m^2: 0.0187 mmol m^-2 s^-1    
# Photosynthetic rate at L = 800 W/m^2: 0.0203 mmol m^-2 s^-1 (increase = 8.6%)

# C_a (Ambient CO2 concentration)
# Photosynthetic rate at C_a = 6.25 mmol/m^3: 0.0114 mmol m^-2 s^-1 (reduction = 39.0%)
# Photosynthetic rate at C_a = 12.5 mmol/m^3: 0.0187 mmol m^-2 s^-1
# Photosynthetic rate at C_a = 25 mmol/m^3: 0.0265 mmol m^-2 s^-1 (increase = 41.7%)

# G_T (Temperature factor)
# Photosynthetic rate at G_T = 1.0: 0.0187 mmol m^-2 s^-1
# Photosynthetic rate at G_T = 0.5: 0.0102 mmol m^-2 s^-1 (reduction = 45.5%)

# So the influence goes from G_T (Temperature factor) > C_a (Ambient CO2 concentration) > L (Light intensity)















