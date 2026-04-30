


#=
==================================================
APPLICATION TO PLANTS
==================================================
=#


# Q1 | Plotting the blackbody radiation law
using Plots
using GLMakie
using CairoMakie
using Unitful
using Printf

# 0-60C in steps of 5C
temps = 0:5:60 #degrees C
sigma = 5.6697e-8 #W/m^2/K^4
f = Figure()
ax = Axis(f[1, 1], title = "Degrees C VS Radiation :)", 
    xlabel = "Celsius", 
    ylabel = "Radiation (W/m²)")
R = sigma * (temps .+ 273.15).^4 
lines!(ax, temps, R )    
f



# Q2 | Convection and leaf dimensions 
k1 = 0.664 # for smooth plates with Re < 5*10^5
V = [0.1, 1, 5] # m/s
D = [0.01, 0.05, 0.1, 0.3] # m
h_c = k1 * (sqrt.(V)' ./ sqrt.(D))
h_c

# Create a table with V as column headers and D as row labels
using PrettyTables
pretty_table(h_c, 
    column_labels = ["V = 0.1 m/s", "V = 1 m/s", "V = 5 m/s"],
    row_labels = ["D = 0.01 m", "D = 0.05 m", "D = 0.1 m", "D = 0.3 m"],
    title = "Convection Coefficients (W/m²K)")
#



# ==============================================
# Q3 Leaf temperature and evaporative water loss (at 25C)
# ==============================================

d_Te = 23   # Saturation density of water vapour in intercellular spaces at a given temp 
d_Ta = 23 # Saturation density of water vapour in the air at a given temp 
h = 0.5   # Relative humidity of the air (0-1)
re = 100  # Internal leaf resistance to water vapour diffusion (s/m)
ra = 50  # Boundary layer resistance to water vapour diffusion (s/m)

water_vap_dens_air = h * d_Ta
water_vap_dens_leaf = d_Te

E = (water_vap_dens_leaf - water_vap_dens_air) / (re + ra)
lambda = 2.45e6 # J/kg, latent heat of vaporization of water
E_energy = E * lambda # W/m², energy lost due to evaporation


# ==================================================
# Q3.1 Leaf temperature and energy lost from evaporation at 30C
# ==================================================
using Printf
using Unitful

# Constants 
Rv = 461.5u"J/kg/K"          # gas constant for water vapour
λ0 = 2.501e6u"J/kg"          # latent heat at 0°C
λ_slope = 2.37e3u"J/kg/K"    # decline per °C

Ta = 30u"°C"                # air temperature
r_total = 200u"s/m"         # r_l + r_a
RHs = [0.1, 0.5, 0.9]       # relative humidities

# Saturation vapour pressure Tetens equation
# T in °C, output in Pa
function saturation_vapour_pressure(T)
    Tc = ustrip(u"°C", T)
    es_kPa = 0.6108 * exp((17.27 * Tc) / (Tc + 237.3))
    return es_kPa * 1000u"Pa"
end

# Saturation vapour density using ideal gas law
function saturation_vapour_density(T)
    es = saturation_vapour_pressure(T)
    Tk = uconvert(u"K", T)
    return uconvert(u"kg/m^3", es / (Rv * Tk))
end

# Latent heat of vaporization at temperature T
function latent_heat(T)
    Tc = ustrip(u"°C", T)
    return (2.501e6 - 2.37e3 * Tc)u"J/kg"
end

# Evaoporation rate
function evaporation_rate(Tl, Ta, h, r_total)
    ρ_leaf = saturation_vapour_density(Tl)
    ρ_air_sat = saturation_vapour_density(Ta)
    return uconvert(u"kg/m^2/s", (ρ_leaf - h * ρ_air_sat) / r_total)
end

# -----------------------------
# Latent energy flux
# λE in W/m^2
# -----------------------------
function latent_flux(Tl, Ta, h, r_total)
    E = evaporation_rate(Tl, Ta, h, r_total)
    λ = latent_heat(Tl)
    return uconvert(u"W/m^2", λ * E)
end


# ===============================================
# Generate data
# ===============================================
using Plots

leaf_temps = collect(20:0.5:45) .* u"°C"

plot()

for h in RHs
    λE = [ustrip(u"W/m^2", latent_flux(Tl, Ta, h, r_total)) for Tl in leaf_temps]
    Tl_vals = ustrip.(u"°C", leaf_temps)

    plot!(
        λE,
        Tl_vals,
        label = "RH = $(Int(h*100))%",
        xlabel = "Latent heat flux, λE (W m⁻²)",
        ylabel = "Leaf temperature (°C)",
        linewidth = 2
    )
end

plot!(
    title = "Evaporative energy loss from a leaf",
    legend = :bottomright
)


# Converting this plot to total Watts convection for a leaf of 4x5cm 
using Unitful

leaf_temps = collect(20:0.5:45) .* u"°C"

V = 5u"m/s"
D = 0.05u"m"
leaf_area = uconvert(u"m^2", 5u"cm" * 4u"cm")

k1 = 9.14u"J/m^2/s^0.5/K"

h_c = k1 * sqrt(V / D)

convection_flux = h_c .* (leaf_temps .- 30u"°C")   # W/m^2
convection_power = uconvert.(u"W", convection_flux .* leaf_area)

# Calculate evaporative power for the same leaf to plot against convection 
evaporation_flux = [ustrip(u"W/m^2", latent_flux(Tl, Ta, 0.5, r_total)) for Tl in leaf_temps]
evaporation_power = evaporation_flux .* ustrip(u"m^2", leaf_area)

# Plot both on same figure with temperature on y-axis
Tl_vals = ustrip.(u"°C", leaf_temps)

plot(
    ustrip.(u"W", convection_power), Tl_vals,
    label = "Convection",
    xlabel = "Power (W)",
    ylabel = "Leaf temperature (°C)",
    linewidth = 2
)

plot!(
    evaporation_power, Tl_vals,
    label = "Evaporation (50% RH)",
    linewidth = 2,
    title = "Heat transfer of 5×4 cm leaf at 30°C Ta and 1 m/s wind",
    legend = :bottomright
)









# ================================================
# Q4 Te versus R + C + E
# ================================================

D = 0.05u"m"
rera = 200u"s/m"
V = 0.1u"m/s"
T_a = 30u"°C"
h = 0.5 # relative humidity

# Calculate reradiated energy
sigma = 5.6697e-8u"W/m^2/K^4"
R = uconvert.(u"W/m^2", sigma .* uconvert.(u"K", T_ls) .^ 4)  # W/m^2

#Leaf temps 
T_ls = collect(20:0.5:45) .* u"°C"

# Calculate convection energy
k1 = 9.14u"J/m^2/s^0.5/K"
h_c = k1 * sqrt(V / D)
C = uconvert.(u"W/m^2", h_c .* (uconvert.(u"K", T_ls) .- uconvert(u"K", T_a)))  # W/m^2

# Calculate evaporative energy
E = [uconvert.(u"W/m^2", latent_flux(Tl, T_a, h, rera)) for Tl in T_ls] # W/m^2
# Total energy loss
total_loss = R .+ C .+ E

total_loss_4x5cm = total_loss .* ustrip(u"m^2", leaf_area) # W

# Plot T_ls vs total loss 
params_str = @sprintf("D=%.2fm, r=%.0fs/m, V=%.1fm/s, Ta=%.0f°C, RH=%.0f%%",
    ustrip(u"m", D), ustrip(u"s/m", rera), ustrip(u"m/s", V), ustrip(u"°C", T_a), h * 100)

plot(
    ustrip.(u"W/m^2", total_loss), Tl_vals,
    label = "Total energy loss",
    xlabel = "Energy loss C + R + E (W/m²)",
    ylabel = "Leaf temperature (°C)",
    linewidth = 2,
    title = "Leaf temperature vs total energy loss\n" * params_str,
    legend = :bottomright
)


# Plot T_ls vs total loss for a 4x5cm leaf
plot(
    ustrip.(u"W", total_loss_4x5cm), Tl_vals,
    label = "Total energy loss (4×5 cm leaf)",
    xlabel = "Energy loss (W/m²)",
    ylabel = "Leaf temperature (°C)",
    linewidth = 2,
    title = "Leaf temperature vs total energy loss \n" * params_str,
    legend = :bottomright
)






# ================================================
# Q5 Equilibrium leaf temperature for 500 W/m² absorbed
# ================================================

D = 0.05u"m"
rera = 200u"s/m"
V = 0.1u"m/s"
T_a = 30u"°C"
h = 0.5 # relative humidity

#Pkg.add("Roots")
using Roots

absorbed = 500.0u"W/m^2"  # W/m²

" ## Find leaf Temp
Energy balance function to find equilibrium leaf temperature for a given absorbed radiation.
All the energy loss minus absorbed should equal zero at equilibrium. The function takes leaf
 temperature as input and calculates the total energy loss (radiation + convection + evaporation) and subtracts the absorbed radiation. We then use a root-finding algorithm to find the temperature at which this balance is zero.    
    "
function energy_balance(Tc)
    Tl = Tc * u"°C"
    R_l = uconvert(u"W/m^2", sigma * uconvert(u"K", Tl)^4)
    C_l = uconvert(u"W/m^2", h_c * (uconvert(u"K", Tl) - uconvert(u"K", T_a)))
    E_l = uconvert(u"W/m^2", latent_flux(Tl, T_a, h, rera))
    return R_l + C_l + E_l - absorbed
end

T_eq = (find_zero(energy_balance, (20.0, 60.0)))u"°C" # Finds the root bewteen 20 and 60C, i.e. the equilibrium leaf temperature where energy balance is zero
@printf("Equilibrium leaf temperature absorbing %.0f W/m² is: %.2f °C\n", ustrip(u"W/m^2", absorbed), ustrip(u"°C", T_eq))

# Proportions lost to each process at equilibrium
function energy_losses(T_eq, absorbed)
    Tl = T_eq
    R_prop = uconvert(u"W/m^2", sigma * uconvert(u"K", Tl)^4)
    C_prop = uconvert(u"W/m^2", h_c * (uconvert(u"K", Tl) - uconvert(u"K", T_a)))
    E_prop = uconvert(u"W/m^2", latent_flux(Tl, T_a, h, rera))
    return ustrip(R_prop / absorbed),
           ustrip(C_prop / absorbed),
           ustrip(E_prop / absorbed)
end

R_loss, C_loss, E_loss = energy_losses(T_eq, absorbed)
@printf("At equilibrium, proportion of energy lost to radiation: %.2f%%\n", R_loss * 100)
@printf("At equilibrium, proportion of energy lost to convection: %.2f%%\n", C_loss * 100)
@printf("At equilibrium, proportion of energy lost to evaporation: %.2f%%\n", E_loss * 100)

#= The -10% is because the air is hotter than the leaf, so convection adds energy 
but condsider that if you add the proportions add up to 0 so the temperature is stable 
at the equilibrium temperature. 
=#

R_loss + C_loss + E_loss  # should be close to 1, i.e. all 
                          # absorbed energy is lost through these processes at equilibrium


# How much water is lost at equilibrium temp in kg/m^2/s?
water_loss = uconvert(u"kg/m^2/s", evaporation_rate(T_eq, T_a, h, rera))
water_loss_g_per_m2_per_day = uconvert(u"g/m^2/d", water_loss)
@printf("At equilibrium, water loss is %.2f\n", water_loss_g_per_m2_per_day)







# =======================================================
# Q6 For h = 6 determine leaf temp 
# =======================================================
absorbed = 700u"W/m^2"
h = 0.1

# Solve for leaf temperature at this absorbed radiation and humidity
T_eq2 = (find_zero(energy_balance, (20.0, 60.0)))u"°C"

# Proportions lost to each process at equilibrium
R_loss2, C_loss2, E_loss2 = energy_losses(T_eq2, absorbed)
@printf("At equilibrium, proportion of energy lost to radiation: %.2f%%\n", R_loss2 * 100)
@printf("At equilibrium, proportion of energy lost to convection: %.2f%%\n", C_loss2 * 100)
@printf("At equilibrium, proportion of energy lost to evaporation: %.2f%%\n", E_loss2 * 100)

R_loss2 + C_loss2 + E_loss2  # should be close to 1
# Reduced h = increased evaoprative loss, so leaf temp is lower at equilibrium and convection
# adds energy to the leaf rather than removing it, so proportion lost to convection is negative.







# ==========================================================
# Q7 Table of boundary layer resistance (ra) for m(resistance measure) and V/D
# ==========================================================

# Q2 | Convection and leaf dimensions 

m = [0.1, 1, 10]u"s/m" # resistance measure
D_V = [1, 10, 100] # V/D ratio
k2 = 200u"s^0.5/m" # proportionality constant for ra
ra = k2 .* (m').^0.2 .* sqrt.(D_V) # 3×3 matrix: rows = V/D, cols = m

# Create a table with V as column headers and D as row labels
using PrettyTables
ra_display = round.(ustrip.(ra), digits=2)
pretty_table(ra_display, 
    column_labels = ["m = 0.1 s/m", "m = 1 s/m", "m = 10 s/m"],
    row_labels = ["D/V = 1", "D/V = 10", "D/V = 100"],
    title = "Boundary Layer Resistance (s/m)")
#



# ==========================================================
# Q8 Determining leaf temperature 
# ==========================================================

# Environmental conditions
D = 0.05u"m" # leaf dimension
W = 0.05u"m" # Width
V = 0.1u"m/s" # wind speed
T_a = 30u"°C" # air temperature
h = 0.5 # relative humidity 
re = 100u"s/m" # internal leaf resistance to water vapour diffusion
ra = 100u"s/m" # boundary layer resistance to water vapour diffusion
rera = re + ra # total resistance to water vapour diffusion
epsilon = 0.96 # leaf emissivity

# Energy absorbed by the leaf
absorbed = 700u"W/m^2"

# Restating the energy balance function
function energy_balance(Tc)
    Tl = Tc * u"°C"
    R_l = uconvert(u"W/m^2", epsilon * sigma * uconvert(u"K", Tl)^4)
    C_l = uconvert(u"W/m^2", h_c * (uconvert(u"K", Tl) - uconvert(u"K", T_a)))
    E_l = uconvert(u"W/m^2", latent_flux(Tl, T_a, h, rera))
    return R_l + C_l + E_l - absorbed
end

T_eq = (find_zero(energy_balance, (20.0, 60.0)))u"°C" # Finds the root bewteen 20 and 60C, i.e. the equilibrium leaf temperature where energy balance is zero
@printf("Equilibrium leaf temperature absorbing %.0f W/m² is: %.2f °C\n", ustrip(u"W/m^2", absorbed), ustrip(u"°C", T_eq))




# ==========================================================
# Q9 Photosynthetic rate P for given conditions
# ==========================================================

# Initial conditions
P_mlt = 0.05u"mmol/m^2/s" # molar transpiration rate
R = 200u"s/m" # total resistance to water vapour diffusion
K = 10u"mmol/m^3" # K is a constant equal to the chloroplast concentration of CO2 at which P = PM/2.
K_l = 100u"W/m^2" # light intensity at which photosynthesis is half of its maximum
L = 400u"W/m^2" # light intensity between 400 and 700 W/m^2, where photosynthesis is saturated
C_a = 12.5u"mmol/m^3" # ambient CO2 concentration in air 
G_T = 1 # efficiency of photosynthetic rate at temp T

# P = 0.0185u"mmole/m^2/s" # photosynthetic rate in this example 

# Get P_M, the maximum photosynthetic rate at light and carbon dioxide saturation and optimum temperature, using the equation:
P_M = (P_mlt * G_T) / (1 + (K_l/L)) # the value of P_M at light and carbon dioxide saturation and optimum temperature.

# Calculate the actual photosynthetic rate P
P = ((C_a + K + R * P_M) - sqrt((C_a + K + R * P_M)^2 - 4 * C_a * R * P_M)) / (2 * R)
@printf("Photosynthetic rate P is: %.4f mmol/m^2/s\n", ustrip(u"mmol/m^2/s", P))

# If K = 5mmol/m^3, how does this affect P?
K = 5u"mmol/m^3"
P = ((C_a + K + R * P_M) - sqrt((C_a + K + R * P_M)^2 - 4 * C_a * R * P_M)) / (2 * R)
@printf("Photosynthetic rate P is: %.4f mmol/m^2/s\n", ustrip(u"mmol/m^2/s", P))

#If PMLT = 0.025 mmole m-2 S-l and all other parameters and conditions are as specified initially?
P_mlt = 0.025u"mmol/m^2/s"
P_M = (P_mlt * G_T) / (1 + (K_l/L))
P = ((C_a + K + R * P_M) - sqrt((C_a + K + R * P_M)^2 - 4 * C_a * R * P_M)) / (2 * R)
@printf("Photosynthetic rate P is: %.4f mmol/m^2/s\n", ustrip(u"mmol/m^2/s", P))

# If PMLT = 0.10 mmole m-2 S-l?
P_mlt = 0.10u"mmol/m^2/s"
P_M = (P_mlt * G_T) / (1 + (K_l/L))
P = ((C_a + K + R * P_M) - sqrt((C_a + K + R * P_M)^2 - 4 * C_a * R * P_M)) / (2 * R)
@printf("Photosynthetic rate P is: %.4f mmol/m^2/s\n", ustrip(u"mmol/m^2/s", P))

# If R = 400 s m-1 and all other parameters and conditions are as specified initially?
R = 400u"s/m"
P_M = (P_mlt * G_T) / (1 + (K_l/L))
P = ((C_a + K + R * P_M) - sqrt((C_a + K + R * P_M)^2 - 4 * C_a * R * P_M)) / (2 * R)
@printf("Photosynthetic rate P is: %.4f mmol/m^2/s\n", ustrip(u"mmol/m^2/s", P))



# ===========================================================
# Q10 Sensitivity of P to changes in R, K and PMLT
# ===========================================================

# Making a function for P_M
# Function to calculate photosynthetic rate with configurable parameters
function calculate_photosynthetic_rate(;
    P_mlt = 0.05u"mmol/m^2/s",
    G_T = 1,
    K_l = 100u"W/m^2",
    L = 400u"W/m^2",
    C_a = 12.5u"mmol/m^3",
    K = 10u"mmol/m^3",
    R = 200u"s/m"
)
    P_M = (P_mlt * G_T) / (1 + (K_l/L))
    P = ((C_a + K + R * P_M) - sqrt(((C_a + K + R * P_M)^2) - (4 * C_a * R * P_M))) / (2 * R)
    return P
end

# Example usage of the function with default parameters
P = calculate_photosynthetic_rate()
P_r = calculate_photosynthetic_rate(R = 400u"s/m")
P_l = calculate_photosynthetic_rate(K = 5u"mmol/m^3")
P_mlt_example = calculate_photosynthetic_rate(P_mlt = 0.10u"mmol/m^2/s")

#=
Doubling R reduces P by 0.003 mmol/m^2/s, which is a 6% decrease from the original P.
Halving K increased P by 0.006 mmol/m^2/s, which is a 12% increase from the original P.
Doubling P_mlt increased P by 0.012 mmol/m^2/s, which is a 24% increase from the original P.

It appears that P is most sensitive to changes in P_mlt, followed by K, and 
least sensitive to changes in R.
=#



# =============================================================
# Q11 Sensitivity of P to changes in L
# =============================================================

# Varying L from 200 to 800 W/m^2 and calculating P for each value
P = calculate_photosynthetic_rate()
P_l200 = calculate_photosynthetic_rate(L = 200u"W/m^2")
P_l800 = calculate_photosynthetic_rate(L = 800u"W/m^2")

# Varying C_a from 6.25 to 25 mmol/m^3 and calculating P for each value
P = calculate_photosynthetic_rate()
P_ca6 = calculate_photosynthetic_rate(C_a = 6.25u"mmol/m^3")
P_ca25 = calculate_photosynthetic_rate(C_a = 25u"mmol/m^3")

# Varying G_T from 0.5 to 1 and calculating P for each value
P = calculate_photosynthetic_rate()
P_gHalf = calculate_photosynthetic_rate(G_T = 0.5)

#=
Doubling L from 400 to 800 W/m^2 increased P by 0.002 mmol/m^2/s, which is a 4% increase from the original P.
Doubling Ca from 12.5 to 25 mmol/m^3 increased P by 0.008 mmol/m^2/s which is a 16% increase from the original P.
Halving G_T from 1 to 0.5 decreased P by 0.008 mmol/m^2/s which is a 16% decrease from the original P.

P appears to be more sensitive to changes in G_T and C_a than to changes in L, with P being most sensitive to changes in G_T.
=# 

 



