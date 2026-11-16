



# 1. Determine the climate-space for a desert iguana, Dipsosaurus dorsalis,
# with the characteristics E = 1.0 and body diameter D = 0.015 m. At
# thermal maximum, M = 10.5 W m-2 , I= 0.005 m2°C W-l, λE = 6.3 W
# m-2, and Tb = 45°C. At thermal minimum, M = 0.1 W m-2, I=0.005
# m2°C W-1, λE = 0.1 W m-2, and Tb = 3°C.
# For air flow speeds of V = 0.1 m S-l and 2.0 m S-l, determine the
# climate-space limits at thermal maximum and thermal minimum in the
# form Qa = A + B*Ta.    Plot these lines on a graph with the same scales as Figs. 4.2 and 4.3.
# Now place on the graph the radiation limits as shown in Fig. 4.3. These
# lines now form an approximate climate-space for the desert iguana.
using Unitful, Plots


D = 0.015u"m" 
M_max = 10.5u"W/m^2"
I_m = 0.005u"m^2*K/W"
λE_max = 6.3u"W/m^2"
Tb_max = 45.0u"°C"

M_min = 0.1u"W/m^2"
I_min = 0.005u"m^2*K/W"
λE_min = 0.1u"W/m^2"
Tb_min = 3.0u"°C"


# State the basic parameters 
ϵ = 1 # Emissivity 
σ = 5.67e-8u"W/m^2/K^4" # Stefan Boltzmann constant 

#At 30°C, 
λ = 2.43e6u"J/kg" # Latent heat of vaporization of water at 30°C in the book
k1 = 3.89u"W*s^0.5/m^2/K" 

# Use equation 4.6 for this question
# Now set some temperatures for the climate-space limits
Tas = [-60.0:0.5:60.0;]u"°C" # Air temperature range from -60 to 60°C
Tas = uconvert.(u"K", Tas) # Convert to Kelvin for calculations
Vs = [0.1, 2.0]u"m/s" # Air flow speeds of 0.1 m/s and 2.0 m/s


# To balance we must have M - λE + Qa = radiation + convection (and + conduction but that comes later)
# So it should take the form of Qa = R + C + λE - M
radiation_max = ϵ * σ * (uconvert(u"K", Tb_max) - I_m*(M_max - λE_max))^4
convection_max = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_max .- reshape(Tas, 1, :) .- I_m*(M_max - λE_max))

thermal_max_Qa = radiation_max .+ convection_max .+ λE_max .- M_max

# Now do the same for thermal minimum
radiation_min = ϵ * σ * (uconvert(u"K", Tb_min) - I_m*(M_min - λE_min))^4
convection_min = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_min .- reshape(Tas, 1, :) .- I_m*(M_min - λE_min))

thermal_min_Qa = radiation_min .+ convection_min .+ λE_min .- M_min



# Now draw the climate-space limits 
# this is the same as applying the radiation law I just converted to Kelvin inside 
# instead of adding 273.15 to C manually.
R_bb = σ * (Tas).^4  # Blackbody radiation law (W/m^2)






# Now plot the results
d_iguana = plot(
    ustrip.(thermal_max_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=0.1 m/s", 
    xlabel= "Qa (W/m^2)", ylabel= "Air Temperature (C)", 
    legend = :outerright, 
    ylims = (-60, 60),
    xlims = (100, 1000),
    title = "Climate-Space for Desert Iguana", 
    aspect_ratio = 10)
plot!(
    ustrip.(thermal_max_Qa[2,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=2.0 m/s"
)
plot!(
    ustrip.(thermal_min_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Min V=0.1 m/s"
)
plot!(
    ustrip.(thermal_min_Qa[2,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Min V=2.0 m/s"
)


# Now add the radiation limits as shown in Fig. 4.3

plot!(
    ustrip.(R_bb), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Blackbody Radiation")

# Note that using these parameters specified in the question, 
# M - λE at thermal maximum is actually = 4.2 W/m^2,
# so we don't see the lines intersect at the BB radiation line 
# Idk why he did it this way



# Now he doesn't give you the values for the sunny day and night sky
# so we're just gonna approximate them here for funsies 
# Just a very rough approximation at this point 

ε_sky = 0.6                  # Assumed effective clear-sky emissivity
α_solar = 0.8                # Solar absorptance used in the text

R_night_clear = 0.5 .* (1 + ε_sky) .* R_bb

S_eff = range(
    ustrip((200u"W/m^2" - R_night_clear[1]) / α_solar),
    ustrip((1000u"W/m^2" - R_night_clear[end]) / α_solar);
    length=length(Tas)
) .* u"W/m^2"

R_day_clear = R_night_clear .+ α_solar .* S_eff


plot!(
    ustrip.(R_night_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Night, clear sky",
)

plot!(
    ustrip.(R_day_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Day, sun, clear sky",
    linestyle=:dash,
)


display(d_iguana)




# Question 2
#  2. Determine the climate-space for a masked shrew (Sorex cinereus) with ε = 1.0 and D = 0.018 m.
# Thermal max: M = 140 W m^-2, I = 0.085 m^2 °C W^-1, λE = 28 W m^-2, T_b = 41 °C.
# Thermal min: M = 349 W m^-2, I = 0.125 m^2 °C W^-1, λE = 3.5 W m^-2, T_b = 37.5 °C.
# For air flow V = 0.1 m s^-1, find limits as: Q_a = A + B * T_a
# Plot these lines and radiation limits (Fig. 4.3) to form the shrew's climate-space (cross-hatch habitat).

D = 0.018u"m" 
M_max = 140.0u"W/m^2"
I_max = 0.085u"m^2*K/W"
λE_max = 28.0u"W/m^2"
Tb_max = 41.0u"°C"

M_min = 349.0u"W/m^2"
I_min = 0.125u"m^2*K/W"
λE_min = 3.5u"W/m^2"
Tb_min = 37.5u"°C"

# State the basic parameters 
ϵ = 1 # Emissivity 
σ = 5.67e-8u"W/m^2/K^4" # Stefan Boltzmann constant 

#At 30°C
λ = 2.43e6u"J/kg" # Latent heat of vaporization of water at 30°C in the book
k1 = 3.89u"W*s^0.5/m^2/K" 

# Now set some temperatures for the climate-space limits
Tas = [-60.0:0.5:60.0;]u"°C" # Air temperature range from -60 to 60°C
Tas = uconvert.(u"K", Tas) # Convert to Kelvin for calculations
Vs = [0.1]u"m/s" # Just for still air 


# To balance we must have M - λE + Qa = radiation + convection (and + conduction but that comes later)
# So it should take the form of Qa = R + C + λE - M
radiation_max = ϵ * σ * (uconvert(u"K", Tb_max) - I_max*(M_max - λE_max))^4
convection_max = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_max .- reshape(Tas, 1, :) .- I_max*(M_max - λE_max))

thermal_max_Qa = radiation_max .+ convection_max .+ λE_max .- M_max

# Now do the same for thermal minimum
radiation_min = ϵ * σ * (uconvert(u"K", Tb_min) - I_min*(M_min - λE_min))^4
convection_min = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_min .- reshape(Tas, 1, :) .- I_min*(M_min - λE_min))

thermal_min_Qa = radiation_min .+ convection_min .+ λE_min .- M_min




# Now plot the results
m_shrew = plot(
    ustrip.(thermal_max_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=0.1 m/s", 
    xlabel= "Qa (W/m^2)", ylabel= "Air Temperature (C)", 
    legend = :outerright, 
    ylims = (-60, 60),
    xlims = (100, 1000),
    title = "Climate-Space for Masked Shrew", 
    aspect_ratio = 10)
plot!(
    ustrip.(thermal_min_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Min V=0.1 m/s"
)

# Add BB radiation line 
plot!(
    ustrip.(R_bb), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Blackbody Radiation")

# Add night and sunny day radiation lines
plot!(
    ustrip.(R_night_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Night, clear sky",
)

plot!(
    ustrip.(R_day_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Day, sun, clear sky",
    linestyle=:dash,
)

display(m_shrew)


# 3. If the shrew with the same characteristics described in Problem 4.2
# is without any convective heat exchange, what will its climate-space
# limits be? What do these calculations suggest to you in terms of the habitat
# in which the shrew actually lives?

Vs = [0.0]u"m/s" # No convection 

# So it should take the form of Qa = R + C + λE - M
radiation_max = ϵ * σ * (uconvert(u"K", Tb_max) - I_max*(M_max - λE_max))^4
convection_max = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_max .- reshape(Tas, 1, :) .- I_max*(M_max - λE_max))

thermal_max_Qa = radiation_max .+ convection_max .+ λE_max .- M_max

# Now do the same for thermal minimum
radiation_min = ϵ * σ * (uconvert(u"K", Tb_min) - I_min*(M_min - λE_min))^4
convection_min = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_min .- reshape(Tas, 1, :) .- I_min*(M_min - λE_min))

thermal_min_Qa = radiation_min .+ convection_min .+ λE_min .- M_min



# Now plot the results
m_shrew_no_conv = plot(
    ustrip.(thermal_max_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=0 m/s", 
    xlabel= "Qa (W/m^2)", ylabel= "Air Temperature (C)", 
    legend = :outerright, 
    ylims = (-60, 60),
    xlims = (-200, 1000),
    title = "Climate-Space for Masked Shrew", 
    aspect_ratio = 10)
plot!(
    ustrip.(thermal_min_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Min V=0 m/s"
)

# Add BB radiation line 
plot!(
    ustrip.(R_bb), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Blackbody Radiation")

# Add night and sunny day radiation lines
plot!(
    ustrip.(R_night_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Night, clear sky",
)

plot!(
    ustrip.(R_day_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Day, sun, clear sky",
    linestyle=:dash,
)

display(m_shrew_no_conv)

# Without convection, the shrew's climate space is tiny 
# the climate space therefore suggests the shrew would live in 
# a habitat with some air movement such as a forest or grassland 



# 4. You are holding a shrew in the palm of your hand. Consider your
# hand to be a blackbody cavity at 32°C. Will the shrew survive a steadystate
# condition in this situation? (Use the climate-space diagram determined
# in Problem 4.2.)

plot(m_shrew, xlims=(100, 1000), ylims=(-60, 60))
plot!(
    [0.0, 1000.0], 
    [32.0, 32.0], 
    label = "Hand Temperature (32°C)", 
    color = :red,
    linestyle = :dashdot
)

# The shew will not survive as this line is above the thermal maximum 
# for black body radiation. RIP. 





# 5. Determine the climate-space for a domestic sheep, Ovis ovis, with 
# the characteristics ε = 1.0 and body diameter D = 0.25 m.
# Thermal max: M = 91 W m^-2, I = 5.1 m^2 °C W^-1, λE = 95 W m^-2, T_b = 41.7 °C.
# Thermal min: M = 70 W m^-2, I = 4.0 m^2 °C W^-1, λE = 13 W m^-2, T_b = 39.5 °C.

D = 0.25u"m" 
M_max = 91.0u"W/m^2"
I_max = 5.1u"m^2*K/W"
λE_max = 95.0u"W/m^2"
Tb_max = 41.7u"°C"

M_min = 70.0u"W/m^2"
I_min = 4.0u"m^2*K/W"
λE_min = 13.0u"W/m^2"
Tb_min = 39.5u"°C"

# State the basic parameters 
ϵ = 1 # Emissivity 
σ = 5.67e-8u"W/m^2/K^4" # Stefan Boltzmann constant 

#At 30°C
λ = 2.43e6u"J/kg" # Latent heat of vaporization of water at 30°C in the book
k1 = 3.89u"W*s^0.5/m^2/K" 

# Now set some temperatures for the climate-space limits
Tas = [-60.0:0.5:60.0;]u"°C" # Air temperature range from -60 to 60°C
Tas = uconvert.(u"K", Tas) # Convert to Kelvin for calculations
Vs = [0.1, 2.0]u"m/s" # For still air and 2.0 m/s. He doesn't specify but i like this comparison.
R_bb = σ * (Tas).^4  # Blackbody radiation law (W/m^2)
R_night_clear = 0.5 .* (1 + ε_sky) .* R_bb

S_eff_sheep = range(
    ustrip((200.0u"W/m^2" - R_night_clear[1]) / α_solar),
    ustrip((1000.0u"W/m^2" - R_night_clear[end]) / α_solar);
    length=length(Tas)
) .* u"W/m^2"

R_day_clear_sheep = R_night_clear .+ α_solar .* S_eff_sheep


# To balance we must have M - λE + Qa = radiation + convection (and + conduction but that comes later)
# So it should take the form of Qa = R + C + λE - M
radiation_max = ϵ * σ * (uconvert(u"K", Tb_max) - I_max*(M_max - λE_max))^4
convection_max = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_max .- reshape(Tas, 1, :) .- I_max*(M_max - λE_max))

thermal_max_Qa = radiation_max .+ convection_max .+ λE_max .- M_max

# Now do the same for thermal minimum
radiation_min = ϵ * σ * (uconvert(u"K", Tb_min) - I_min*(M_min - λE_min))^4
convection_min = k1 .* (sqrt.(reshape(Vs, :, 1)) ./ sqrt(D)) .*(Tb_min .- reshape(Tas, 1, :) .- I_min*(M_min - λE_min))

thermal_min_Qa = radiation_min .+ convection_min .+ λE_min .- M_min


# Now for plotting 
# currently I have no idea what's going on here



sheep = plot(
    ustrip.(thermal_max_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=0.1 m/s", 
    xlabel= "Qa (W/m^2)", ylabel= "Air Temperature (C)", 
    legend = :outerright, 
    ylims = (-60, 80),
    xlims = (0, 1000),
    title = "Climate-Space for a Domestic Sheep", 
    aspect_ratio = 10)
plot!(
    ustrip.(thermal_max_Qa[2,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=2.0 m/s"
)
plot!(
    ustrip.(thermal_min_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Min V=0.1 m/s"
)
plot!(
    ustrip.(thermal_min_Qa[2,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Min V=2.0 m/s"
)

# Add BB radiation line 
plot!(
    ustrip.(R_bb), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Blackbody Radiation")

# Add night and sunny day radiation lines
plot!(
    ustrip.(R_night_clear),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Night, clear sky",
)

plot!(
    ustrip.(R_day_clear_sheep),
    ustrip.(uconvert.(u"°C", Tas)),
    label="Day, sun, clear sky",
    linestyle=:dash,
)

display(sheep)












