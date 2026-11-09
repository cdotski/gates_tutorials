



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
plot(
    ustrip.(thermal_max_Qa[1,:]), 
    ustrip.(uconvert.(u"°C", Tas)),
    label = "Thermal Max V=0.1 m/s", 
    xlabel= "Qa (W/m^2)", ylabel= "Air Temperature (C)", 
    legend = :bottomright, 
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






