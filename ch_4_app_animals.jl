# ===========================================================
# Q 0.1 plotting radiation on clear night, BB radiation, and clear day
# ===========================================================














# ===========================================================
# Q1 Determine the climate-space for a desert iguana, Dipsosaurus dorsalis
# ===========================================================
using Unitful

# Define the parameters for the desert iguana
epsilon = 1.0 # Emissivity
sigma = 5.67e-8u"W/m^2/K^4" # Stefan-Boltzmann constant
D = 0.015u"m" # Body diameter
I = 0.005u"m^2*K/W" # Insulation quality
M_max = 10.5u"W/m^2" # Metrate at Thermal maximum 
M_min = 0.1u"W/m^2" # Metrate at Thermal minimum
lambda_E_max = 6.3u"W/m^2" # Evaporative cooling rate
lambda_E_min = 0.1u"W/m^2" # Evaporative cooling rate
Tb_max = uconvert(u"K", 45.0u"°C") # Body temperature
Tb_min = uconvert(u"K", 3.0u"°C") # Body temperature

# Define the environment parameters
Vs = [0.1, 2]u"m/s" # Wind speed
Tas = uconvert.(u"K", [-60:1:60;]u"°C") # Air temperature range
#Ta = 25 # Air temperature

# Maximum temp
k1 = 3.89u"W/m^2*s^(1/2)/K"

# Qa at max body temp at 0.1m/s
# Just applying the @. makes it a FOR LOOP for every Ta in Tas without having to do it manually 
Qa_MAX_free = @. (epsilon * sigma * (Tb_max - I*(M_max - lambda_E_max))^4) +
            k1 * (Vs[1]^0.5 / D^0.5) * (Tb_max - Tas - I*(M_max - lambda_E_max)) - (M_max + lambda_E_max)
            
# Qa max with 2m/s
Qa_MAX_2 = @. (epsilon * sigma * (Tb_max - I*(M_max - lambda_E_max))^4) +
            k1 * (Vs[2]^0.5 / D^0.5) * (Tb_max - Tas - I*(M_max - lambda_E_max)) - (M_max + lambda_E_max)

# Qa min with 0.1m/s
Qa_MIN_free = @. (epsilon * sigma * (Tb_min - I*(M_min - lambda_E_min))^4) +
            k1 * (Vs[1]^0.5 / D^0.5) * (Tb_min - Tas - I*(M_min - lambda_E_min)) - (M_min + lambda_E_min)

# Qa min with 0.1m/s
Qa_MIN_2 = @. (epsilon * sigma * (Tb_min - I*(M_min - lambda_E_min))^4) +
            k1 * (Vs[2]^0.5 / D^0.5) * (Tb_min - Tas - I*(M_min - lambda_E_min)) - (M_min + lambda_E_min)



# Hand calculations for max and min Qa at V = 0.1 m/s and V = 2 m/s
#Qa_free_max = [1031.249 - 10.04 * Ta for Ta in Tas] #with V = 0.1 m/s
#Qa_2_max = [2599.98 - 44.91 * Ta for Ta in Tas] #with V = 2 m/s
#Qa_free_min = [508.04 - 10.04 * Ta for Ta in Tas] #with V = 0.1 m/s
#Qa_2_min = [642.77 - 44.91 * Ta for Ta in Tas] #with V = 2 m/s  



# Environmental radiation lines
alpha_s = 0.8   # Solar absorptance
alpha_bb = 1.0  # Blackbody absorptance
S_solar = 600.0u"W/m^2" # Direct + diffuse solar for clear day (W/m²)
S_solar_night = 0.0u"W/m^2" # No solar for clear night

# Clear night arbitrarily set  
Ra_clear = @. 1.22 * sigma * Tas^4 - 171u"W/m^2"
Rg_ground = @. sigma * Tas^4
Qa_night = @. 0.5 * Ra_clear + 0.5 * Rg_ground

# Clear day with sun: same longwave + absorbed solar
Qa_day = @. ((sigma * (Tas)^4) + alpha_s * S_solar)

# Blackbody: all surroundings at Ta (maximum possible longwave environment)
Qa_blackbody = @. (alpha_bb * sigma * (Tas)^4) 

using Plots


# With the quick calcs
Tas_C = ustrip.(u"°C", Tas)

plot(ustrip.(u"W/m^2", Qa_MAX_free), Tas_C,
    label="V = 0.1 m/s (Max)",
    xlims=(0, 1200),
    ylims=(-20, 70),
    xlabel="Radiation (W/m²)",
    ylabel="Air Temperature (°C)",
    title="Climate Space for Desert Iguana",
    legend=:bottomright)

plot!(ustrip.(u"W/m^2", Qa_MAX_2), Tas_C, label="V = 2 m/s (Max)")
plot!(ustrip.(u"W/m^2", Qa_MIN_free), Tas_C, label="V = 0.1 m/s (Min)", linestyle=:dash)
plot!(ustrip.(u"W/m^2", Qa_MIN_2), Tas_C, label="V = 2 m/s (Min)", linestyle=:dash)
plot!(ustrip.(u"W/m^2", Qa_night), Tas_C, label="Clear night", color=:black, linestyle=:dot)
plot!(ustrip.(u"W/m^2", Qa_day), Tas_C, label="Clear day (sun)", color=:orange, linestyle=:dot)
plot!(ustrip.(u"W/m^2", Qa_blackbody), Tas_C, label="Blackbody", color=:red, linestyle=:dot)






# ============================================================ 
# Q2 Determine the climate-space for a masked shrew, Sarex cinereus,
# with the characteristics E = 1.0 and body diameter D = 0.018 m. At
# thermal maximum, M = 140 W m-2, 1= 0.085 m20C W-l, 1\£ = 28 \v m-2,
# and Tb = 41°C. At thermal minimum, M = 349 W m-2, 1= 0.125 m2°C
# W-l, 1\£ = 3.5 W m-2, and Tb = 37.5°C.
# ============================================================

# Define the parameters for the masked shrew
epsilon = 1 # Emissivity
sigma = 5.67e-8 # Stefan-Boltzmann constant
D = 0.018 # Body diameter
I_hot = 0.085 # Insulation quality
I_cold = 0.125 # Insulation quality
M_max = 140 # Metrate at Thermal maximum 
M_min = 349 # Metrate at Thermal minimum
lambda_E_max = 28 # Evaporative cooling rate
lambda_E_min = 3.5 # Evaporative cooling rate
Tb_max = 41 # Body temperature
Tb_min = 37.5 # Body temperature
Tas = [-60:1:60;] # Air temperature range


# Smarter calcs for max and min Qa at V = 0.1 m/s and V = 2 m/s
Qa_MAX_free = [(epsilon * sigma * (Tb_max + 273 - I*(M_max - lambda_E_max))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_max - Ta - I*(M_max - lambda_E_max)) for Ta in Tas]
Qa_MIN_free = [(epsilon * sigma * (Tb_min + 273 - I*(M_min - lambda_E_min))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_min - Ta - I*(M_min - lambda_E_min)) for Ta in Tas]


            

# Environmental radiation lines
alpha_s = 0.8   # Solar absorptance
S_solar = 400.0 # Direct + diffuse solar for clear day (W/m²)

# Clear night: sky radiates as blackbody ~20°C cooler than air, ground at Ta
# Animal receives ~half from sky, half from ground
Qa_night = [0.5 * sigma * (Ta + 273.0)^4 + 0.5 * sigma * (Ta + 273.0)^4 for Ta in Tas]

# Clear day with sun: same longwave + absorbed solar
Qa_day = [0.5 * sigma * (Ta + 273.0)^4 + 0.5 * sigma * (Ta + 273.0)^4 + alpha_s * S_solar for Ta in Tas]

# Blackbody: all surroundings at Ta (maximum possible longwave environment)
Qa_blackbody = [sigma * (Ta + 273.0)^4 for Ta in Tas]

using Plots


plot(Qa_MAX_free, Tas, label="V = 0.1 m/s (Max)", xlims=(0, 1200), ylims=(-10, 70),
    xlabel="Radiation (W/m²)", 
    ylabel="Air Temperature (°C)", 
    title="Climate Space for the masked shrew, (Sarex cinereus)", 
    legend=:bottomright)
    plot!(Qa_MIN_free, Tas, label="V = 0.1 m/s (Min)")
    plot!(Qa_night, Tas, label="Clear night", color=:black, linestyle=:dot)
    plot!(Qa_day, Tas, label="Clear day (sun)", color=:orange, linestyle=:dot)

#end





# ============================================================ 
# Q3, what if the shrew had no airflow not even free convection?
# ============================================================

# Define the parameters for the masked shrew
epsilon = 1 # Emissivity
sigma = 5.67e-8 # Stefan-Boltzmann constant
D = 0.018 # Body diameter
I_hot = 0.085 # Insulation quality
I_cold = 0.125 # Insulation quality
M_max = 140 # Metrate at Thermal maximum 
M_min = 349 # Metrate at Thermal minimum
lambda_E_max = 28 # Evaporative cooling rate
lambda_E_min = 3.5 # Evaporative cooling rate
Tb_max = 41u"°C" # Body temperature
Tb_max = uconvert(u"K", Tb_max)
Tb_min = 37.5 # Body temperature
Tas = [-60:1:60;] # Air temperature range
V = 0.0 # No airflow

# Smarter calcs for max and min Qa at V = 0.0 m/s
Qa_MAX_free = [(epsilon * sigma * (Tb_max + 273 - I*(M_max - lambda_E_max))^4) + 
            ((k1 * (V^0.5 / D^0.5)) * (Tb_max - Ta - I*(M_max - lambda_E_max))) for Ta in Tas]
Qa_MIN_free = [(epsilon * sigma * (Tb_min + 273 - I*(M_min - lambda_E_min))^4) + 
            ((k1 * (V^0.5 / D^0.5)) * (Tb_min - Ta - I*(M_min - lambda_E_min))) for Ta in Tas]





# Environmental radiation lines
alpha_s = 0.8   # Solar absorptance
S_solar = 400.0 # Direct + diffuse solar for clear day (W/m²)

# Clear night: sky radiates as blackbody ~20°C cooler than air, ground at Ta
# Animal receives ~half from sky, half from ground
Qa_night = [0.5 * sigma * (Ta + 273.0)^4 + 0.5 * sigma * (Ta + 273.0)^4 for Ta in Tas]

# Clear day with sun: same longwave + absorbed solar
Qa_day = [0.5 * sigma * (Ta + 273.0)^4 + 0.5 * sigma * (Ta + 273.0)^4 + alpha_s * S_solar for Ta in Tas]

# Blackbody: all surroundings at Ta (maximum possible longwave environment)
Qa_blackbody = [sigma * (Ta + 273.0)^4 for Ta in Tas]

using Plots


plot(Qa_MAX_free, Tas, label="V = 0.0 m/s (Max)", xlims=(0, 1200), ylims=(-10, 70),
    xlabel="Radiation (W/m²)", 
    ylabel="Air Temperature (°C)", 
    title="Climate Space for the masked shrew, (Sarex cinereus)", 
    legend=:bottomright)
    plot!(Qa_MIN_free, Tas, label="V = 0.0 m/s (Min)")
    plot!(Qa_night, Tas, label="Clear night", color=:black, linestyle=:dot)
    plot!(Qa_day, Tas, label="Clear day (sun)", color=:orange, linestyle=:dot)
    plot!(Qa_blackbody, Tas, label="Blackbody", colour=:red, linestyle=:dot)
#end





