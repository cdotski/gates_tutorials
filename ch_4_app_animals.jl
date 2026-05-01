

# ===========================================================
# Q1 Determine the climate-space for a desert iguana, Dipsosaurus dorsalis
# ===========================================================
using Unitful

# Define the parameters for the desert iguana
epsilon = 1 # Emissivity
sigma = 5.67e-8 # Stefan-Boltzmann constant
D = 0.015 # Body diameter
I = 0.005 # Insulation quality
M_max = 10.5 # Metrate at Thermal maximum 
M_min = 0.1 # Metrate at Thermal minimum
lambda_E_max = 6.3 # Evaporative cooling rate
lambda_E_min = 0.1 # Evaporative cooling rate
Tb_max = 45 # Body temperature
Tb_min = 3 # Body temperature

# Define the environment parameters
Vs = [0.1, 2] # Wind speed
Tas = [-60:1:60;] # Air temperature range
#Ta = 25 # Air temperature

# Maximum temp
k1 = 3.89

Qa_MAX = (epsilon * sigma * (Tb_max + 273 - I*(M_max - lambda_E_max))^4) + 
            k1 * (2^0.5 / D^0.5) * (Tb_max - Ta - I*(M_max - lambda_E_max))
# Simplifies to 579.66 + 10.04 * (44.979 - Ta)
        # --> 1031.249 - 10.04 * Ta with V = 0.1 m/s
        # --> 2599.98 - 44.91 * Ta with V = 2 m/s
Qa_MIN = (epsilon * sigma * (Tb_min + 273 - I*(M_min - lambda_E_min))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_min - Ta - I*(M_min - lambda_E_min))
# Simplifies to 477.92 + 10.04 * (3 - Ta)
        # --> 508.04 - 10.04 * Ta with V = 0.1 m/s
        # --> 642.77 - 44.91 * Ta with V = 2 m/s

# Hand calculations for max and min Qa at V = 0.1 m/s and V = 2 m/s
Qa_free = [1031.249 - 10.04 * Ta for Ta in Tas] #with V = 0.1 m/s
Qa_2 = [2599.98 - 44.91 * Ta for Ta in Tas] #with V = 2 m/s

Qa_free_min = [508.04 - 10.04 * Ta for Ta in Tas] #with V = 0.1 m/s
Qa_2_min = [642.77 - 44.91 * Ta for Ta in Tas] #with V = 2 m/s  






# Smarter calcs for max and min Qa at V = 0.1 m/s and V = 2 m/s
Qa_MAX_free = [(epsilon * sigma * (Tb_max + 273 - I*(M_max - lambda_E_max))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_max - Ta - I*(M_max - lambda_E_max)) for Ta in Tas]
Qa_MAX_2 = [(epsilon * sigma * (Tb_max + 273 - I*(M_max - lambda_E_max))^4) + 
            k1 * (2^0.5 / D^0.5) * (Tb_max - Ta - I*(M_max - lambda_E_max)) for Ta in Tas]
Qa_MIN_free = [(epsilon * sigma * (Tb_min + 273 - I*(M_min - lambda_E_min))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_min - Ta - I*(M_min - lambda_E_min)) for Ta in Tas]
Qa_MIN_2 = [(epsilon * sigma * (Tb_min + 273 - I*(M_min - lambda_E_min))^4) + 
            k1 * (2^0.5 / D^0.5) * (Tb_min - Ta - I*(M_min - lambda_E_min)) for Ta in Tas]


# Environmental radiation lines
alpha_s = 0.8   # Solar absorptance
S_solar = 400.0 # Direct + diffuse solar for clear day (W/m²)

# Clear night: sky radiates as blackbody ~20°C cooler than air, ground at Ta
# Animal receives ~half from sky, half from ground
Qa_night = [0.5 * sigma * (Ta + 253.0)^4 + 0.5 * sigma * (Ta + 273.0)^4 for Ta in Tas]

# Clear day with sun: same longwave + absorbed solar
Qa_day = [0.5 * sigma * (Ta + 253.0)^4 + 0.5 * sigma * (Ta + 273.0)^4 + alpha_s * S_solar for Ta in Tas]

# Blackbody: all surroundings at Ta (maximum possible longwave environment)
Qa_blackbody = [sigma * (Ta + 273.0)^4 for Ta in Tas]

using Plots


plot(Qa_MAX_free, Tas, label="V = 0.1 m/s (Max)", xlims=(0, 1200), ylims=(-70, 70),
    xlabel="Radiation (W/m²)", 
    ylabel="Air Temperature (°C)", 
    title="Climate Space for Desert Iguana", 
    legend=:bottomright)
    plot!(Qa_MAX_2, Tas, label="V = 2 m/s (Max)")
    plot!(Qa_MIN_free, Tas, label="V = 0.1 m/s (Min)", linestyle=:dash)
    plot!(Qa_MIN_2, Tas, label="V = 2 m/s (Min)", linestyle=:dash)
    plot!(Qa_night, Tas, label="Clear night", color=:black, linestyle=:dot)
    plot!(Qa_day, Tas, label="Clear day (sun)", color=:orange, linestyle=:dot)
    plot!(Qa_blackbody, Tas, label="Blackbody", color=:red, linestyle=:dot)




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
I = 0.085 # Insulation quality
M_max = 140 # Metrate at Thermal maximum 
M_min = 349 # Metrate at Thermal minimum
lambda_E_max = 28 # Evaporative cooling rate
lambda_E_min = 3.5 # Evaporative cooling rate
Tb_max = 41 # Body temperature
Tb_min = 37.5 # Body temperature


# Smarter calcs for max and min Qa at V = 0.1 m/s and V = 2 m/s
Qa_MAX_free = [(epsilon * sigma * (Tb_max + 273 - I*(M_max - lambda_E_max))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_max - Ta - I*(M_max - lambda_E_max)) for Ta in Tas]

Qa_MIN_free = [(epsilon * sigma * (Tb_min + 273 - I*(M_min - lambda_E_min))^4) + 
            k1 * (0.1^0.5 / D^0.5) * (Tb_min - Ta - I*(M_min - lambda_E_min)) for Ta in Tas]


plot(Qa_MAX_free, Tas, label="V = 0.1 m/s (Max)", xlims=(0, 1200), ylims=(-70, 70),
    xlabel="Radiation (W/m²)", 
    ylabel="Air Temperature (°C)", 
    title="Climate Space for Masked Shrew", 
    legend=:bottomright)
    plot!(Qa_MIN_free, Tas, label="V = 0.1 m/s (Min)", linestyle=:dash)
    plot!(Qa_night, Tas, label="Clear night", color=:black, linestyle=:dot)
    plot!(Qa_day, Tas, label="Clear day (sun)", color=:orange, linestyle=:dot)
    plot!(Qa_blackbody, Tas, label="Blackbody", color=:red, linestyle=:dot)





