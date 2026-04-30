

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
Tb_min = 30 # Body temperature

# Define the environment parameters
Vs = [0.1, 2] # Wind speed
#Ta = 25 # Air temperature

# Maximum temp
k1 = 3.89
A = epsilon * sigma * ((Tb_max + 273) - I*(M_max - lambda_E_max))^4
B1 = k1 * (sqrt.(Vs)/sqrt(D))
B2 = (Tb_max - I*(M_max - lambda_E_max))


function linear_Qa(; σ = 5.67e-8, Tb = 42.5, V = 0.10, D = 0.05, k1 = 3.89)
    # Radiation term
    A_rad = σ * (Tb + 273)^4

    # Convection coefficient
    h = k1 * sqrt(V / D)

    # Qa = A_rad + h * (Tb - Ta)
    # Qa = (A_rad + h*Tb) + (-h)*Ta
    A = A_rad + h * Tb
    B = -h

    return A, B
end

A, B = linear_Qa()

println("Qa = $(round(A, digits = 2)) + $(round(B, digits = 2))Ta")





















