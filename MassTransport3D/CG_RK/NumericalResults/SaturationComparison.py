# COMPARISON OF WATER SATURATION PRESSURES VALUES
import math

Temperature = 25 # Celsius

# SCIENCE DIRECT
Ps_SD = 100 * pow(10, 0.66077 + (7.5 * Temperature) / (237.3 + Temperature)) # Pa
Ps_SD2 = 1000 * (100 * (0.004516 + 0.0007178 * Temperature - 2.649 * pow(10,-6.0) * pow(Temperature, 2.0) + 6.944 * pow(10, -7.0) * pow(Temperature, 3.0)))

# PDF PROPIEDADES GASES

T = Temperature + 273.15
Ps_FOR = math.exp(-5.8002206 * pow(10, 3.0) * pow(T, -1.0) + 1.3914993 - 4.8640239 * pow(10, -2.0) * T + 4.1764768 * pow(10, -5.0) * pow(T, 2.0) - 1.4452093 * pow(10, -8.0) * pow(T, 3.0) + 6.5459673 * math.log(T))


print('Water saturation pressures: \n')
print('Pressure Science Direct 1: {0} Pa'.format(Ps_SD))
print('Pressure Science Direct 2: {0} Pa'.format(Ps_SD2))
print('Pressure Formulae: {0} Pa'.format(Ps_FOR))
print('Relative Difference: {0} %'.format(100 * abs(Ps_FOR - Ps_SD)/Ps_FOR))
