#(100 * pow(10, 0.66077 + (7.5 * T.Bottom[BOTTOM(i,0,k)]) / (237.3 + T.Bottom[BOTTOM(i,0,k)])) * (MW_Water / 1000)) / (8.314 * (273.15 + T.Bottom[BOTTOM(i,0,k)]));

Temp_a = 22.4;
Temp_b = 20.0;

a = (100 * pow(10, 0.66077 + (7.5 * Temp_a) / (237.3 + Temp_a)) * (18 / 1000)) / (8.314 * (273.15 + Temp_a))
b = (100 * pow(10, 0.66077 + (7.5 * Temp_b) / (237.3 + Temp_b)) * (18 / 1000)) / (8.314 * (273.15 + Temp_b))

print(a)
print(b)