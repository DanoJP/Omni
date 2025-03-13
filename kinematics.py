import numpy as np

# Parámetros del robot
R = 0.05  # Radio de las ruedas (m)
L = 0.11  # Distancia del centro a las ruedas (m)

# Matriz de cinemática inversa
E = (1/R) * np.array([
    [np.sqrt(3)/2, -0.5, -L],
    [0, 1, -L],
    [-np.sqrt(3)/2, -0.5, -L]
])
