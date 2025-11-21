import numpy as np

# -----------------------------
# Datos de la red
# -----------------------------
Y = np.array([
    [-41.684j,  24.876j,   0j,       4.355j,   12.453j],
    [ 24.876j, -31.099j,   0j,       6.223j,    0j     ],
    [  0j,       0j,     -8.71j,     4.355j,    4.355j ],
    [  4.355j,   6.223j,   4.355j, -14.933j,    0j     ],
    [ 12.453j,   0j,       4.355j,    0j,     -16.808j ]
], dtype=complex)

num_nodos = Y.shape[0]

# Definición nodos (fijos en esta red)
slack = 0
pv = num_nodos - 1
pq_nodes = np.arange(1, pv)     # 1, 2, 3

# -----------------------------
# Datos eléctricos
# -----------------------------
V_slack = np.array([1.0977 + 0j])
V_pq = np.ones(len(pq_nodes), dtype=complex)
V_pv = np.array([1 + 0j])

V = np.concatenate((V_slack, V_pq, V_pv))

S_pq = np.array([-0.125 - 0.07628j,
                 -0.66  - 0.36688j,
                 -0.25  - 0.15j])

P_pv = 0.8

numero_iteraciones = 7

# -----------------------------
# Funciones
# -----------------------------
def calcular_S_pv(V, Y, P_pv):
    """Calcula Q del nodo PV a partir del balance."""
    I = Y[pv, :] @ V
    S_pv = P_pv - (-np.imag(V[pv] * np.conjugate(I))) * 1j
    return np.array([S_pv])

def calcular_S_slack(V, Y):
    I = Y[slack, :] @ V
    return V[slack] * np.conjugate(I)

def actualizar_voltajes(V, S_pq, S_pv, Y):
    """Gauss-Seidel para nodos PQ y PV (excepto slack)."""

    S = np.concatenate((S_pq, S_pv))

    for idx, nodo in enumerate(pq_nodes):
        suma = (Y[nodo, :] * V).sum() - Y[nodo, nodo] * V[nodo]
        V[nodo] = (np.conjugate(S[idx]) / np.conjugate(V[nodo]) - suma) / Y[nodo, nodo]

    # Nodo PV
    suma = (Y[pv, :] * V).sum() - Y[pv, pv] * V[pv]
    V[pv] = (np.conjugate(S[-1]) / np.conjugate(V[pv]) - suma) / Y[pv, pv]

    return V

# -----------------------------
# Iteraciones
# -----------------------------
V_prev = V.copy()

for n in range(numero_iteraciones + 1):

    # -----------------------------
    # Cálculo de potencias
    # -----------------------------
    S_pv = calcular_S_pv(V, Y, P_pv)
    S_slack = calcular_S_slack(V, Y)
    S_tot = np.concatenate((S_pq, S_pv))

    # -----------------------------
    # Mostrar resultados
    # -----------------------------
    print(f"\nIteración {n}")

    print(f"Nodo Slack  V = {V[slack]:.4f}   S = {S_slack:.4f}")
    for i, nodo in enumerate(pq_nodes):
        print(f"Nodo PQ {nodo+1}  V = {V[nodo]:.4f}   S = {S_pq[i]:.4f}")
    print(f"Nodo PV {pv+1}   V = {V[pv]:.4f}   S = {S_pv[0]:.4f}")

    # Error
    if n > 0:
        error = np.sum(np.abs(V - V_prev))
        print(f"Error = {error:.5f}")

    V_prev = V.copy()

    # -----------------------------
    # Actualización Gauss-Seidel
    # -----------------------------
    V = actualizar_voltajes(V, S_pq, S_pv, Y)

    # Corrección de módulo del PV
    V[pv] = V_pv[0] * V[pv] / np.abs(V[pv])
