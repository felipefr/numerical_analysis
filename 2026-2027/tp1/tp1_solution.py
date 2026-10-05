import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# TP1 - Python scientifique
# Simulation et Modélisation en Mécanique
# ============================================================


# ============================================================
# 1. Maillage 1D
# ============================================================

print("\n" + "=" * 60)
print("1. MAILLAGE 1D")
print("=" * 60)

L = 1.0
N = 21
dx = L / (N - 1)

x = np.zeros(N)
u = np.zeros(N)

for i in range(N):
    x[i] = i * dx
    u[i] = np.sin(2 * np.pi * x[i])

print("L =", L)
print("N =", N)
print("dx =", dx)

plt.figure()
plt.plot(x, u, '-o')
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Maillage 1D : u(x) = sin(2πx)")
plt.grid()
plt.show()


# ============================================================
# 2. Opérations vectorisées
# ============================================================

print("\n" + "=" * 60)
print("2. OPÉRATIONS VECTORISÉES")
print("=" * 60)

N_vectorise = 41
x_vectorise = np.linspace(0, L, N_vectorise)
u_vectorise = np.sin(2 * np.pi * x_vectorise)

plt.figure()
plt.plot(x, u, '-o', label="N = 21")
plt.plot(x_vectorise, u_vectorise, '-x', label="N = 41")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Comparaison de deux maillages")
plt.legend()
plt.grid()
plt.savefig("maillage_1D.png")
plt.savefig("maillage_1D.pdf")
plt.show()


# ============================================================
# 3. Résolution d'un système linéaire
# ============================================================

print("\n" + "=" * 60)
print("3. RÉSOLUTION D'UN SYSTÈME LINÉAIRE")
print("=" * 60)

K = np.array([
    [4, -1, 0],
    [-1, 4, -1],
    [0, -1, 4]
], dtype=float)

f = np.array([1, 2, 1], dtype=float)

print("K =")
print(K)
print("\nf =")
print(f)

print("\nDimensions :")
print("shape(K) =", np.shape(K))
print("shape(f) =", np.shape(f))

# Résolution directe, sans inversion de K
u_systeme = np.linalg.solve(K, f)

print("\nSolution u =")
print(u_systeme)

# Résidu
residu = np.linalg.norm(K @ u_systeme - f)
print("\nRésidu ||Ku - f|| =", residu)

# Inverse de K
K_inv = np.linalg.inv(K)

print("\nK^(-1) =")
print(K_inv)

# Vérification KK^(-1) = K^(-1)K = I
I3 = np.eye(3)

print("\n||KK^(-1) - I3|| =", np.linalg.norm(K @ K_inv - I3))
print("||K^(-1)K - I3|| =", np.linalg.norm(K_inv @ K - I3))

# Vérification de la deuxième forme de résolution
u_inverse = K_inv @ f
residu_inverse = np.linalg.norm(u_systeme - u_inverse)

print("\nRésidu ||u - K^(-1)f|| =", residu_inverse)


# ============================================================
# 4. Oscillateur masse-ressort et subplots
# ============================================================

print("\n" + "=" * 60)
print("4. OSCILLATEUR MASSE-RESSORT")
print("=" * 60)

m = 1.0
k = 100.0
x0 = 0.01

omega = np.sqrt(k / m)

# Plusieurs périodes
T = 2 * np.pi / omega
t = np.linspace(0, 3 * T, 1000)

x_osc = x0 * np.cos(omega * t)
v_osc = -x0 * omega * np.sin(omega * t)
a_osc = -x0 * omega**2 * np.cos(omega * t)

Ec = 0.5 * m * v_osc**2
Ep = 0.5 * k * x_osc**2
Em = Ec + Ep

print("m =", m, "kg")
print("k =", k, "N/m")
print("x0 =", x0, "m")
print("omega =", omega, "rad/s")
print("Énergie mécanique minimale =", np.min(Em))
print("Énergie mécanique maximale =", np.max(Em))
print("Variation maximale d'énergie =", np.max(Em) - np.min(Em))

fig, axes = plt.subplots(2, 2, figsize=(10, 7))

axes[0, 0].plot(t, x_osc)
axes[0, 0].set_title("Position")
axes[0, 0].set_xlabel("t (s)")
axes[0, 0].set_ylabel("x(t) (m)")
axes[0, 0].grid()

axes[0, 1].plot(t, v_osc)
axes[0, 1].set_title("Vitesse")
axes[0, 1].set_xlabel("t (s)")
axes[0, 1].set_ylabel("v(t) (m/s)")
axes[0, 1].grid()

axes[1, 0].plot(t, a_osc)
axes[1, 0].set_title("Accélération")
axes[1, 0].set_xlabel("t (s)")
axes[1, 0].set_ylabel("a(t) (m/s²)")
axes[1, 0].grid()

axes[1, 1].plot(t, Ec, label="Ec")
axes[1, 1].plot(t, Ep, label="Ep")
axes[1, 1].plot(t, Em, label="Em")
axes[1, 1].set_title("Énergies")
axes[1, 1].set_xlabel("t (s)")
axes[1, 1].set_ylabel("Énergie (J)")
axes[1, 1].legend()
axes[1, 1].grid()

fig.suptitle("Oscillateur masse-ressort")
fig.tight_layout()
plt.show()


# ============================================================
# 5. Visualisation d'un champ 2D
# ============================================================

print("\n" + "=" * 60)
print("5. VISUALISATION D'UN CHAMP 2D")
print("=" * 60)

x_2d = np.linspace(0, 1, 100)
y_2d = np.linspace(0, 1, 100)
X, Y = np.meshgrid(x_2d, y_2d)

U = np.sin(np.pi * X) * np.sin(np.pi * Y)

plt.figure()
contours = plt.contourf(X, Y, U, levels=30)
plt.colorbar(contours, label="u(x,y)")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Champ scalaire u(x,y)")
plt.axis("equal")
plt.show()


# ============================================================
# 6. Dérivation numérique
# ============================================================

print("\n" + "=" * 60)
print("6. DÉRIVATION NUMÉRIQUE")
print("=" * 60)


N = 101
dx = L / (N - 1)
x = np.linspace(0, L, N)

u = np.sin(2 * np.pi * x)
du_exact = 2 * np.pi * np.cos(2 * np.pi * x)

du_A = np.zeros(N)
du_B = np.zeros(N)
du_C = np.zeros(N)

# Méthode i) : différence avant, valable pour i = 0,...,N-2
du_A[:-1] = (u[1:] - u[:-1]) / dx
du_A[-1] = du_exact[-1]  # On utilise la dérivée exacte au dernier point

# Méthode ii) : différence arrière, valable pour i = 1,...,N-1
du_B[1:] = (u[1:] - u[:-1]) / dx
du_B[0] = du_exact[0]  # On utilise la dérivée exacte au premier point

# Méthode iii) : différence centrée, valable pour i = 1,...,N-2
du_C[1:-1] = (u[2:] - u[:-2]) / (2 * dx)
du_C[-1] = du_exact[-1]  # On utilise la dérivée exacte au dernier point
du_C[0] = du_exact[0]  # On utilise la dérivée exacte au premier point

print("Méthode i)  : i = 0,...,N-2")
print("Méthode ii) : i = 1,...,N-1")
print("Méthode iii): i = 1,...,N-2")

plt.figure(figsize=(10, 6))
plt.plot(x, du_A, '-o', markevery=5, label="Différence avant")
plt.plot(x, du_B, '-x', markevery=5, label="Différence arrière")
plt.plot(x, du_C, '-+', markevery=5, label="Différence centrée")
plt.plot(x, du_exact, '--', label="Analytique")

plt.xlabel("x")
plt.ylabel("u'(x)")
plt.title("Approximation numérique de la dérivée")
plt.legend()
plt.grid()
plt.show()

# ============================================================
# 7. Étude de convergence
# ============================================================

print("\n" + "=" * 60)
print("7. ÉTUDE DE CONVERGENCE")
print("=" * 60)

N_list = [21, 41, 81, 161, 321]

erreur_A_list = []
erreur_C_list = []
dx_list = []

for N in N_list:
    dx = L / (N - 1)
    x = np.linspace(0, L, N)

    u = np.sin(2 * np.pi * x)
    du_exact = 2 * np.pi * np.cos(2 * np.pi * x)

    du_A = np.zeros(N)
    du_B = np.zeros(N)
    du_C = np.zeros(N)

    # Méthode i) : différence avant, valable pour i = 0,...,N-2
    du_A[:-1] = (u[1:] - u[:-1]) / dx
    du_A[-1] = du_exact[-1]  # On utilise la dérivée exacte au dernier point

    # Méthode iii) : différence centrée, valable pour i = 1,...,N-2
    du_C[1:-1] = (u[2:] - u[:-2]) / (2 * dx)
    du_C[-1] = du_exact[-1]  # On utilise la dérivée exacte au dernier point
    du_C[0] = du_exact[0]  # On utilise la dérivée exacte au premier point


    erreur_A = np.max(np.abs(du_A - du_exact))
    erreur_C = np.max(np.abs(du_C - du_exact))

    dx_list.append(dx)
    erreur_A_list.append(erreur_A)
    erreur_C_list.append(erreur_C)

print("\nRésultats de convergence :")
print(" N       dx             erreur A          erreur C")

for N, dx_conv, err_A, err_C in zip(
    N_list, dx_list, erreur_A_list, erreur_C_list
):
    print(
        f"{N:3d}   {dx_conv:.6e}   "
        f"{err_A:.6e}   {err_C:.6e}"
    )

# Estimation de l'ordre de convergence
ordre_A = np.log(erreur_A_list[-1]/erreur_A_list[-2])/np.log(dx_list[-1]/dx_list[-2])
ordre_C = np.log(erreur_C_list[-1]/erreur_C_list[-2])/np.log(dx_list[-1]/dx_list[-2])

print("\nOrdre de convergence estimé :")
print("Méthode i)   :", ordre_A)
print("Méthode iii) :", ordre_C)

plt.figure(figsize=(8, 6))
plt.loglog(
    dx_list,
    erreur_A_list,
    '-o',
    label="Différence avant"
)
plt.loglog(
    dx_list,
    erreur_C_list,
    '-x',
    label="Différence centrée"
)

plt.xlabel("dx")
plt.ylabel("Erreur maximale")
plt.title("Étude de convergence")
plt.legend()
plt.grid(True, which="both")
plt.show()


# ============================================================
# 8. Méthode de la bissection
# ============================================================

print("\n" + "=" * 60)
print("8. MÉTHODE DE LA BISSECTION")
print("=" * 60)

# Méthode de la bissection (CM1) modifiée
def bissection(f, a, b, tol=1e-10, nmax=100):
    fa = f(a)
    fb = f(b)
    hist = []
    if fa * fb > 0:
        raise ValueError("Pas de changement de signe")

    for k in range(nmax):
        m = (a + b) / 2
        fm = f(m)
        hist.append(m)

        if abs(fm) < tol or (b-a)/2 < tol:
            return m, hist

        if fa * fm < 0:
            b, fb = m, fm
        else:
            a, fa = m, fm

    m = (a + b) / 2
    hist.append(m)

    return m, hist


def f(x):
    return x**3 - x - 2


racine, hist = bissection(f, 1, 2)


hist = np.array(hist)

# Erreur par rapport à la solution finale obtenue
erreur_bissection = np.abs(hist[0:-1] - racine)

print("Solution finale    :", racine)
print("Nombre d'itérations:", len(hist))
print("Dernière erreur     :", erreur_bissection[-1])

plt.figure()
plt.semilogy(erreur_bissection, '-o')
plt.xlabel("Itération")
plt.ylabel("Erreur")
plt.title("Convergence de la méthode de la bissection")
plt.grid()
plt.show()


print("\nTP terminé.")
