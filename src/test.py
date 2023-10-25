import numpy as np

def S_elem(A, B, K, w=None):
    dim = A.shape[0]
    coul = 0
    D = A + B
    R = np.linalg.inv(D)
    M0 = (np.pi**dim / np.linalg.det(D))**(3/2)
    trace = np.trace(B @ K @ A @ R)
    if w is not None:
        for k in range(len(w)):
            beta = 1 / (w[k].T @ R @ w[k])
            if k == 2:
                coul += 2 * np.sqrt(beta / np.pi) * M0
            else:
                coul -= 2 * np.sqrt(beta / np.pi) * M0
        return M0, trace, coul
    else:
        return M0, trace

def transform_list(alphas):
    g_new = [np.ones((1, 1)) * alphas[i] for i in range(len(alphas))]
    return g_new

def S_wave(alphas, K, w=None):
    length = len(alphas)
    alphas = transform_list(alphas)
    kinetic = np.zeros((length, length))
    overlap = np.zeros((length, length))
    coulomb = np.zeros((length, length))
    for i in range(length):
        for j in range(length):
            if j <= i:
                A = alphas[i]
                B = alphas[j]
                M0, trace, coul = S_elem(A, B, K, w)
                R = np.linalg.inv(A + B)
                overlap[i, j] = M0
                overlap[j, i] = overlap[i, j]
                kinetic[i, j] = 6 * trace * M0
                kinetic[j, i] = kinetic[i, j]
                coulomb[i, j] = coul
                coulomb[j, i] = coulomb[i, j]
    return overlap, kinetic, coulomb

# Example usage
A = np.array([[1.0, 0.0], [0.0, 2.0]])
B = np.array([[0.0, 0.5], [0.5, 0.0]])
K = np.array([[2.0, 0.0], [0.0, 1.0]])
w = [np.array([0.5, 0.3]), np.array([0.2, 0.1])]

overlap, kinetic, coulomb = S_wave([1.0, 2.0], K, w)
print("Overlap:")
print(overlap)
print("Kinetic:")
print(kinetic)
print("Coulomb:")
print(coulomb)
