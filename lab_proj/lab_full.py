import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Параметри системи
r1 = 16.2  # Ом
r21 = r22 = 15.8  # Ом
RH = 5.7e3  # Ом
alpha1 = 5.0  # Гн^-1
alpha21 = alpha22 = 10.0  # Гн^-1
C = 30e-6  # Ф
Um = 311  # В
omega = 314.1593  # рад/с
T = 0.02  # період, с
dt = 0.00001  # крок інтегрування, с

# Крива намагнічення
psi1 = 0.3
psi2 = 0.9
phi_psi1 = 0.06
phi_psi2 = 2.2
a0 = 6.8
a1 = 0.2
a2 = 10.0

# Апроксимація кривої намагнічення
def phi(psi):
    psi = np.abs(psi)
    return np.where(
        psi <= psi1,
        phi_psi1,
        np.where(
            psi <= psi2,
            a0 + a1 * psi + a2 * psi**2,
            phi_psi2
        )
    )

# Розрахунок коефіцієнтів
def calc_coefficients(psi, k1, k2):
    alpha_prime = phi(psi) / psi if psi != 0 else phi_psi1
    delta = alpha_prime + alpha1 + k1 * alpha21 + k2 * alpha22
    g1 = alpha1 / delta
    g21 = k1 * alpha21 / delta
    g22 = -k2 * alpha22 / delta
    return g1, g21, g22

# Система рівнянь
def dydt(t, y, k1, k2):
    psi, i21, i22, uC = y
    u1 = Um * np.sin(omega * t)
    g1, g21, g22 = calc_coefficients(psi, k1, k2)

    dpsi_dt = g1 * (u1 - r1 * (phi(psi) * psi - i21 - i22)) \
              + g21 * (-uC - r21 * i21) \
              + g22 * (-uC - r22 * i22)
    di21_dt = -alpha21 * g1 * (u1 - r1 * (phi(psi) * psi - i21 - i22)) \
              + alpha21 * (1 - g21) * (-uC - r21 * i21)
    di22_dt = alpha22 * g1 * (u1 - r1 * (phi(psi) * psi - i21 - i22)) \
              + alpha22 * (1 + g22) * (-uC - r22 * i22)
    duC_dt = (i21 + i22 - uC / RH) / C

    return np.array([dpsi_dt, di21_dt, di22_dt, duC_dt])

# Метод Рунге-Кутта
def runge_kutta_step(t, y, dt, k1, k2):
    k1_vals = dydt(t, y, k1, k2)
    k2_vals = dydt(t + dt / 2, y + dt * k1_vals / 2, k1, k2)
    k3_vals = dydt(t + dt / 2, y + dt * k2_vals / 2, k1, k2)
    k4_vals = dydt(t + dt, y + dt * k3_vals, k1, k2)
    return y + dt * (k1_vals + 2 * k2_vals + 2 * k3_vals + k4_vals) / 6

# Інтегрування рівнянь
time = np.arange(0, 4 * T, dt)
y = np.zeros((len(time), 4))  # Массив для значень [psi, i21, i22, uC]
y[0] = [0.0, 0.0, 0.0, 0.0]  # Початкові умови

k1, k2 = 1, 0  # Початковий стан вентилів

for i in range(1, len(time)):
    if k1 == 1 and y[i - 1, 1] <= 0:
        k1 = 0
    elif k1 == 0 and -y[i - 1, 0] - y[i - 1, 3] > 0:
        k1 = 1
    if k2 == 1 and y[i - 1, 2] <= 0:
        k2 = 0
    elif k2 == 0 and y[i - 1, 0] - y[i - 1, 3] > 0:
        k2 = 1

    y[i] = runge_kutta_step(time[i - 1], y[i - 1], dt, k1, k2)

# Візуалізація результатів
plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
plt.plot(time, y[:, 1], label="i21")
plt.plot(time, y[:, 2], label="i22")
plt.title("Вторинні струми")
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(time, y[:, 0], label="Ψ")
plt.title("Потокозчеплення")
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(time, y[:, 3], label="uC")
plt.title("Напруга на конденсаторі")
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(time, phi(y[:, 0]) * y[:, 0], label="i1")
plt.title("Первинний струм")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig("results.png")
