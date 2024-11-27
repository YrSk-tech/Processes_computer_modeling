import matplotlib
matplotlib.use('Agg')  # Використовуємо неінтерактивний бекенд
import matplotlib.pyplot as plt

# Завантаження даних з файлів
def load_data(filename):
    time = []
    values = []
    with open(filename, 'r') as file:
        for line in file:
            t, val = map(float, line.split())
            time.append(t)
            values.append(val)
    return time, values

# Завантаження даних
time_i1, i1 = load_data("i1.txt")
time_i21, i21 = load_data("i21.txt")
time_i22, i22 = load_data("i22.txt")
time_uc, uc = load_data("uC.txt")

# Побудова графіків
plt.figure(figsize=(12, 8))

# Графік i1
plt.subplot(2, 2, 1)
plt.plot(time_i1, i1, label="i1(t)", color="blue")
plt.xlabel("Time [s]")
plt.ylabel("i1 [A]")
plt.title("Primary Current (i1)")
plt.grid(True)
plt.legend()

# Графік i21
plt.subplot(2, 2, 2)
plt.plot(time_i21, i21, label="i21(t)", color="red")
plt.xlabel("Time [s]")
plt.ylabel("i21 [A]")
plt.title("Secondary Current (i21)")
plt.grid(True)
plt.legend()

# Графік i22
plt.subplot(2, 2, 3)
plt.plot(time_i22, i22, label="i22(t)", color="green")
plt.xlabel("Time [s]")
plt.ylabel("i22 [A]")
plt.title("Secondary Current (i22)")
plt.grid(True)
plt.legend()

# Графік uC
plt.subplot(2, 2, 4)
plt.plot(time_uc, uc, label="uC(t)", color="purple")
plt.xlabel("Time [s]")
plt.ylabel("uC [V]")
plt.title("Capacitor Voltage (uC)")
plt.grid(True)
plt.legend()

# Збереження графіків у файл
output_file = "graphs.png"
plt.tight_layout()
plt.savefig(output_file, dpi=300)

print(f"Графіки збережено у файл: {output_file}")
