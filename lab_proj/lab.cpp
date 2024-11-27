#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#define PI 3.14159265359
#define N 4 // Кількість змінних стану

using namespace std;

// Константи
const double R1 = 16.2;       // Опір первинної обмотки
const double R21 = 15.8;      // Опір першої половини вторинної обмотки
const double R22 = 15.8;      // Опір другої половини вторинної обмотки
const double RH = 5700.0;     // Опір навантаження
const double C = 30e-6;       // Ємність конденсатора
const double alpha1 = 5.0;    // Індуктивність первинної обмотки
const double alpha21 = 10.0;  // Індуктивність вторинної обмотки (половина 1)
const double alpha22 = 10.0;  // Індуктивність вторинної обмотки (половина 2)
const double Um = 311.0;      // Амплітудна напруга
const double omega = 314.1593;// Кутова частота

// Функція для моделювання рівнянь стану
void dydt(double t, double y[], double dy[], int k1, int k2) {
    double g1 = k1 * alpha1;
    double g21 = k1 * alpha21;
    double g22 = k2 * alpha22;

    double i1 = g1 * y[0] - y[1];
    double i21 = g21 * y[0] - y[2];
    double i22 = g22 * y[0] - y[3];

    dy[0] = Um * sin(omega * t) - R1 * i1;
    dy[1] = i21 - i22 - (y[1] / RH);
    dy[2] = -y[2] / C + (i21 / alpha21);
    dy[3] = -y[3] / C + (i22 / alpha22);
}

// Функція інтегрування методом Рунге-Кутти
void runge_kutta(double h, double t_end) {
    double t = 0.0;
    double y[N] = {0.0, 0.0, 0.0, 0.0};
    double dy[N] = {0.0};

    ofstream file1("i1.txt"), file2("i21.txt"), file3("i22.txt"), file4("uC.txt");

    while (t < t_end) {
        int k1 = (y[1] > 0) ? 1 : 0;
        int k2 = (y[3] > 0) ? 1 : 0;

        dydt(t, y, dy, k1, k2);

        double k1_vals[N], k2_vals[N], k3_vals[N], k4_vals[N];
        double temp_y[N];

        for (int i = 0; i < N; i++) {
            k1_vals[i] = h * dy[i];
            temp_y[i] = y[i] + 0.5 * k1_vals[i];
        }

        dydt(t + 0.5 * h, temp_y, dy, k1, k2);
        for (int i = 0; i < N; i++) {
            k2_vals[i] = h * dy[i];
            temp_y[i] = y[i] + 0.5 * k2_vals[i];
        }

        dydt(t + 0.5 * h, temp_y, dy, k1, k2);
        for (int i = 0; i < N; i++) {
            k3_vals[i] = h * dy[i];
            temp_y[i] = y[i] + k3_vals[i];
        }

        dydt(t + h, temp_y, dy, k1, k2);
        for (int i = 0; i < N; i++) {
            k4_vals[i] = h * dy[i];
            y[i] += (k1_vals[i] + 2.0 * k2_vals[i] + 2.0 * k3_vals[i] + k4_vals[i]) / 6.0;
        }

        file1 << setw(12) << t << setw(12) << y[0] << endl;
        file2 << setw(12) << t << setw(12) << y[1] << endl;
        file3 << setw(12) << t << setw(12) << y[2] << endl;
        file4 << setw(12) << t << setw(12) << y[3] << endl;

        t += h;
    }

    file1.close();
    file2.close();
    file3.close();
    file4.close();
}

int main() {
    double h = 1e-6;     // Крок інтегрування
    double t_end = 0.02; // Кінець моделювання (1 період)

    runge_kutta(h, t_end);

    cout << "Simulation completed. Results saved to files." << endl;
    return 0;
}