#include <iostream>
#include <cmath>
#include <stdio.h>

//using namespace std;

double arctan(double x, int terms) {
    double rezultat = 0.0;
    for (int n = 0; n < terms; ++n) {
        rezultat += pow(-1, n) * pow(x, 2 * n + 1) / (2 * n + 1);
    }
    return rezultat;
}

double integral(double x) {
    return exp(3 * x) * arctan(x / 2.0, 20);
}

double trapeznaMetoda(double a, double b, int n) {
    double h = (b - a) / n;
    double rezultat = 0.5 * (integral(a) + integral(b));

    for (int i = 1; i < n; ++i) {
        double x_i = a + i * h;
        rezultat += integral(x_i);
    }

    return h * rezultat;
}

int main() {
    double a = 0.0;
    double b = M_PI / 4.0;
    int n = 1000;

    double priblizekIntegrala = trapeznaMetoda(a, b, n);

    std::cout << "Priblizna vrednost integrala: " << priblizekIntegrala << std::endl;

    return 0;
}