#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <locale>

int load_net(int xwn, std::vector<double> Xw, int ywn, std::vector<double> Yw, int L, std::vector<double> X_kord, std::vector<double> Y_kord) {
    
    return 0;
}

int main() {
    setlocale(LC_ALL, "Russian");

    // Параметры сетки
    int xwn = 3;                              // Количество узлов по оси x
    std::vector<double> Xw = {0.0, 1.0, 2.0}; // Координаты узлов по оси x
    int ywn = 3;                              // Количество узлов по оси y
    std::vector<double> Yw = {0.0, 1.0, 2.0}; // Координаты узлов по оси y
    int L = 1;                                // Количество прямоугольных эл.
    std::vector<double> X_kord = {0.0, 2.0};   // диапозон сетки по x
    std::vector<double> Y_kord = {0.0, 2.0};   // диапозон сетки по y

    load_net(xwn, Xw, ywn, Yw, L, X_kord, Y_kord);
    
    // Параметры краевых условий
    int dirihle_count = 2;
    // тут нужен вектор значений для каждого index (их 2)

    int neiman_count = 1;
    // тут нужен вектор значений для каждого index (их 1)

    int robin_count = 1;
    // тут нужен вектор значений для каждого index (их 1)

    return 0;
}