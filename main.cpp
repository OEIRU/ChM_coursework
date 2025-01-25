#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <set>  // Используется для эффективного хранения уникальных индексов столбцов
#include <clocale>

// Константы, определяющие параметры задачи
constexpr int BIQUAD_NODES = 9;      // Количество узлов в биквадратичном элементе
constexpr int BILIN_NODES = 4;       // Количество узлов в билинейном элементе
constexpr int GAUSS_POINTS = 3;      // Количество точек интегрирования Гаусса
constexpr double EPS = 1e-12;        // Малая величина для проверки на нуль

// Точки и веса квадратуры Гаусса (3 точки)
const double gauss_points[GAUSS_POINTS] = { -0.7745966692, 0.0, 0.7745966692 };
const double gauss_weights[GAUSS_POINTS] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

// Структура, описывающая узел сетки
struct Node {
    double x, y;            // Координаты узла
    int bc_type = 0;        // Тип граничного условия (0 - внутренний узел)
    double bc_value = 0.0;  // Значение граничного условия
    bool is_active = true;  // Флаг активности узла (используется для Дирихле)
};

// Структура, описывающая конечный элемент
struct Element {
    int nodes[BIQUAD_NODES]; // Индексы узлов элемента
    double gamma[BILIN_NODES]; // Параметры материала (например, теплопроводность)
};

// Структура для хранения граничных условий на ребре
struct EdgeBC {
    int node1, node2;       // Индексы узлов, образующих ребро
    int func_id;            // Идентификатор функции для вычисления значения
    double beta = 0.0;      // Параметр для условий Робина (beta)
    double value = 0.0;     // Значение граничного условия
};

// Параметры solve  и допустимая погрешность)
struct SolverParams {
    int max_iter = 100000;  // Макс. итераций
    double tolerance = 1e-8;// Допустимая погрешность
};

// Структура для хранения разреженной матрицы в формате CSR
struct SparseMatrix {
    std::vector<double> diag;    // Диагональные элементы
    std::vector<double> lower;   // Нижние внедиагональные элементы
    std::vector<int> row_ptr;    // Указатели на начало строк в col_idx
    std::vector<int> col_idx;    // Индексы столбцов для нижних элементов
    int size;                    // Размер матрицы (size x size)
};

// Класс, реализующий МКЭ-решатель
class FEMSolver {
private:
    std::vector<Node> nodes;         // Список узлов
    std::vector<Element> elements;   // Список элементов
    std::vector<EdgeBC> bc_first,    // Граничные условия 1-го рода (Дирихле)
        bc_second,    // Граничные условия 2-го рода (Неймана)
        bc_third;     // Граничные условия 3-го рода (Робин)
    SolverParams params;             // Параметры решателя
    SparseMatrix matrix;             // Разреженная матрица СЛАУ
    std::vector<double> rhs;         // Вектор правой части
    std::vector<double> solution;    // Вектор решения

    // Билинейные базисные функции (для 4-узловых элементов)
    double bilinear_basis(double xi, double eta, int i) {
        const double n[BILIN_NODES][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1} };
        return 0.25 * (1 + xi * n[i][0]) * (1 + eta * n[i][1]);
    }

    // Биквадратичные базисные функции (для 9-узловых элементов)
    double biquadratic_basis(double xi, double eta, int i) {
        const double n[BIQUAD_NODES][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1},
                                         {0,-1}, {1,0}, {0,1}, {-1,0}, {0,0} };
        if (i < 4) // Угловые узлы
            return 0.25 * (xi * n[i][0] + 1) * (eta * n[i][1] + 1) * (xi * n[i][0] + eta * n[i][1] - 1);
        if (i < 8) // Средние узлы ребер
            return 0.5 * (1 - xi * xi) * (eta * n[i - 4][1] + 1) + 0.5 * (1 - eta * eta) * (xi * n[i - 4][0] + 1);
        return (1 - xi * xi) * (1 - eta * eta); // Центральный узел
    }

    // Вычисление градиентов биквадратичных функций
    void biquadratic_grad(double xi, double eta, int i, double& dx, double& dy) {
        const double n[BIQUAD_NODES][2] = { {-1,-1}, {1,-1}, {1,1}, {-1,1},
                                          {0,-1}, {1,0}, {0,1}, {-1,0}, {0,0} };

        if (i < 4) { // Угловые узлы
            double xi_i = n[i][0], eta_i = n[i][1];
            dx = 0.25 * (eta_i + 1) * (2 * xi * xi_i + eta_i * eta - 1);
            dy = 0.25 * (xi_i + 1) * (2 * eta * eta_i + xi_i * xi - 1);
        }
        else if (i < 8) { // Средние узлы ребер
            int k = i - 4;
            if (k % 2 == 0) { // Вертикальные ребра (xi = const)
                dx = -xi * (1 + n[k][1] * eta);
                dy = 0.5 * (1 - xi * xi) * n[k][1];
            }
            else { // Горизонтальные ребра (eta = const)
                dx = 0.5 * (1 - eta * eta) * n[k][0];
                dy = -eta * (1 + n[k][0] * xi);
            }
        }
        else { // Центральный узел
            dx = -2 * xi * (1 - eta * eta);
            dy = -2 * eta * (1 - xi * xi);
        }

        // Проверка на NaN (некорректные вычисления)
        if (std::isnan(dx) || std::isnan(dy)) {
            std::cerr << "Ошибка: Нечисловое значение в градиенте! Узел: " << i
                << " Координаты: (" << xi << ", " << eta << ")\n";
            exit(1);
        }
    }

    // Поиск узлов, принадлежащих ребру (n1, n2)
    std::vector<int> get_edge_nodes(int n1, int n2) {
        const int edge_indices[4][3] = { {0,4,1}, {1,5,2}, {2,6,3}, {3,7,0} };

        // Поиск в стандартных ребрах элементов
        for (const auto& elem : elements) {
            for (int e = 0; e < 4; ++e) {
                int a = elem.nodes[edge_indices[e][0]];
                int mid = elem.nodes[edge_indices[e][1]];
                int b = elem.nodes[edge_indices[e][2]];

                // Проверка ребра в обоих направлениях
                if ((a == n1 && b == n2) || (a == n2 && b == n1)) {
                    return { a, mid, b }; // Возврат узлов ребра
                }
            }

            // Дополнительный поиск по всем парам узлов элемента
            for (int i = 0; i < BIQUAD_NODES; ++i) {
                for (int j = i + 1; j < BIQUAD_NODES; ++j) {
                    int a = elem.nodes[i], b = elem.nodes[j];
                    if ((a == n1 && b == n2) || (a == n2 && b == n1)) {
                        return { a, -1, b }; // -1 означает отсутствие среднего узла
                    }
                }
            }
        }
        return {}; // Ребро не найдено
    }

    // Вычисление длины ребра
    double calculate_edge_length(int n1, int n2) {
        return std::hypot(nodes[n2].x - nodes[n1].x, nodes[n2].y - nodes[n1].y);
    }

    // Функция для получения значения граничного условия (пример)
    double get_bc_value(int func_id, double x, double y) {
        // Реальные функции, соответствующие вашей задаче
        switch (func_id) {
        case 1: return x + y;         
        case 2: return 2 * 3.14 * 3.14 * sin(3.14 * x) * sin(3.14 * y);
        default: return 0.0;
        }
    }

public:
    // Инициализация структуры разреженной матрицы
    void init_sparsity() {
        std::vector<std::set<int>> temp(nodes.size()); // Временное хранение индексов столбцов

        // Перебор всех элементов для заполнения шаблона матрицы
        for (const auto& elem : elements) {
            for (int i = 0; i < BIQUAD_NODES; ++i) {
                int row = elem.nodes[i];
                if (!nodes[row].is_active) continue;

                // Сбор индексов столбцов для текущей строки
                for (int j = 0; j < BIQUAD_NODES; ++j) {
                    int col = elem.nodes[j];
                    if (nodes[col].is_active && col <= row) {
                        temp[row].insert(col);
                    }
                }
            }
        }

        // Преобразование шаблона в формат CSR
        matrix.size = nodes.size();
        matrix.row_ptr.push_back(0);
        for (int i = 0; i < matrix.size; ++i) {
            matrix.diag.push_back(0.0);
            // Добавление индексов столбцов для строки i
            matrix.col_idx.insert(matrix.col_idx.end(), temp[i].begin(), temp[i].end());
            matrix.row_ptr.push_back(matrix.col_idx.size());
        }
        matrix.lower.resize(matrix.col_idx.size(), 0.0);
    }

    // Загрузка узлов из файла
    void load_nodes(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Ошибка открытия файла узлов\n";
            return;
        }

        double dummy;
        int n_nodes;
        file >> dummy >> dummy >> dummy >> n_nodes; 

        nodes.resize(n_nodes);
        for (int i = 0; i < n_nodes; ++i) {
            file >> nodes[i].x >> nodes[i].y >> nodes[i].bc_type >> nodes[i].bc_value;
        }
    }

    // Загрузка элементов из файла
    void load_elems(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Ошибка открытия файла элементов\n";
            return;
        }

        int n_elems;
        file >> n_elems;
        elements.resize(n_elems);

        for (auto& elem : elements) {
            for (int i = 0; i < BIQUAD_NODES; ++i) {
                file >> elem.nodes[i];
                if (elem.nodes[i] < 0 || elem.nodes[i] >= nodes.size()) {
                    std::cerr << "Неверный индекс узла: " << elem.nodes[i] << "\n";
                    exit(1);
                }
            }
            // Чтение параметров gamma
            for (int i = 0; i < BILIN_NODES; ++i) file >> elem.gamma[i];
        }
    }

    // Загрузка граничных условий из файла
    void load_bc(const std::string& filename, int bc_type) {
        std::ifstream file(filename);
        if (!file) return;

        int n_edges;
        file >> n_edges;

        for (int i = 0; i < n_edges; ++i) {
            EdgeBC bc;
            if (bc_type == 1 || bc_type == 2) {
                file >> bc.node1 >> bc.node2 >> bc.func_id;
            }
            else if (bc_type == 3) {
                file >> bc.node1 >> bc.node2 >> bc.func_id >> bc.beta;
            }
            // Добавление в соответствующий список условий
            auto& target = (bc_type == 1) ? bc_first :
                (bc_type == 2) ? bc_second : bc_third;
            target.push_back(bc);
        }
    }

    // Применение граничных условий
    void apply_boundary_conditions() {
        solution.resize(nodes.size(), 0.0);
        rhs.resize(nodes.size(), 0.0);

        // Условия Дирихле (используем bc_value из узлов)
        for (const auto& bc : bc_first) {
            auto edge_nodes = get_edge_nodes(bc.node1, bc.node2);
            if (edge_nodes.empty()) {
                std::cerr << "Предупреждение: Ребро " << bc.node1 << "-" << bc.node2 << " не найдено\n";
                continue;
            }

            // Фиксируем узлы и задаем значения из nodes.bc_value
            for (int node : edge_nodes) {
                if (node < 0 || node >= nodes.size()) continue;
                nodes[node].is_active = false;
                solution[node] = nodes[node].bc_value; // Прямое значение из файла
            }
        }

        // Условия Неймана (интегрируем нагрузку)
        for (const auto& bc : bc_second) {
            auto edge_nodes = get_edge_nodes(bc.node1, bc.node2);
            if (edge_nodes.empty()) continue;

            double length = calculate_edge_length(edge_nodes[0], edge_nodes[2]);
            for (int gi = 0; gi < GAUSS_POINTS; ++gi) {
                double xi = gauss_points[gi];
                double weight = gauss_weights[gi];
                double N[3] = { 0.5 * xi * (xi - 1), 1.0 - xi * xi, 0.5 * xi * (xi + 1) };

                for (int i = 0; i < 3; ++i) {
                    int node = edge_nodes[i];
                    if (node >= nodes.size() || !nodes[node].is_active) continue;
                    // Вычисляем значение нагрузки в координатах узла
                    double x = nodes[node].x;
                    double y = nodes[node].y;
                    rhs[node] += get_bc_value(bc.func_id, x, y) * N[i] * weight * (length / 2.0);
                }
            }
        }

        // Условия Робина (матрица + правая часть)
        for (const auto& bc : bc_third) {
            auto edge_nodes = get_edge_nodes(bc.node1, bc.node2);
            if (edge_nodes.empty()) continue;

            double length = calculate_edge_length(edge_nodes[0], edge_nodes[2]);
            for (int gi = 0; gi < GAUSS_POINTS; ++gi) {
                double xi = gauss_points[gi];
                double weight = gauss_weights[gi];
                double N[3] = { 0.5 * xi * (xi - 1), 1.0 - xi * xi, 0.5 * xi * (xi + 1) };

                // Правая часть
                for (int i = 0; i < 3; ++i) {
                    int row = edge_nodes[i];
                    if (row >= nodes.size() || !nodes[row].is_active) continue;
                    double x = nodes[row].x;
                    double y = nodes[row].y;
                    rhs[row] += get_bc_value(bc.func_id, x, y) * N[i] * weight * (length / 2.0);
                }

                // Матрица (beta * N_i * N_j)
                for (int i = 0; i < 3; ++i) {
                    int row = edge_nodes[i];
                    if (row >= nodes.size() || !nodes[row].is_active) continue;

                    for (int j = 0; j < 3; ++j) {
                        int col = edge_nodes[j];
                        if (col >= nodes.size() || !nodes[col].is_active || col > row) continue;

                        double val = bc.beta * N[i] * N[j] * weight * (length / 2.0);
                        if (row == col) {
                            matrix.diag[row] += val;
                        }
                        else {
                            auto it = std::lower_bound(
                                matrix.col_idx.begin() + matrix.row_ptr[row],
                                matrix.col_idx.begin() + matrix.row_ptr[row + 1],
                                col
                            );
                            if (it != matrix.col_idx.end() && *it == col) {
                                int idx = it - matrix.col_idx.begin();
                                matrix.lower[idx] += val;
                            }
                        }
                    }
                }
            }
        }
    }
    // Сборка глобальной матрицы и правой части
    void assemble_system() {
        rhs.assign(nodes.size(), 0.0);
        matrix.diag.assign(nodes.size(), 0.0);
        matrix.lower.assign(matrix.lower.size(), 0.0);

        // Перебор всех элементов
        for (const auto& elem : elements) {
            double Ke[BIQUAD_NODES][BIQUAD_NODES] = { 0 }; // Локальная матрица
            double Fe[BIQUAD_NODES] = { 0 };               // Локальный вектор

            // Интегрирование по Гауссу
            for (int gi = 0; gi < GAUSS_POINTS; ++gi) {
                for (int gj = 0; gj < GAUSS_POINTS; ++gj) {
                    double xi = gauss_points[gi], eta = gauss_points[gj];
                    double w = gauss_weights[gi] * gauss_weights[gj];

                    // Вычисление якобиана преобразования
                    double J[2][2] = { 0 };
                    for (int k = 0; k < BIQUAD_NODES; ++k) {
                        double dx, dy;
                        biquadratic_grad(xi, eta, k, dx, dy);
                        J[0][0] += dx * nodes[elem.nodes[k]].x;
                        J[0][1] += dx * nodes[elem.nodes[k]].y;
                        J[1][0] += dy * nodes[elem.nodes[k]].x;
                        J[1][1] += dy * nodes[elem.nodes[k]].y;
                    }

                    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
                    detJ = std::abs(detJ); // Модуль определителя

                    if (detJ < EPS) { // Пропуск точек с малым якобианом
                        std::cerr << "Пропуск точки интегрирования (" << xi << ", " << eta
                            << ") из-за малого detJ: " << detJ << "\n";
                        continue;
                    }

                    // Вычисление градиентов базисных функций
                    double dN[BIQUAD_NODES][2];
                    for (int k = 0; k < BIQUAD_NODES; ++k) {
                        double dxi, deta;
                        biquadratic_grad(xi, eta, k, dxi, deta);
                        // Преобразование градиентов в глобальные координаты
                        dN[k][0] = (J[1][1] * dxi - J[0][1] * deta) / detJ;
                        dN[k][1] = (-J[1][0] * dxi + J[0][0] * deta) / detJ;
                    }

                    // Вклад в локальную матрицу и вектор
                    for (int i = 0; i < BIQUAD_NODES; ++i) {
                        // Вычисляем глобальные координаты (x, y) для точки интегрирования
                        double x = 0.0, y = 0.0;
                        for (int k = 0; k < BIQUAD_NODES; ++k) {
                            double Nk = biquadratic_basis(xi, eta, k);
                            x += Nk * nodes[elem.nodes[k]].x;
                            y += Nk * nodes[elem.nodes[k]].y;
                        }
                        for (int j = 0; j < BIQUAD_NODES; ++j) {
                            double gamma_avg = 0.0;
                            for (int k = 0; k < BILIN_NODES; ++k) {
                                gamma_avg += elem.gamma[k] * bilinear_basis(xi, eta, k);
                            }
                            Ke[i][j] += gamma_avg * (dN[i][0] * dN[j][0] + dN[i][1] * dN[j][1]) * w * detJ;
                        }

                        Fe[i] += get_bc_value(2, x, y) * biquadratic_basis(xi, eta, i) * w * detJ;
                    }
                }

            }

            // Перенос локальных вкладов в глобальную систему
            for (int i = 0; i < BIQUAD_NODES; ++i) {
                int row = elem.nodes[i];
                if (!nodes[row].is_active) continue;  // Пропуск неактивных

                rhs[row] += Fe[i]; // Добавление в правую часть
                for (int j = 0; j < BIQUAD_NODES; ++j) {
                    int col = elem.nodes[j];
                    if (col > row || !nodes[col].is_active) continue;

                    if (col == row) {
                        matrix.diag[row] += Ke[i][j]; // Диагональный элемент
                    }
                    else {
                        // Поиск позиции в разреженной матрице
                        auto it = std::lower_bound(
                            matrix.col_idx.begin() + matrix.row_ptr[row],
                            matrix.col_idx.begin() + matrix.row_ptr[row + 1],
                            col
                        );
                        if (it != matrix.col_idx.end() && *it == col) {
                            int idx = it - matrix.col_idx.begin();
                            matrix.lower[idx] += Ke[i][j]; // Внедиагональный элемент
                        }
                    }
                }
            }
        }
    }

    // Решение СЛАУ методом сопряженных градиентов
    void solve() {
        solution.resize(nodes.size(), 0.0);
        std::vector<double> r = rhs; // Вектор невязки
        std::vector<double> p = r;   // Направление поиска
        std::vector<double> Ap(nodes.size(), 0.0); // Вектор A*p

        double rr_old = 0.0;         // Квадрат нормы невязки
        bool converged = false;      // Флаг сходимости

        // Вычисление начальной невязки
        for (size_t i = 0; i < r.size(); ++i) {
            if (std::isnan(r[i]) || std::isinf(r[i])) {
                std::cerr << "Ошибка: Некорректное значение в невязке (r)!\n";
                exit(1);
            }
            rr_old += r[i] * r[i];
        }

        if (std::sqrt(rr_old) < params.tolerance) {
            std::cout << "Решение уже сошлось (начальная невязка).\n";
            return;
        }

        // Итерации метода сопряженных градиентов
        for (int iter = 0; iter < params.max_iter; ++iter) {
            // Вычисление Ap = A * p
            std::fill(Ap.begin(), Ap.end(), 0.0);
            for (int i = 0; i < matrix.size; ++i) {
                if (!nodes[i].is_active) continue;
                Ap[i] += matrix.diag[i] * p[i];
                for (int j = matrix.row_ptr[i]; j < matrix.row_ptr[i + 1]; ++j) {
                    Ap[i] += matrix.lower[j] * p[matrix.col_idx[j]];
                    // Симметричная часть для верхнего треугольника
                    if (matrix.col_idx[j] < i && nodes[matrix.col_idx[j]].is_active) {
                        Ap[matrix.col_idx[j]] += matrix.lower[j] * p[i];
                    }
                }
            }

            // Вычисление шага alpha
            double pAp = 0.0;
            for (size_t i = 0; i < p.size(); ++i) {
                if (!nodes[i].is_active) continue;
                pAp += p[i] * Ap[i];
            }

            if (std::abs(pAp) < EPS) {
                std::cerr << "Ошибка: Деление на ноль (pAp = " << pAp << ")\n";
                break;
            }
            double alpha = rr_old / pAp;

            // Обновление решения и невязки
            for (size_t i = 0; i < solution.size(); ++i) {
                if (!nodes[i].is_active) continue;
                solution[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
                if (std::isnan(solution[i])) {
                    std::cerr << "NaN в решении на итерации " << iter << "\n";
                    exit(1);
                }
            }

            // Проверка нормы новой невязки
            double rr_new = 0.0;
            for (double val : r) rr_new += val * val;

            if (std::sqrt(rr_new) < params.tolerance) {
                std::cout << "Сходимость достигнута за " << iter << " итераций.\n";
                converged = true;
                break;
            }

            // Вычисление коэффициента beta
            if (rr_old < EPS) {
                std::cerr << "Ошибка: Деление на ноль (rr_old = " << rr_old << ")\n";
                break;
            }
            double beta = rr_new / rr_old;

            // Обновление направления поиска
            for (size_t i = 0; i < p.size(); ++i) {
                if (!nodes[i].is_active) continue;
                p[i] = r[i] + beta * p[i];
            }

            rr_old = rr_new;
        }

        if (!converged) {
            std::cerr << "Сходимость не достигнута за " << params.max_iter << " итераций.\n";
        }
    }

    // Сохранение результатов в CSV-файл
    void save_results(const std::string& filename) {
        std::ofstream file(filename);
        file << "x,y,u\n";
        for (size_t i = 0; i < nodes.size(); ++i) {
            file << nodes[i].x << "," << nodes[i].y << "," << solution[i] << "\n";
        }
    }
};

int main() {
    setlocale(LC_ALL, "RU"); 
    FEMSolver solver;
    solver.load_nodes("nodes.txt");     // Загрузка узлов
    solver.load_elems("elems.txt");     // Загрузка элементов
    solver.load_bc("first.txt", 1);     // Условия Дирихле
    solver.load_bc("second.txt", 2);    // Условия Неймана
    solver.load_bc("third.txt", 3);     // Условия Робина
    solver.init_sparsity();             // Инициализация разреженности
    solver.apply_boundary_conditions(); // Применение граничных условий
    solver.assemble_system();           // Сборка системы
    solver.solve();                     // Решение СЛАУ
    solver.save_results("results.csv"); // Сохранение результатов
    return 0;
}