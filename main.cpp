#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <locale>
#include <string>
#include <algorithm>
#include <stdexcept>

// Константы
constexpr int NUM_NODES_PER_ELEMENT = 9; // Для биквадратичных функций
constexpr int MAX_ITER = 100000;
constexpr double EPSILON = 1e-12;

// Структура для представления разреженной матрицы в формате CSR
struct CSRMatrix {
    int n; // Размерность матрицы
    int nnz; // Количество ненулевых элементов
    std::vector<double> values; // Значения ненулевых элементов
    std::vector<int> col_idx; // Индексы столбцов
    std::vector<int> row_ptr; // Указатели на начало каждой строки

    CSRMatrix(int size = 0) : n(size), nnz(0), row_ptr(size + 1, 0) {}
};

// Структура для представления сетки конечных элементов
struct Mesh {
    int xwn = 0; // Количество узлов по X
    int ywn = 0; // Количество узлов по Y
    std::vector<double> Xw; // Координаты узлов по X
    std::vector<double> Yw; // Координаты узлов по Y
    int L = 0; // Число подобластей
    std::vector<std::vector<int>> W; // Координаты подобластей (левый, правый, нижний, верхний)

    // Глобальные номера узлов для каждого КЭ
    std::vector<int> nvkat; // Номера подобластей для каждого КЭ
    std::vector<int> nvtr; // Глобальные номера узлов для каждого КЭ (одномерный массив)
    std::vector<bool> fict; // Флаг фиктивных узлов
    int n = 0; // Общее количество узлов
    int k = 0; // Количество конечных элементов

    // Координаты узлов
    std::vector<double> xy_x;
    std::vector<double> xy_y;

    // Конструктор по умолчанию
    Mesh() = default;
};

// Структуры для краевых условий
struct BoundaryCondition {
    // Краевые условия 1-го рода (Дирихле)
    std::vector<int> kt1;
    std::vector<std::pair<int, int>> l1;

    // Краевые условия 2-го рода (Неймана)
    std::vector<int> nvk2;
    std::vector<std::pair<int, int>> nvr2;

    // Краевые условия 3-го рода (Робин)
    std::vector<int> nvk3;
    std::vector<std::pair<int, int>> nvr3;
};

// Класс для представления конечного элемента
class FiniteElement {
public:
    std::vector<double> G; // Матрица жесткости
    std::vector<double> C; // Матрица массы
    std::vector<double> Ck; // Матрица для краевых условий
    std::vector<double> Ak; // Локальная матрица
    std::vector<double> Bk; // Локальный вектор

    FiniteElement()
        : G(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        C(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        Ck(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        Ak(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        Bk(NUM_NODES_PER_ELEMENT, 0.0) {
        initialize_matrices();
    }

    void initialize_matrices() {
        // Заполняем G, C и Ck для биквадратичных функций
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
                G[i * NUM_NODES_PER_ELEMENT + j] = calculate_G(i, j);
                C[i * NUM_NODES_PER_ELEMENT + j] = calculate_C(i, j);
                Ck[i * NUM_NODES_PER_ELEMENT + j] = calculate_Ck(i, j);
            }
        }
    }

private:
    // Функции для расчета элементов матриц
    double calculate_G(int i, int j) {
        // Пример: интеграл от произведения градиентов базисных функций
        double lambda = 1.0; // Коэффициент теплопроводности

        // Градиенты базисных функций (пример для биквадратичных функций)
        double grad_phi_i_x = (i % 3 == 0) ? -1.0 : ((i % 3 == 1) ? 0.0 : 1.0);
        double grad_phi_i_y = (i / 3 == 0) ? -1.0 : ((i / 3 == 1) ? 0.0 : 1.0);
        double grad_phi_j_x = (j % 3 == 0) ? -1.0 : ((j % 3 == 1) ? 0.0 : 1.0);
        double grad_phi_j_y = (j / 3 == 0) ? -1.0 : ((j / 3 == 1) ? 0.0 : 1.0);

        // Вычисляем скалярное произведение градиентов
        double grad_dot = grad_phi_i_x * grad_phi_j_x + grad_phi_i_y * grad_phi_j_y;

        // Интеграл по элементу (для прямоугольника площадь = hx * hy)
        double hx = 1.0; // Шаг по x (зависит от сетки)
        double hy = 1.0; // Шаг по y (зависит от сетки)
        double area = hx * hy;

        return lambda * grad_dot * area;
    }

    double calculate_C(int i, int j) {
        // Пример: интеграл от произведения базисных функций
        // Для простоты используем константу gamma = 1.0
        double gamma = 1.0; // Коэффициент теплоотдачи

        // Значения базисных функций (пример для биквадратичных функций)
        double phi_i = (i % 3 == 0) ? 1.0 : ((i % 3 == 1) ? 0.0 : 1.0);
        double phi_j = (j % 3 == 0) ? 1.0 : ((j % 3 == 1) ? 0.0 : 1.0);

        // Интеграл по элементу (для прямоугольника площадь = hx * hy)
        double hx = 1.0; // Шаг по x (зависит от сетки)
        double hy = 1.0; // Шаг по y (зависит от сетки)
        double area = hx * hy;

        return gamma * phi_i * phi_j * area;
    }

    double calculate_Ck(int i, int j) {
        // Пример: интеграл от произведения базисных функций с учетом коэффициентов
        // Для простоты используем константу beta = 1.0
        double beta = 1.0; // Коэффициент теплообмена

        // Значения базисных функций (пример для биквадратичных функций)
        double phi_i = (i % 3 == 0) ? 1.0 : ((i % 3 == 1) ? 0.0 : 1.0);
        double phi_j = (j % 3 == 0) ? 1.0 : ((j % 3 == 1) ? 0.0 : 1.0);

        // Интеграл по элементу (для прямоугольника площадь = hx * hy)
        double hx = 1.0; // Шаг по x (зависит от сетки)
        double hy = 1.0; // Шаг по y (зависит от сетки)
        double area = hx * hy;

        return beta * phi_i * phi_j * area;
    }
};

// Функции параметров уравнения
static double lambda_func(double x, double y) {
    return 1.0; // Пример: коэффициент теплопроводности
}

static double gamma_func(double x, double y) {
    return 0.0; // Пример: коэффициент теплоотдачи
}

static double f_func(double x, double y) {
    return -6.0 * x - 6.0 * y; // Пример: правая часть уравнения
}

// Краевые условия 1-го рода (Дирихле)
static double ug(double x, double y, int index) {
    return std::pow(x, 3) + std::pow(y, 3); // Пример: граничное условие
}

// Краевые условия 2-го рода (Неймана)
static double teta(double x, double y, int index) {
    return 0.0; // Пример: граничное условие
}

// Краевые условия 3-го рода (Робин)
static double beta_func(double x, double y, int index) {
    return 1.0; // Пример: коэффициент теплообмена
}

static double u_beta_func(double x, double y, int index) {
    return 0.0; // Пример: граничное условие
}

// Класс для управления сеткой
class Grid {
public:
    Mesh mesh;

    void loadnet(const std::string& filename) {
        std::ifstream fp(filename);
        if (!fp.is_open()) {
            throw std::runtime_error("Ошибка открытия файла " + filename);
        }

        mesh = Mesh(); // Сбрасываем mesh

        // Чтение xwn и Xw
        fp >> mesh.xwn;
        mesh.Xw.resize(mesh.xwn);
        for (auto& x : mesh.Xw) {
            fp >> x;
        }

        // Чтение ywn и Yw
        fp >> mesh.ywn;
        mesh.Yw.resize(mesh.ywn);
        for (auto& y : mesh.Yw) {
            fp >> y;
        }

        // Чтение L и W
        fp >> mesh.L;
        mesh.W.resize(4, std::vector<int>(mesh.L));
        for (int i = 0; i < mesh.L; ++i) {
            fp >> mesh.W[0][i] >> mesh.W[1][i] >> mesh.W[2][i] >> mesh.W[3][i];
            // Корректировка индексации (0-индексация)
            for (int p = 0; p < 4; ++p) {
                if (mesh.W[p][i] > 0) { // Если индекс начинается с 1
                    mesh.W[p][i]--;     // Переводим в 0-индексацию
                }
            }
        }
        fp.close();

        // Построение сетки
        mesh.n = (3 * mesh.xwn - 2) * (3 * mesh.ywn - 2);
        mesh.fict.resize(mesh.n, true); // Инициализация вектора fict

        // Проверка размера fict
        if (mesh.fict.size() != mesh.n) {
            throw std::runtime_error("Некорректный размер вектора fict.");
        }

        // Выделение памяти для координат узлов
        mesh.xy_x.resize(mesh.n, 0.0);
        mesh.xy_y.resize(mesh.n, 0.0);

        // Формирование координат узлов
        int z_idx = 0;
        for (int i = 0; i < mesh.ywn; ++i) {
            for (int j = 0; j < mesh.xwn; ++j) {
                double hx = (j < mesh.xwn - 1) ? (mesh.Xw[j + 1] - mesh.Xw[j]) : (mesh.Xw[j] - mesh.Xw[j - 1]);
                double hy = (i < mesh.ywn - 1) ? (mesh.Yw[i + 1] - mesh.Yw[i]) : (mesh.Yw[i] - mesh.Yw[i - 1]);

                // Центральные узлы
                for (int mi = 0; mi < 3; ++mi) {
                    for (int mj = 0; mj < 3; ++mj) {
                        if (z_idx >= mesh.n) break;
                        mesh.xy_x[z_idx] = mesh.Xw[j] + (mj / 2.0) * hx;
                        mesh.xy_y[z_idx] = mesh.Yw[i] + (mi / 2.0) * hy;
                        z_idx++;
                    }
                }
            }
        }

        // Вычисление количества конечных элементов
        mesh.k = 0;
        for (int p = 0; p < mesh.L; ++p) {
            int x_start = mesh.W[0][p];
            int x_end = mesh.W[1][p];
            int y_start = mesh.W[2][p];
            int y_end = mesh.W[3][p];

            // Проверка корректности индексов
            if (x_start >= x_end || y_start >= y_end) {
                throw std::runtime_error("Некорректные индексы подобласти.");
            }

            // Количество конечных элементов в подобласти
            int p_count = (x_end - x_start) * (y_end - y_start); // Исправлено
            mesh.k += p_count;
        }

        if (mesh.k <= 0) {
            throw std::runtime_error("Количество конечных элементов (k) должно быть больше 0. Проверьте входные данные.");
        }

        // Выделение памяти для nvkat и nvtr
        mesh.nvkat.resize(mesh.k, 0);
        mesh.nvtr.resize(mesh.k * NUM_NODES_PER_ELEMENT, -1);

        // Заполнение массивов nvkat и nvtr
        int elem_idx = 0;
        for (int p = 0; p < mesh.L; ++p) {
            for (int i = mesh.W[2][p]; i < mesh.W[3][p]; ++i) { // По строкам
                for (int j = mesh.W[0][p]; j < mesh.W[1][p]; ++j) { // По столбцам
                    if (elem_idx >= mesh.k) break;
                    mesh.nvkat[elem_idx] = p;
                    // Заполнение локальных номеров узлов
                    for (int m = 0; m < NUM_NODES_PER_ELEMENT; ++m) {
                        int row = i * 3 + (m / 3);
                        int col = j * 3 + (m % 3);

                        // Проверка на отрицательные индексы
                        if (row < 0 || col < 0) {
                            throw std::runtime_error("Некорректные индексы row или col: row=" + std::to_string(row) + ", col=" + std::to_string(col));
                        }

                        int global_index = row * (3 * mesh.xwn - 2) + col;

                        // Проверка на выход за границы
                        if (global_index < 0 || global_index >= mesh.n) {
                            throw std::runtime_error("Некорректный индекс узла: " + std::to_string(global_index));
                        }

                        if (row >= (3 * mesh.ywn - 2) || col >= (3 * mesh.xwn - 2)) {
                            mesh.nvtr[elem_idx * NUM_NODES_PER_ELEMENT + m] = -1; // Фиктивный узел
                            mesh.fict[global_index] = false;
                        }
                        else {
                            mesh.nvtr[elem_idx * NUM_NODES_PER_ELEMENT + m] = global_index;
                        }
                    }
                    elem_idx++;
                }
            }
        }

        std::cout << "Сетка загружена: n=" << mesh.n << ", k=" << mesh.k << std::endl;
        std::cout << "Проверка nvtr: первый элемент = " << mesh.nvtr[0] << ", последний элемент = " << mesh.nvtr.back() << std::endl;
    }
};

// Класс для управления краевыми условиями
class BoundaryConditionsLoader {
public:
    BoundaryCondition bc;

    void loadbc(const std::string& ku1, const std::string& ku2, const std::string& ku3) {
        load_boundary_condition(ku1, bc.kt1, bc.l1);
        load_boundary_condition(ku2, bc.nvk2, bc.nvr2);
        load_boundary_condition(ku3, bc.nvk3, bc.nvr3);
    }

private:
    void load_boundary_condition(const std::string& filename, std::vector<int>& indices, std::vector<std::pair<int, int>>& ranges) {
        std::ifstream fp(filename);
        if (!fp.is_open()) {
            throw std::runtime_error("Ошибка открытия файла " + filename);
        }

        int count;
        fp >> count;
        indices.resize(count);
        ranges.resize(count, { -1, -1 });

        for (int i = 0; i < count; ++i) {
            int index, start, end;
            fp >> index >> start >> end;
            indices[i] = index - 1; // Корректировка индексации
            ranges[i] = { start - 1, end - 1 };
        }

        fp.close();
    }
};

// Класс для управления разреженной матрицей и операциями над ней
class SparseMatrixCSR {
public:
    CSRMatrix A;

    SparseMatrixCSR(int size = 0) : A(size) {}

    void initialize(int size) {
        A.n = size;
        A.row_ptr.assign(size + 1, 0);
    }

    void allocate_structure(const Mesh& mesh) {
        std::vector<int> temp(mesh.n, 0);
        for (int ielem = 0; ielem < mesh.k; ielem++) {
            for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
                int node_i = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
                if (node_i == -1 || !mesh.fict[node_i]) continue;
                for (int j = 0; j < NUM_NODES_PER_ELEMENT; j++) {
                    int node_j = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + j];
                    if (node_j == -1 || !mesh.fict[node_j] || node_j > node_i) continue;
                    temp[node_i]++;
                }
            }
        }

        A.row_ptr[0] = 0;
        for (int i = 0; i < mesh.n; i++) {
            A.row_ptr[i + 1] = A.row_ptr[i] + temp[i];
        }
        A.nnz = A.row_ptr[mesh.n];

        A.col_idx.assign(A.nnz, 0);
        A.values.assign(A.nnz, 0.0);

        std::vector<int> current(mesh.n, 0);
        for (int ielem = 0; ielem < mesh.k; ielem++) {
            for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
                int row = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
                if (row == -1 || !mesh.fict[row]) continue;
                for (int j = 0; j < NUM_NODES_PER_ELEMENT; j++) {
                    int col = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + j];
                    if (col == -1 || !mesh.fict[col] || col > row) continue;
                    int pos = A.row_ptr[row] + current[row]++;
                    A.col_idx[pos] = col;
                }
            }
        }

        for (int i = 0; i < mesh.n; ++i) {
            int start = A.row_ptr[i];
            int end = A.row_ptr[i + 1];
            std::sort(A.col_idx.begin() + start, A.col_idx.begin() + end);
        }
    }

    void add_local_to_global(int ielem, const FiniteElement& fe, const Mesh& mesh, std::vector<double>& b) {
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
            int row = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (row == -1 || !mesh.fict[row]) continue;

            for (int j = 0; j < NUM_NODES_PER_ELEMENT; j++) {
                int col = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + j];
                if (col == -1 || !mesh.fict[col]) continue;

                int start = A.row_ptr[row];
                int end = A.row_ptr[row + 1];
                auto it = std::lower_bound(A.col_idx.begin() + start, A.col_idx.begin() + end, col);
                if (it != A.col_idx.begin() + end && *it == col) {
                    int pos = std::distance(A.col_idx.begin(), it);
                    A.values[pos] += fe.Ak[i * NUM_NODES_PER_ELEMENT + j];
                }
            }
        }

        for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
            int row = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (row == -1 || !mesh.fict[row]) continue;
            b[row] += fe.Bk[i];
        }
    }

    std::vector<double> multiply(const std::vector<double>& v) const {
        std::vector<double> res(A.n, 0.0);
        for (int i = 0; i < A.n; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                res[i] += A.values[j] * v[A.col_idx[j]];
            }
        }
        return res;
    }
};

// Класс для решения СЛАУ методом сопряжённых градиентов (CG)
class ConjugateGradientSolver {
public:
    void solve(const CSRMatrix& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
        int n = A.n;
        std::vector<double> r(n, 0.0);
        std::vector<double> p(n, 0.0);
        std::vector<double> Ap(n, 0.0);

        std::vector<double> Ax = multiply(A, x);
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - Ax[i];
        }

        p = r;

        double rsold = dot_product(r, r);

        for (int iter = 0; iter < max_iter; ++iter) {
            Ap = multiply(A, p);

            double pAp = dot_product(p, Ap);
            if (std::abs(pAp) < 1e-12) { // Проверка на нулевое значение
                throw std::runtime_error("pAp близко к нулю. Возможно, матрица A вырождена.");
            }
            double alpha = rsold / pAp;

            for (int i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
            }

            for (int i = 0; i < n; ++i) {
                r[i] -= alpha * Ap[i];
            }

            double rsnew = dot_product(r, r);
            if (std::sqrt(rsnew) < tol) {
                std::cout << "Сходимость достигнута за " << iter + 1 << " итераций." << std::endl;
                break;
            }

            double beta = rsnew / rsold;
            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + beta * p[i];
            }

            rsold = rsnew;

            if (iter == max_iter - 1) {
                std::cout << "Максимальное количество итераций (" << max_iter << ") достигнуто." << std::endl;
            }
        }
    }

private:
    std::vector<double> multiply(const CSRMatrix& A, const std::vector<double>& v) const {
        std::vector<double> res(A.n, 0.0);
        for (int i = 0; i < A.n; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                res[i] += A.values[j] * v[A.col_idx[j]];
            }
        }
        return res;
    }

    double dot_product(const std::vector<double>& a, const std::vector<double>& b) const {
        double res = 0.0;
        for (int i = 0; i < a.size(); ++i) {
            res += a[i] * b[i];
        }
        return res;
    }
};

// Класс для решения СЛАУ методом ЛОС с неполной факторизацией
class LOS_Solver {
public:
    void solve(const CSRMatrix& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
        int n = A.n;
        std::vector<double> r(n, 0.0);
        std::vector<double> z(n, 0.0);
        std::vector<double> p(n, 0.0);
        std::vector<double> Ap(n, 0.0);

        // Начальное приближение
        std::vector<double> Ax = multiply(A, x);
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - Ax[i];
        }

        z = r;
        p = multiply(A, z);

        double alpha, beta, rr;

        for (int iter = 0; iter < max_iter; ++iter) {
            rr = dot_product(r, r);
            alpha = rr / dot_product(p, p);

            for (int i = 0; i < n; ++i) {
                x[i] += alpha * z[i];
                r[i] -= alpha * p[i];
            }

            double rr_new = dot_product(r, r);
            if (std::sqrt(rr_new) < tol) {
                std::cout << "Сходимость достигнута за " << iter + 1 << " итераций." << std::endl;
                break;
            }

            beta = rr_new / rr;
            for (int i = 0; i < n; ++i) {
                z[i] = r[i] + beta * z[i];
            }

            p = multiply(A, z);

            if (iter == max_iter - 1) {
                std::cout << "Максимальное количество итераций (" << max_iter << ") достигнуто." << std::endl;
            }
        }
    }

private:
    std::vector<double> multiply(const CSRMatrix& A, const std::vector<double>& v) const {
        std::vector<double> res(A.n, 0.0);
        for (int i = 0; i < A.n; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                res[i] += A.values[j] * v[A.col_idx[j]];
            }
        }
        return res;
    }

    double dot_product(const std::vector<double>& a, const std::vector<double>& b) const {
        double res = 0.0;
        for (int i = 0; i < a.size(); ++i) {
            res += a[i] * b[i];
        }
        return res;
    }
};

// Основной класс, объединяющий все компоненты
class FiniteElementSolver {
public:
    Mesh mesh;
    BoundaryCondition bc;
    CSRMatrix A;
    std::vector<double> b_vector;
    std::vector<double> solution;
    FiniteElement fe;

    FiniteElementSolver() : A(0) {}

    void initialize() {
        fe.initialize_matrices();
    }

    void load_data(const std::string& mesh_file, const std::string& ku1, const std::string& ku2, const std::string& ku3) {
        Grid grid;
        grid.loadnet(mesh_file);

        // Перемещаем сетку из grid в mesh
        mesh = std::move(grid.mesh);

        // Проверка количества конечных элементов
        if (mesh.k <= 0) {
            throw std::runtime_error("Количество конечных элементов (k) должно быть больше 0. Проверьте входные данные.");
        }

        // Загрузка краевых условий
        BoundaryConditionsLoader loader;
        loader.loadbc(ku1, ku2, ku3);
        bc = loader.bc;

        std::cout << "Данные загружены: n=" << mesh.n << ", k=" << mesh.k << std::endl;
        std::cout << "Проверка краевых условий: количество условий Дирихле = " << bc.kt1.size() << std::endl;
    }

    void assemble_system() {
        std::cout << "Начало сборки системы." << std::endl;

        SparseMatrixCSR matrix(mesh.n);
        matrix.allocate_structure(mesh);

        b_vector.assign(mesh.n, 0.0);

        for (int ielem = 0; ielem < mesh.k; ++ielem) {
            compile_local_element(ielem);
            matrix.add_local_to_global(ielem, fe, mesh, b_vector);
        }

        apply_boundary_conditions(matrix);

        A = matrix.A;

        std::cout << "Сборка системы завершена." << std::endl;

        // Отладочный вывод матрицы A
        std::cout << "Матрица A:" << std::endl;
        for (int i = 0; i < A.n; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                std::cout << "A[" << i << ", " << A.col_idx[j] << "] = " << A.values[j] << std::endl;
            }
        }

        // Отладочный вывод вектора b
        std::cout << "Вектор b:" << std::endl;
        for (int i = 0; i < b_vector.size(); ++i) {
            std::cout << "b[" << i << "] = " << b_vector[i] << std::endl;
        }

    }

    void solve_system() {
        solution.assign(mesh.n, 0.0);
        ConjugateGradientSolver solver;
        solver.solve(A, b_vector, solution, MAX_ITER, EPSILON);
    }

    void save_solution(const std::string& filename) const {
        std::ofstream fp(filename);
        if (!fp.is_open()) {
            throw std::runtime_error("Ошибка открытия файла " + filename + " для записи.");
        }
        for (const auto& val : solution) {
            fp << val << "\t";
        }
        fp.close();
        std::cout << "Решение записано в файл " << filename << std::endl;
    }

private:
    void compile_local_element(int ielem) {
        std::vector<double> element_x(NUM_NODES_PER_ELEMENT, 0.0);
        std::vector<double> element_y(NUM_NODES_PER_ELEMENT, 0.0);
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node == -1 || !mesh.fict[node]) continue;
            if (node < 0 || node >= mesh.n) {
                throw std::runtime_error("Некорректный индекс узла в compile_local_element.");
            }
            element_x[i] = mesh.xy_x[node];
            element_y[i] = mesh.xy_y[node];
        }

        double lambdak = lambda_func(element_x[8], element_y[8]);
        double gammak = 0.0;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            if (mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i] == -1) continue;
            gammak += gamma_func(element_x[i], element_y[i]);
        }
        gammak /= NUM_NODES_PER_ELEMENT;

        double hx = element_x[2] - element_x[0];
        double hy = element_y[6] - element_y[0];
        double k1 = lambdak * hy / hx;
        double k2 = lambdak * hx / hy;
        double k3 = gammak * hx * hy;

        for (int i = 0; i < NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT; ++i) {
            fe.Ak[i] = k1 * fe.G[i] + k2 * fe.C[i] + k3 * fe.Ck[i];
        }

        std::vector<double> fk(NUM_NODES_PER_ELEMENT, 0.0);
        double k_area = hx * hy;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node == -1) continue;
            fk[i] = f_func(element_x[i], element_y[i]);
        }

        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            fe.Bk[i] = fk[i] * k_area;
        }
        // Отладочный вывод локальных матриц
        std::cout << "Локальная матрица Ak для элемента " << ielem << ":" << std::endl;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
                std::cout << fe.Ak[i * NUM_NODES_PER_ELEMENT + j] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Локальный вектор Bk для элемента " << ielem << ":" << std::endl;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            std::cout << fe.Bk[i] << " ";
        }
        std::cout << std::endl;
    }

    void apply_boundary_conditions(SparseMatrixCSR& matrix) {
        for (size_t i = 0; i < bc.kt1.size(); ++i) {
            int ind = bc.l1[i].first;
            if (ind < 0 || ind >= mesh.n) {
                throw std::runtime_error("Некорректный индекс узла для краевого условия Дирихле.");
            }

            double u_val = ug(mesh.xy_x[ind], mesh.xy_y[ind], bc.kt1[i]);

            // Обнуляем строку и столбец для узла с условием Дирихле
            for (int j = matrix.A.row_ptr[ind]; j < matrix.A.row_ptr[ind + 1]; ++j) {
                matrix.A.values[j] = 0.0;
            }

            // Устанавливаем диагональный элемент в 1
            bool diag_found = false;
            for (int j = matrix.A.row_ptr[ind]; j < matrix.A.row_ptr[ind + 1]; ++j) {
                if (matrix.A.col_idx[j] == ind) {
                    matrix.A.values[j] = 1.0;
                    diag_found = true;
                    break;
                }
            }
            if (!diag_found) {
                throw std::runtime_error("Диагональный элемент не найден.");
            }

            // Устанавливаем значение в векторе b
            b_vector[ind] = u_val;
        }
    }
};

int main() {
    try {
        std::setlocale(LC_ALL, "Russian");

        FiniteElementSolver solver;
        solver.initialize();

        solver.load_data("st.txt", "ku1.txt", "ku2.txt", "ku3.txt");

        solver.assemble_system();
        solver.solve_system();
        solver.save_solution("q.txt");

        std::cout << "Решение записано в файл q.txt" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return 0;
}

