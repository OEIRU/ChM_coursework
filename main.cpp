#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <locale>
#include <string>
#include <algorithm> // Для std::sort и std::lower_bound

// Константы
constexpr int NUM_NODES_PER_ELEMENT = 16;
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
    int xwn, ywn; // Количество узлов по X и Y
    std::vector<double> Xw; // Координаты узлов по X
    std::vector<double> Yw; // Координаты узлов по Y
    int L; // Число подобластей
    std::vector<std::vector<int>> W; // Координаты подобластей (левый, правый, нижний, верхний)

    // Глобальные номера узлов для каждого КЭ
    std::vector<int> nvkat; // Номера подобластей для каждого КЭ
    std::vector<int> nvtr; // Глобальные номера узлов для каждого КЭ (одномерный массив)
    std::vector<bool> fict; // Флаг фиктивных узлов
    int n; // Общее количество узлов
    int k; // Количество конечных элементов

    // Координаты узлов
    std::vector<double> xy_x;
    std::vector<double> xy_y;
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
    std::vector<double> G;
    std::vector<double> C;
    std::vector<double> Ck;
    std::vector<double> Ak;
    std::vector<double> Bk;

    FiniteElement()
        : G(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        C(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        Ck(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        Ak(NUM_NODES_PER_ELEMENT* NUM_NODES_PER_ELEMENT, 0.0),
        Bk(NUM_NODES_PER_ELEMENT, 0.0) {
    }

    void initialize_matrices() {
        // Для простоты заполняем G, C и Ck как единичные матрицы
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            G[i * NUM_NODES_PER_ELEMENT + i] = 1.0;
            C[i * NUM_NODES_PER_ELEMENT + i] = 1.0;
            Ck[i * NUM_NODES_PER_ELEMENT + i] = 1.0;
        }
    }
};

// Функции параметров уравнения
static double lambda_func(double x, double y) {
    return 1.0;
}

static double gamma_func(double x, double y) {
    // Разложение по билинейным базисным функциям
    return 0.0;
}

static double f_func(double x, double y) {
    // Тест 1.1: u = x^3 + y^3 => f = -6x -6y
    return -6.0 * x - 6.0 * y;
}

// Краевые условия 1-го рода (Дирихле)
static double ug(double x, double y, int index) {
    return std::pow(x, 3) + std::pow(y, 3);
}

// Краевые условия 2-го рода (Неймана)
static double teta(double x, double y, int index) {
    return 0.0;
}

// Краевые условия 3-го рода (Робин)
static double beta_func(double x, double y, int index) {
    return 1.0;
}

static double u_beta_func(double x, double y, int index) {
    return 0.0;
}

// Класс для управления сеткой
class Grid {
public:
    Mesh mesh;

    void loadnet(const std::string& filename) {
        std::ifstream fp(filename);
        if (!fp.is_open()) {
            std::cerr << "Ошибка открытия файла " << filename << std::endl;
            exit(EXIT_FAILURE);
        }

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
                mesh.W[p][i]--;
            }
        }
        fp.close();

        // Построение сетки
        mesh.n = (3 * mesh.xwn - 2) * (3 * mesh.ywn - 2);
        mesh.fict.resize(mesh.n, true);

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
            int p_count = (mesh.W[1][p] - mesh.W[0][p] - 1) * (mesh.W[3][p] - mesh.W[2][p] - 1);
            mesh.k += p_count;
        }

        // Выделение памяти для nvkat и nvtr
        mesh.nvkat.resize(mesh.k, 0);
        mesh.nvtr.resize(mesh.k * NUM_NODES_PER_ELEMENT, -1);

        // Заполнение массивов nvkat и nvtr
        int elem_idx = 0;
        for (int p = 0; p < mesh.L; ++p) {
            for (int i = mesh.W[2][p]; i < mesh.W[3][p] - 1; ++i) { // По строкам
                for (int j = mesh.W[0][p]; j < mesh.W[1][p] - 1; ++j) { // По столбцам
                    if (elem_idx >= mesh.k) break;
                    mesh.nvkat[elem_idx] = p;
                    // Заполнение локальных номеров узлов (пример)
                    for (int m = 0; m < NUM_NODES_PER_ELEMENT; ++m) {
                        int row = i * 3 + (m / 4);
                        int col = j * 3 + (m % 4);
                        if (row >= (3 * mesh.ywn - 2) || col >= (3 * mesh.xwn - 2)) {
                            mesh.nvtr[elem_idx * NUM_NODES_PER_ELEMENT + m] = -1; // Фиктивный узел
                            mesh.fict[row * (3 * mesh.xwn - 2) + col] = false;
                        }
                        else {
                            mesh.nvtr[elem_idx * NUM_NODES_PER_ELEMENT + m] = row * (3 * mesh.xwn - 2) + col;
                        }
                    }
                    elem_idx++;
                }
            }
        }

        // Инициализация остальных массивов уже выполнена
        std::cout << "Сетка загружена: n=" << mesh.n << ", k=" << mesh.k << std::endl;
    }
};

// Класс для управления краевыми условиями
class BoundaryConditionsLoader {
public:
    BoundaryCondition bc;

    void loadbc(const std::string& ku1, const std::string& ku2, const std::string& ku3) {
        // Загрузка краевых условий 1-го рода
        load_boundary_condition(ku1, bc.kt1, bc.l1, "ку1.txt");

        // Загрузка краевых условий 2-го рода
        load_boundary_condition(ku2, bc.nvk2, bc.nvr2, "ку2.txt");

        // Загрузка краевых условий 3-го рода
        load_boundary_condition(ku3, bc.nvk3, bc.nvr3, "ку3.txt");
    }

private:
    void load_boundary_condition(const std::string& filename, std::vector<int>& indices, std::vector<std::pair<int, int>>& ranges, const std::string& condition_name = "") {
        std::ifstream fp(filename);
        if (!fp.is_open()) {
            std::cerr << "Ошибка открытия файла " << filename << std::endl;
            exit(EXIT_FAILURE);
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
        // Подсчёт ненулевых элементов в каждой строке
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

        // Формирование row_ptr
        A.row_ptr[0] = 0;
        for (int i = 0; i < mesh.n; i++) {
            A.row_ptr[i + 1] = A.row_ptr[i] + temp[i];
        }
        A.nnz = A.row_ptr[mesh.n];

        // Выделение памяти для col_idx и values
        A.col_idx.assign(A.nnz, 0);
        A.values.assign(A.nnz, 0.0);

        // Временный массив для текущей позиции вставки в каждой строке
        std::vector<int> current(mesh.n, 0);

        // Заполнение col_idx
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

        // Сортировка столбцов внутри каждой строки для ускорения доступа
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
                if (col == -1 || !mesh.fict[col] || col > row) continue;

                // Поиск позиции (row, col) в CSR формате
                int start = A.row_ptr[row];
                int end = A.row_ptr[row + 1];
                // Используем бинарный поиск, так как col_idx отсортирован
                auto it = std::lower_bound(A.col_idx.begin() + start, A.col_idx.begin() + end, col);
                if (it != A.col_idx.begin() + end && *it == col) {
                    int pos = std::distance(A.col_idx.begin(), it);
                    A.values[pos] += fe.Ak[i * NUM_NODES_PER_ELEMENT + j];
                }
                else {
                    std::cerr << "Ошибка: Элемент (" << row << ", " << col << ") не найден в матрице." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }

        // Добавление локального вектора Bk в глобальный вектор b
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
            int row = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (row == -1 || !mesh.fict[row]) continue;
            b[row] += fe.Bk[i];
        }
    }

    // Функция для умножения матрицы на вектор
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

        // r = b - A*x
        std::vector<double> Ax = multiply(A, x);
        for (int i = 0; i < n; ++i) {
            r[i] = b[i] - Ax[i];
        }

        // p = r
        p = r;

        double rsold = dot_product(r, r);

        for (int iter = 0; iter < max_iter; ++iter) {
            // Ap = A * p
            Ap = multiply(A, p);

            double pAp = dot_product(p, Ap);
            if (pAp == 0.0) {
                std::cerr << "Ошибка: Деление на ноль при вычислении alpha." << std::endl;
                exit(EXIT_FAILURE);
            }
            double alpha = rsold / pAp;

            // x = x + alpha * p
            for (int i = 0; i < n; ++i) {
                x[i] += alpha * p[i];
            }

            // r = r - alpha * Ap
            for (int i = 0; i < n; ++i) {
                r[i] -= alpha * Ap[i];
            }

            double rsnew = dot_product(r, r);
            if (std::sqrt(rsnew) < tol) {
                std::cout << "Сходимость достигнута за " << iter + 1 << " итераций." << std::endl;
                break;
            }

            // p = r + (rsnew / rsold) * p
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
    // Вспомогательные функции
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
        mesh = grid.mesh;

        // Проверка количества конечных элементов
        if (mesh.k <= 0) {
            std::cerr << "Ошибка: Количество конечных элементов (k) должно быть больше 0. Проверьте входные данные." << std::endl;
            exit(EXIT_FAILURE);
        }

        load_boundary_conditions(ku1, ku2, ku3);
    }

    void assemble_system() {
        // Инициализация матрицы
        SparseMatrixCSR matrix(mesh.n);

        // Подсчёт ненулевых элементов и выделение структуры
        matrix.allocate_structure(mesh);

        // Инициализация вектора правой части
        b_vector.assign(mesh.n, 0.0);

        // Сборка матрицы и вектора правой части
        for (int ielem = 0; ielem < mesh.k; ++ielem) {
            compile_local_element(ielem);
            matrix.add_local_to_global(ielem, fe, mesh, b_vector);
        }

        // Применение краевых условий
        apply_boundary_conditions(matrix);

        // Перенос собранной матрицы
        A = matrix.A;
    }

    void solve_system() {
        solution.assign(mesh.n, 0.0);
        ConjugateGradientSolver solver;
        solver.solve(A, b_vector, solution, MAX_ITER, EPSILON);
    }

    void save_solution(const std::string& filename) const {
        std::ofstream fp(filename);
        if (!fp.is_open()) {
            std::cerr << "Ошибка открытия файла " << filename << " для записи." << std::endl;
            exit(EXIT_FAILURE);
        }
        for (const auto& val : solution) {
            fp << val << "\t";
        }
        fp.close();
        std::cout << "Решение записано в файл " << filename << std::endl;
    }

private:
    void load_boundary_conditions(const std::string& ku1, const std::string& ku2, const std::string& ku3) {
        BoundaryConditionsLoader loader;
        loader.loadbc(ku1, ku2, ku3);
        bc = loader.bc;
    }

    void compile_local_element(int ielem) {
        // Получение координат узлов элемента
        std::vector<double> element_x(NUM_NODES_PER_ELEMENT, 0.0);
        std::vector<double> element_y(NUM_NODES_PER_ELEMENT, 0.0);
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node == -1) continue; // Фиктивный узел
            element_x[i] = mesh.xy_x[node];
            element_y[i] = mesh.xy_y[node];
        }

        // Формирование локальной матрицы жесткости Ak
        double lambdak = lambda_func(element_x[8], element_y[8]);
        double gammak = 0.0;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            if (mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i] == -1) continue;
            gammak += gamma_func(element_x[i], element_y[i]);
        }
        gammak /= NUM_NODES_PER_ELEMENT;

        double hx = element_x[15] - element_x[0];
        double hy = element_y[15] - element_y[0];
        double k1 = lambdak * hy / hx;
        double k2 = lambdak * hx / hy;
        double k3 = gammak * hx * hy;

        for (int i = 0; i < NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT; ++i) {
            fe.Ak[i] = k1 * fe.G[i] + k2 * fe.C[i] + k3 * fe.Ck[i];
        }

        // Формирование локального вектора правой части Bk
        std::vector<double> fk(NUM_NODES_PER_ELEMENT, 0.0);
        double k_area = hx * hy;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node == -1) continue;
            fk[i] = f_func(element_x[i], element_y[i]);
        }

        for (int m = 0; m < NUM_NODES_PER_ELEMENT; ++m) {
            fe.Bk[m] = 0.0;
            for (int n = 0; n < NUM_NODES_PER_ELEMENT; ++n) {
                fe.Bk[m] += fe.Ck[m * NUM_NODES_PER_ELEMENT + n] * fk[n];
            }
            fe.Bk[m] *= k_area;
        }
    }

    void apply_boundary_conditions(SparseMatrixCSR& matrix) {
        // Учет краевых условий 1-го рода (Дирихле)
        for (size_t i = 0; i < bc.kt1.size(); ++i) {
            int ind = bc.l1[i].first;
            // Установка значения вектора решения
            double u_val = ug(mesh.xy_x[ind], mesh.xy_y[ind], bc.kt1[i]);
            solution[ind] = u_val;

            // Обнуление строки матрицы и установка диагонального элемента
            int row_start = matrix.A.row_ptr[ind];
            int row_end = matrix.A.row_ptr[ind + 1];
            for (int j = row_start; j < row_end; ++j) {
                matrix.A.values[j] = 0.0;
            }

            // Поиск диагонального элемента
            bool diag_found = false;
            for (int j = row_start; j < row_end; ++j) {
                if (matrix.A.col_idx[j] == ind) {
                    matrix.A.values[j] = 1.0;
                    diag_found = true;
                    break;
                }
            }
            if (!diag_found) {
                std::cerr << "Ошибка: Диагональный элемент (" << ind << ", " << ind << ") не найден." << std::endl;
                exit(EXIT_FAILURE);
            }

            // Установка соответствующего значения вектора правой части
            b_vector[ind] = u_val;
        }

        // Учет краевых условий 2-го рода (Неймана)
        for (size_t i = 0; i < bc.nvk2.size(); ++i) {
            int ind_start = bc.nvr2[i].first;
            int ind_end = bc.nvr2[i].second;
            double theta = teta(mesh.xy_x[ind_start], mesh.xy_y[ind_start], bc.nvk2[i]);

            for (int ind = ind_start; ind <= ind_end; ++ind) {
                b_vector[ind] += theta;
            }
        }

        // Учет краевых условий 3-го рода (Робин)
        for (size_t i = 0; i < bc.nvk3.size(); ++i) {
            int ind_start = bc.nvr3[i].first;
            int ind_end = bc.nvr3[i].second;
            double beta_val = beta_func(mesh.xy_x[ind_start], mesh.xy_y[ind_start], bc.nvk3[i]);
            double u_beta_val = u_beta_func(mesh.xy_x[ind_start], mesh.xy_y[ind_start], bc.nvk3[i]);

            for (int ind = ind_start; ind <= ind_end; ++ind) {
                // Добавление вклада в диагональные элементы матрицы
                int row_start = matrix.A.row_ptr[ind];
                int row_end = matrix.A.row_ptr[ind + 1];
                bool diag_found = false;
                for (int j = row_start; j < row_end; ++j) {
                    if (matrix.A.col_idx[j] == ind) {
                        matrix.A.values[j] += beta_val;
                        diag_found = true;
                        break;
                    }
                }
                if (!diag_found) {
                    std::cerr << "Ошибка: Диагональный элемент (" << ind << ", " << ind << ") не найден для Робин условий." << std::endl;
                    exit(EXIT_FAILURE);
                }
                // Добавление вклада вектора правой части
                b_vector[ind] += beta_val * u_beta_val;
            }
        }
    }
};

int main() {
    // Установка локали для корректного отображения русских символов
    std::setlocale(LC_ALL, "Russian");

    // Создание экземпляра решателя
    FiniteElementSolver solver;
    solver.initialize();

    // Загрузка данных
    solver.load_data("st.txt", "ku1.txt", "ku2.txt", "ku3.txt");

    // Проверка количества конечных элементов
    if (solver.mesh.k == 0) {
        std::cerr << "Ошибка: Количество конечных элементов (k) равно 0. Проверьте входные данные." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Сборка системы
    solver.assemble_system();

    // Решение системы
    solver.solve_system();

    // Сохранение решения
    solver.save_solution("q.txt");

    std::cout << "Решение записано в файл q.txt" << std::endl;

    // Открытие файла с решением (для Windows)
    // Важно: Команда system("notepad q.txt") зависит от операционной системы и может быть небезопасной.
    // Рекомендуется использовать более безопасные методы для открытия файлов.
    system("notepad q.txt");

    return 0;
}
