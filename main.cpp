#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <locale>
#include <unordered_map>

constexpr int NUM_NODES_PER_ELEMENT = 9;
constexpr int MAX_ITER = 100000;
constexpr double EPSILON = 1e-8;  // Увеличено для улучшения устойчивости

struct CSRMatrix {
    int n, nnz;
    std::vector<double> values;
    std::vector<int> col_idx, row_ptr;

    CSRMatrix(int size = 0) : n(size), nnz(0), row_ptr(size + 1, 0) {}
};

struct Mesh {
    int xwn, ywn, L, n, k;
    std::vector<double> Xw, Yw, xy_x, xy_y;
    std::vector<std::vector<int>> W;
    std::vector<int> nvkat, nvtr;
    std::vector<bool> fict;
};

struct BoundaryCondition {
    std::vector<int> kt1, nvk2, nvk3;
    std::vector<std::pair<int, int>> l1, nvr2, nvr3;
    std::vector<double> kt1_values, nvk2_values, nvk3_values;
};

double lambda_func(double x, double y) { return 1.0; }
double gamma_func(double x, double y) { return 0.0; }
double f_func(double x, double y) { return -6.0 * x - 6.0 * y; }
double ug(double x, double y, int index) { return std::pow(x, 3) + std::pow(y, 3); }
double teta(double x, double y, int index) { return 0.0; }
double beta_func(double x, double y, int index) { return 1.0; }
double u_beta_func(double x, double y, int index) { return 0.0; }

double bilinear_basis_function(double xi, double eta, int i) {
    switch (i) {
        case 0: return 0.25 * (1 - xi) * (1 - eta);
        case 1: return 0.25 * (1 + xi) * (1 - eta);
        case 2: return 0.25 * (1 + xi) * (1 + eta);
        case 3: return 0.25 * (1 - xi) * (1 + eta);
        default: throw std::runtime_error("Invalid basis function index.");
    }
}

double calculate_G(int i, int j, double x, double y) {
    double grad_phi_i_x = (i % 3 == 0) ? -1.0 : ((i % 3 == 1) ? 0.0 : 1.0);
    double grad_phi_i_y = (i / 3 == 0) ? -1.0 : ((i / 3 == 1) ? 0.0 : 1.0);
    double grad_phi_j_x = (j % 3 == 0) ? -1.0 : ((j % 3 == 1) ? 0.0 : 1.0);
    double grad_phi_j_y = (j / 3 == 0) ? -1.0 : ((j / 3 == 1) ? 0.0 : 1.0);
    return lambda_func(x, y) * (grad_phi_i_x * grad_phi_j_x + grad_phi_i_y * grad_phi_j_y);
}

double calculate_C(int i, int j, double x, double y) {
    double phi_i = (i % 3 == 0) ? 1.0 : ((i % 3 == 1) ? 0.0 : 1.0);
    double phi_j = (j % 3 == 0) ? 1.0 : ((j % 3 == 1) ? 0.0 : 1.0);
    return gamma_func(x, y) * phi_i * phi_j;
}

double calculate_Ck(int i, int j, double x, double y) {
    double phi_i = (i % 3 == 0) ? 1.0 : ((i % 3 == 1) ? 0.0 : 1.0);
    double phi_j = (j % 3 == 0) ? 1.0 : ((j % 3 == 1) ? 0.0 : 1.0);
    return beta_func(x, y, 0) * phi_i * phi_j;
}


void loadnet(const std::string& filename, Mesh& mesh) {
    std::ifstream fp(filename);
    if (!fp) throw std::runtime_error("Ошибка открытия файла " + filename);

    // Чтение количества узлов по x и y
    fp >> mesh.xwn;
    fp >> mesh.ywn;

    std::cout << "Количество узлов по x: " << mesh.xwn << "\n";
    std::cout << "Количество узлов по y: " << mesh.ywn << "\n";

    // Чтение координат узлов по x
    mesh.Xw.resize(mesh.xwn);
    for (int i = 0; i < mesh.xwn; ++i) {
        fp >> mesh.Xw[i];
        std::cout << "Xw[" << i << "]: " << mesh.Xw[i] << "\n";
    }

    // Чтение координат узлов по y
    mesh.Yw.resize(mesh.ywn);
    for (int i = 0; i < mesh.ywn; ++i) {
        fp >> mesh.Yw[i];
        std::cout << "Yw[" << i << "]: " << mesh.Yw[i] << "\n";
    }

    // Чтение количества элементов
    fp >> mesh.L;
    std::cout << "Количество элементов: " << mesh.L << "\n";

    // Чтение индексов угловых узлов элементов
    mesh.W.resize(4, std::vector<int>(mesh.L));
    for (int i = 0; i < mesh.L; ++i) {
        for (int p = 0; p < 4; ++p) {
            fp >> mesh.W[p][i];
            if (mesh.W[p][i] > 0) mesh.W[p][i]--;
            std::cout << "W[" << p << "][" << i << "]: " << mesh.W[p][i] << "\n";
        }
    }

    // Количество узлов равно mesh.xwn * mesh.ywn
    mesh.n = mesh.xwn * mesh.ywn;
    std::cout << "Количество узлов (mesh.n): " << mesh.n << "\n";  // Отладочное сообщение

    if (mesh.n <= 0) throw std::runtime_error("Некорректное количество узлов.");
    mesh.fict.assign(mesh.n, true);
    mesh.xy_x.resize(mesh.n);
    mesh.xy_y.resize(mesh.n);

    // Вычисляем шаг сетки
    double hx = (mesh.xwn > 1) ? (mesh.Xw[1] - mesh.Xw[0]) : 1.0;
    double hy = (mesh.ywn > 1) ? (mesh.Yw[1] - mesh.Yw[0]) : 1.0;

    std::cout << "hx: " << hx << ", hy: " << hy << "\n";  // Отладочное сообщение

    if (hx <= 0 || hy <= 0) throw std::runtime_error("Нулевой или отрицательный шаг сетки.");

    int z_idx = 0;
    for (int i = 0; i < mesh.ywn; ++i) {
        for (int j = 0; j < mesh.xwn; ++j) {
            // Вычисляем координаты узла
            mesh.xy_x[z_idx] = mesh.Xw[j];
            mesh.xy_y[z_idx] = mesh.Yw[i];
            std::cout << "xy_x[" << z_idx << "]: " << mesh.xy_x[z_idx] << ", xy_y[" << z_idx << "]: " << mesh.xy_y[z_idx] << "\n";  // Отладочное сообщение
            z_idx++;
        }
    }

    mesh.k = 0;
    for (int p = 0; p < mesh.L; ++p) {
        int x_start = mesh.W[0][p], x_end = mesh.W[1][p];
        int y_start = mesh.W[2][p], y_end = mesh.W[3][p];
        mesh.k += (x_end - x_start) * (y_end - y_start);
    }

    mesh.nvkat.resize(mesh.k);
    mesh.nvtr.resize(mesh.k * NUM_NODES_PER_ELEMENT, -1);
    int elem_idx = 0;
    for (int p = 0; p < mesh.L; ++p) {
        for (int i = mesh.W[2][p]; i < mesh.W[3][p]; ++i) {
            for (int j = mesh.W[0][p]; j < mesh.W[1][p]; ++j) {
                std::cout << "elem_idx: " << elem_idx << ", mesh.k: " << mesh.k << "\n";  // Отладочное сообщение
                if (elem_idx >= mesh.k) throw std::runtime_error("Выход за пределы массива.");
                mesh.nvkat[elem_idx] = p;
                for (int m = 0; m < NUM_NODES_PER_ELEMENT; ++m) {
                    int row = i + (m / 3);  // Убираем умножение на 3
                    int col = j + (m % 3);  // Убираем умножение на 3
                    int global_index = row * mesh.xwn + col;  // Используем mesh.xwn вместо (3 * mesh.xwn - 2)
                    std::cout << "global_index: " << global_index << ", mesh.n: " << mesh.n << "\n";  // Отладочное сообщение
                    if (row >= mesh.ywn || col >= mesh.xwn) {
                        mesh.nvtr[elem_idx * NUM_NODES_PER_ELEMENT + m] = -1;
                        mesh.fict[global_index] = false;
                    } else {
                        mesh.nvtr[elem_idx * NUM_NODES_PER_ELEMENT + m] = global_index;
                    }
                }
                elem_idx++;
            }
        }
    }
}

void loadbc(const std::string& filename, std::vector<int>& indices, std::vector<std::pair<int, int>>& ranges, std::vector<double>& values) {
    std::ifstream fp(filename);
    if (!fp) throw std::runtime_error("Ошибка открытия файла " + filename);

    int count;
    fp >> count;
    indices.resize(count);
    ranges.resize(count);
    values.resize(count);
    for (int i = 0; i < count; ++i) {
        int index, start, end;
        double value;
        fp >> index >> start >> end >> value;
        indices[i] = index - 1;
        ranges[i] = {start - 1, end - 1};
        values[i] = value;
    }
}

void allocate_structure(const Mesh& mesh, CSRMatrix& A) {
    std::vector<int> temp(mesh.n, 0);
    for (int ielem = 0; ielem < mesh.k; ++ielem) {
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node_i = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node_i == -1 || !mesh.fict[node_i]) continue;
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
                int node_j = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + j];
                if (node_j == -1 || !mesh.fict[node_j] || node_j > node_i) continue;
                temp[node_i]++;
            }
        }
    }

    A.row_ptr.resize(mesh.n + 1);
    A.row_ptr[0] = 0;
    for (int i = 0; i < mesh.n; ++i) A.row_ptr[i + 1] = A.row_ptr[i] + temp[i];
    A.nnz = A.row_ptr[mesh.n];

    A.col_idx.resize(A.nnz);
    A.values.resize(A.nnz);

    std::vector<int> current(mesh.n, 0);
    for (int ielem = 0; ielem < mesh.k; ++ielem) {
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int row = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (row == -1 || !mesh.fict[row]) continue;
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
                int col = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + j];
                if (col == -1 || !mesh.fict[col] || col > row) continue;
                int pos = A.row_ptr[row] + current[row]++;
                A.col_idx[pos] = col;
            }
        }
    }
}

void add_local_to_global(int ielem, const std::vector<double>& Ak, const std::vector<double>& Bk, const Mesh& mesh, CSRMatrix& A, std::vector<double>& b) {
    for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
        int row = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
        if (row == -1 || !mesh.fict[row]) continue;
        for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
            int col = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + j];
            if (col == -1 || !mesh.fict[col]) continue;
            auto it = std::lower_bound(A.col_idx.begin() + A.row_ptr[row], A.col_idx.begin() + A.row_ptr[row + 1], col);
            if (it != A.col_idx.begin() + A.row_ptr[row + 1] && *it == col) {
                A.values[std::distance(A.col_idx.begin(), it)] += Ak[i * NUM_NODES_PER_ELEMENT + j];
            }
        }
        b[row] += Bk[i];
    }

    // Дополнительные отладочные сообщения
    if (ielem == 0) {
        std::cout << "Первые 5 строк глобальной матрицы A:\n";
        for (int i = 0; i < 5; ++i) {
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                std::cout << A.values[j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "Первые 5 значений вектора правой части b:\n";
        for (int i = 0; i < 5; ++i) {
            std::cout << b[i] << " ";
        }
        std::cout << "\n";
    }
}


void LOS_solve(const CSRMatrix& A, const std::vector<double>& b, std::vector<double>& x, int max_iter, double tol) {
    int n = A.n;
    std::vector<double> r(n, 0.0), z(n, 0.0), p(n, 0.0), Ap(n, 0.0), Az(n, 0.0);
    x.assign(n, 0.0);
    r = b;
    z = r;
    p = r;
    double r_norm = 0.0;
    for (int i = 0; i < n; ++i) r_norm += r[i] * r[i];
    r_norm = std::sqrt(r_norm);
    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            Ap[i] = 0.0;
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                Ap[i] += A.values[j] * p[A.col_idx[j]];
            }
        }
        double alpha = 0.0, pAp = 0.0;
        for (int i = 0; i < n; ++i) {
            alpha += r[i] * z[i];
            pAp += p[i] * Ap[i];
        }
        if (std::abs(pAp) < EPSILON) {
            throw std::runtime_error("pAp близко к нулю.");
        }
        alpha /= pAp;
        for (int i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        double r_new_norm = 0.0;
        for (int i = 0; i < n; ++i) r_new_norm += r[i] * r[i];
        r_new_norm = std::sqrt(r_new_norm);
        std::cout << "Итерация " << iter + 1 << ": r_norm = " << r_new_norm << "\n";
        if (r_new_norm < tol * r_norm) {
            std::cout << "Сходимость за " << iter + 1 << " итераций.\n";
            break;
        }
        for (int i = 0; i < n; ++i) {
            Az[i] = 0.0;
            for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
                Az[i] += A.values[j] * z[A.col_idx[j]];
            }
        }
        double beta = 0.0, AzAz = 0.0;
        for (int i = 0; i < n; ++i) {
            beta += r[i] * Az[i];
            AzAz += Az[i] * Az[i];
        }
        if (std::abs(AzAz) < EPSILON) {
            throw std::runtime_error("AzAz близко к нулю.");
        }
        beta /= AzAz;
        for (int i = 0; i < n; ++i) {
            z[i] = r[i] - beta * Az[i];
            p[i] = z[i] + beta * p[i];
        }
    }
}
void assemble_system(const Mesh& mesh, const BoundaryCondition& bc, CSRMatrix& A, std::vector<double>& b) {
    allocate_structure(mesh, A);
    b.assign(mesh.n, 0.0);

    for (int ielem = 0; ielem < mesh.k; ++ielem) {
        std::vector<double> element_x(NUM_NODES_PER_ELEMENT), element_y(NUM_NODES_PER_ELEMENT);
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node == -1 || !mesh.fict[node]) continue;
            element_x[i] = mesh.xy_x[node];
            element_y[i] = mesh.xy_y[node];
        }

        // Вычисляем шаг сетки для текущего элемента
        double hx = element_x[2] - element_x[0];
        double hy = element_y[6] - element_y[0];
        if (hx <= 0 || hy <= 0) throw std::runtime_error("Нулевой или отрицательный шаг сетки.");

        double lambdak = lambda_func(element_x[4], element_y[4]);
        double gammak = 0.0;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            if (mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i] == -1) continue;
            double xi = (element_x[i] - element_x[0]) / hx, eta = (element_y[i] - element_y[0]) / hy;
            gammak += gamma_func(element_x[i], element_y[i]) * bilinear_basis_function(xi, eta, i % 4);
        }
        gammak /= NUM_NODES_PER_ELEMENT;

        double k1 = lambdak * hy / hx, k2 = lambdak * hx / hy, k3 = gammak * hx * hy;

        std::vector<double> Ak(NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT);
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
                Ak[i * NUM_NODES_PER_ELEMENT + j] = k1 * calculate_G(i, j, element_x[4], element_y[4]) +
                                                     k2 * calculate_C(i, j, element_x[4], element_y[4]) +
                                                     k3 * calculate_Ck(i, j, element_x[4], element_y[4]);
            }
        }

        std::vector<double> Bk(NUM_NODES_PER_ELEMENT);
        double k_area = hx * hy;
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            int node = mesh.nvtr[ielem * NUM_NODES_PER_ELEMENT + i];
            if (node == -1) continue;
            Bk[i] = f_func(element_x[i], element_y[i]) * k_area;
        }

        // Отладочное сообщение: выводим локальную матрицу и вектор правой части
        std::cout << "Локальная матрица Ak для элемента " << ielem << ":\n";
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; ++j) {
                std::cout << Ak[i * NUM_NODES_PER_ELEMENT + j] << " ";
            }
            std::cout << "\n";
        }

        std::cout << "Локальный вектор правой части Bk для элемента " << ielem << ":\n";
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; ++i) {
            std::cout << Bk[i] << " ";
        }
        std::cout << "\n";

        add_local_to_global(ielem, Ak, Bk, mesh, A, b);
    }

    // Применяем граничные условия
    for (size_t i = 0; i < bc.kt1.size(); ++i) {
        int ind = bc.l1[i].first;
        if (ind < 0 || ind >= mesh.n) throw std::runtime_error("Некорректный индекс узла.");
        double u_val = bc.kt1_values[i];

        // Отладочное сообщение: выводим изменения в глобальной матрице и векторе правой части
        std::cout << "Применяем граничное условие Дирихле для узла " << ind << " со значением " << u_val << "\n";

        for (int j = A.row_ptr[ind]; j < A.row_ptr[ind + 1]; ++j) A.values[j] = 0.0;
        for (int j = A.row_ptr[ind]; j < A.row_ptr[ind + 1]; ++j) {
            if (A.col_idx[j] == ind) {
                A.values[j] = 1.0;
                break;
            }
        }
        b[ind] = u_val;
    }
}
int main() {
    try {
        std::setlocale(LC_ALL, "Russian");
        Mesh mesh;
        BoundaryCondition bc;
        loadnet("st.txt", mesh);
        loadbc("ku1.txt", bc.kt1, bc.l1, bc.kt1_values);
        loadbc("ku2.txt", bc.nvk2, bc.nvr2, bc.nvk2_values);
        loadbc("ku3.txt", bc.nvk3, bc.nvr3, bc.nvk3_values);
        CSRMatrix A(mesh.n);
        std::vector<double> b(mesh.n, 0.0);
        assemble_system(mesh, bc, A, b);
        std::vector<double> x(mesh.n, 0.0);
        LOS_solve(A, b, x, MAX_ITER, EPSILON); 
        std::ofstream fp("q.txt");
        if (!fp) throw std::runtime_error("Ошибка открытия файла q.txt.");
        for (const auto& val : x) fp << val << "\t";
        fp.close();
        std::cout << "Решение записано в файл q.txt\n";
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << "\n";
        return EXIT_FAILURE;
    }
    return 0;
}