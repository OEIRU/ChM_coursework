#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <array> // Добавить в заголовки


constexpr int BIQUAD_NODES = 9;
constexpr int BILIN_NODES = 4;
constexpr int GAUSS_POINTS = 3;
constexpr double EPS = 1e-12;

const double gauss_points[GAUSS_POINTS] = {-0.7745966692, 0.0, 0.7745966692};
const double gauss_weights[GAUSS_POINTS] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

struct Node {
    double x, y;
    int bc_type;        // Тип краевого условия (целое)
    double bc_value;    // Значение краевого условия
    bool is_active;
    int id;
    
    // Конструктор для правильной инициализации
    Node(double x_, double y_, int bt, double bv, bool ia, int id_) 
        : x(x_), y(y_), bc_type(bt), bc_value(bv), is_active(ia), id(id_) {}
};

struct Element {
    std::array<int, BIQUAD_NODES> nodes;
    double gamma[BILIN_NODES]; // Коэффициенты γ в билинейных узлах
};

struct SparseMatrix {
    std::vector<double> diag;
    std::vector<double> lower;
    std::vector<int> row_ptr;
    std::vector<int> col_idx;
    int size;
};

class FEMSolver {
private:
    std::vector<Node> nodes;
    std::vector<Element> elements;
    SparseMatrix matrix;
    std::vector<double> rhs;
    std::vector<double> solution;

    // Билинейные базисные функции
    double bilinear_basis(double xi, double eta, int i) const {
        const double n[BILIN_NODES][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};
        return 0.25*(1 + xi*n[i][0])*(1 + eta*n[i][1]);
    }

    // Биквадратичные базисные функции
    double biquadratic_basis(double xi, double eta, int i) const {
        const double n[BIQUAD_NODES][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}, 
                                         {0,-1}, {1,0}, {0,1}, {-1,0}, {0,0}};
        double x = n[i][0], y = n[i][1];
        
        if(i < 4) { // Угловые узлы
            return 0.25*(1 + xi*x)*(1 + eta*y)*(xi*x + eta*y - 1);
        } 
        else if(i < 8) { // Реберные узлы
            if(i % 2 == 0) // Вертикальные
                return 0.5*(1 - xi*xi)*(1 + eta*y);
            else // Горизонтальные
                return 0.5*(1 - eta*eta)*(1 + xi*x);
        } 
        else { // Центральный узел
            return (1 - xi*xi)*(1 - eta*eta);
        }
    }

    // Градиенты биквадратичных функций
    void biquadratic_grad(double xi, double eta, int i, double& dx, double& dy) const {
        const double n[BIQUAD_NODES][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}, 
                                         {0,-1}, {1,0}, {0,1}, {-1,0}, {0,0}};
        double x = n[i][0], y = n[i][1];
        
        if(i < 4) {
            dx = 0.25*x*(1 + eta*y)*(2*xi*x + eta*y);
            dy = 0.25*y*(1 + xi*x)*(2*eta*y + xi*x);
        }
        else if(i < 8) {
            if(i % 2 == 0) { // Вертикальные
                dx = -xi*(1 + eta*y);
                dy = 0.5*(1 - xi*xi)*y;
            } else { // Горизонтальные
                dx = 0.5*(1 - eta*eta)*x;
                dy = -eta*(1 + xi*x);
            }
        }
        else { // Центральный
            dx = -2*xi*(1 - eta*eta);
            dy = -2*eta*(1 - xi*xi);
        }
    }

    // Вычисление якобиана
    void jacobian(double xi, double eta, const Element& elem, 
                 double J[2][2], double& detJ) const {
        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
        
        for(int k = 0; k < BIQUAD_NODES; ++k) {
            double dxi, deta;
            biquadratic_grad(xi, eta, k, dxi, deta);
            J[0][0] += dxi * nodes[elem.nodes[k]].x;
            J[0][1] += dxi * nodes[elem.nodes[k]].y;
            J[1][0] += deta * nodes[elem.nodes[k]].x;
            J[1][1] += deta * nodes[elem.nodes[k]].y;
        }
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    }

public:
    void init_sparsity() {
        std::vector<std::unordered_map<int, bool>> conn(nodes.size());
        
        // Построение связей
        for(const auto& elem : elements) {
            for(int i = 0; i < BIQUAD_NODES; ++i) {
                int row = elem.nodes[i];
                for(int j = 0; j < BIQUAD_NODES; ++j) {
                    int col = elem.nodes[j];
                    if(col <= row) conn[row].insert({col, true});
                }
            }
        }
        
        // Формирование CSR формата
        matrix.size = nodes.size();
        matrix.row_ptr.push_back(0);
        for(int i = 0; i < matrix.size; ++i) {
            matrix.diag.push_back(0.0);
            std::vector<int> cols;
            for(auto& p : conn[i]) cols.push_back(p.first);
            std::sort(cols.begin(), cols.end());
            matrix.col_idx.insert(matrix.col_idx.end(), cols.begin(), cols.end());
            matrix.row_ptr.push_back(matrix.col_idx.size());
        }
        matrix.lower.resize(matrix.col_idx.size(), 0.0);
    }

    void apply_boundary_conditions() {
        // Проставляем флаги активных узлов
        for(auto& node : nodes) {
            node.is_active = (node.bc_type != 1);
        }
    }

    void enforce_dirichlet() {
        for(size_t i = 0; i < nodes.size(); ++i) {
            if(nodes[i].bc_type == 1) {
                // Обнуляем строку
                matrix.diag[i] = 1.0;
                for(int j = matrix.row_ptr[i]; j < matrix.row_ptr[i+1]; ++j)
                    matrix.lower[j] = 0.0;
                
                // Правая часть
                rhs[i] = nodes[i].bc_value;
            }
        }
    }

    void assemble_system() {
        rhs.assign(nodes.size(), 0.0);
        matrix.diag.assign(nodes.size(), 0.0);
        matrix.lower.assign(matrix.lower.size(), 0.0);

        for(const auto& elem : elements) {
            double Ke[BIQUAD_NODES][BIQUAD_NODES] = {0};
            double Fe[BIQUAD_NODES] = {0};

            // Интегрирование по Гауссу
            for(int gi = 0; gi < GAUSS_POINTS; ++gi) {
                for(int gj = 0; gj < GAUSS_POINTS; ++gj) {
                    double xi = gauss_points[gi];
                    double eta = gauss_points[gj];
                    double w = gauss_weights[gi] * gauss_weights[gj];

                    // Якобиан
                    double J[2][2], detJ;
                    jacobian(xi, eta, elem, J, detJ);
                    detJ = fabs(detJ);

                    // Градиенты базисных функций
                    double dN[BIQUAD_NODES][2];
                    for(int k = 0; k < BIQUAD_NODES; ++k) {
                        double dxi, deta;
                        biquadratic_grad(xi, eta, k, dxi, deta);
                        dN[k][0] = (J[1][1]*dxi - J[0][1]*deta)/detJ;
                        dN[k][1] = (-J[1][0]*dxi + J[0][0]*deta)/detJ;
                    }

                    // Коэффициент γ (билинейная интерполяция)
                    double gamma = 0.0;
                    for(int k = 0; k < BILIN_NODES; ++k)
                        gamma += elem.gamma[k] * bilinear_basis(xi, eta, k);

                    // Вклад в матрицу и вектор
                    for(int i = 0; i < BIQUAD_NODES; ++i) {
                        for(int j = 0; j < BIQUAD_NODES; ++j) {
                            Ke[i][j] += (dN[i][0]*dN[j][0] + dN[i][1]*dN[j][1]
                                      + gamma * biquadratic_basis(xi,eta,i) 
                                              * biquadratic_basis(xi,eta,j))
                                      * w * detJ;
                        }
                        Fe[i] += (-4.0) * biquadratic_basis(xi,eta,i) * w * detJ;
                    }
                }
            }

            // Сборка в глобальную систему
            for(int i = 0; i < BIQUAD_NODES; ++i) {
                int row = elem.nodes[i];
                if(!nodes[row].is_active) continue;

                rhs[row] += Fe[i];
                
                for(int j = 0; j < BIQUAD_NODES; ++j) {
                    int col = elem.nodes[j];
                    if(col > row || !nodes[col].is_active) continue;

                    if(col == row) {
                        matrix.diag[row] += Ke[i][j];
                    } else {
                        auto it = std::lower_bound(
                            matrix.col_idx.begin() + matrix.row_ptr[row],
                            matrix.col_idx.begin() + matrix.row_ptr[row+1],
                            col
                        );
                        int idx = it - matrix.col_idx.begin();
                        matrix.lower[idx] += Ke[i][j];
                    }
                }
            }
        }
    }

    void solve() {
        int n = nodes.size();
        solution.assign(n, 0.0);
        std::vector<double> r = rhs;
        std::vector<double> p = r;
        std::vector<double> Ap(n, 0.0);
        double rsold = 0.0;

        for(int i = 0; i < n; ++i) rsold += r[i]*r[i];

        for(int iter = 0; iter < 1000; ++iter) {
            // Ap = A*p
            std::fill(Ap.begin(), Ap.end(), 0.0);
            for(int i = 0; i < n; ++i) {
                Ap[i] += matrix.diag[i] * p[i];
                for(int j = matrix.row_ptr[i]; j < matrix.row_ptr[i+1]; ++j) {
                    Ap[i] += matrix.lower[j] * p[matrix.col_idx[j]];
                    Ap[matrix.col_idx[j]] += matrix.lower[j] * p[i];
                }
            }

            // alpha = rsold / (p'*Ap)
            double alpha = 0.0, pAp = 0.0;
            for(int i = 0; i < n; ++i) {
                alpha += r[i] * r[i];
                pAp += p[i] * Ap[i];
            }
            if(fabs(pAp) < 1e-12) break;
            alpha /= pAp;

            // Обновление решения и невязки
            for(int i = 0; i < n; ++i) {
                solution[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }

            // Проверка сходимости
            double rsnew = 0.0;
            for(int i = 0; i < n; ++i) rsnew += r[i] * r[i];
            if(sqrt(rsnew) < 1e-8) break;

            // Обновление направления
            double beta = rsnew / rsold;
            rsold = rsnew;
            for(int i = 0; i < n; ++i) 
                p[i] = r[i] + beta * p[i];
        }
    }

    void save_results(const std::string& filename) {
        std::ofstream file(filename);
        file << "x,y,u\n";
        for(const auto& node : nodes) {
            file << node.x << "," << node.y << "," << solution[node.id] << "\n";
        }
    }

    // Загрузка тестовых данных
    void load_test_data() {
    nodes = {
        Node(0.0, 0.0, 1, 0.0, false, 0),
        Node(2.0, 0.0, 1, 8.0, false, 1),
        Node(2.0, 2.0, 1, 16.0, false, 2),
        Node(0.0, 2.0, 1, 8.0, false, 3),
        Node(1.0, 0.0, 0, 2.0, true, 4),
        Node(2.0, 1.0, 0, 10.0, true, 5),
        Node(1.0, 2.0, 0, 10.0, true, 6),
        Node(0.0, 1.0, 3, -7.21762, true, 7), 
        Node(1.0, 1.0, 0, 11.0671, true, 8)
    };

    Element elem;
    elem.nodes = {0,1,2,3,4,5,6,7,8}; // Теперь работает с std::array
    for(int i = 0; i < BILIN_NODES; ++i) elem.gamma[i] = 1.0;
    elements.push_back(elem);
    }
};

int main() {
    FEMSolver solver;
    solver.load_test_data();
    
    solver.init_sparsity();
    solver.apply_boundary_conditions();
    solver.assemble_system();
    solver.enforce_dirichlet();
    solver.solve();
    solver.save_results("results.csv");

    std::cout << "Решение сохранено в results.csv" << std::endl;
    return 0;
}