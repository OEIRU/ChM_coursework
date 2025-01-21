#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>

constexpr int BIQUAD_NODES = 9;
constexpr int BILIN_NODES = 4;
constexpr int GAUSS_POINTS = 3;
constexpr double EPS = 1e-12;

const double gauss_points[GAUSS_POINTS] = {-0.7745966692, 0.0, 0.7745966692};
const double gauss_weights[GAUSS_POINTS] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

struct Node {
    double x, y;
    int bc_type = 0;
    double bc_value = 0.0;
    bool is_active = true;
};

struct Element {
    int nodes[BIQUAD_NODES];
    double gamma[BILIN_NODES];
};

struct EdgeBC {
    int node1, node2;
    int func_id;
    double beta = 0.0;
};

struct SolverParams {
    int max_iter = 1000;
    double tolerance = 1e-8;
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
    std::vector<EdgeBC> bc_first, bc_second, bc_third;
    SolverParams params;
    SparseMatrix matrix;
    std::vector<double> rhs;
    std::vector<double> solution;

    double bilinear_basis(double xi, double eta, int i) {
        const double n[BILIN_NODES][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};
        return 0.25*(1 + xi*n[i][0])*(1 + eta*n[i][1]);
    }

    double biquadratic_basis(double xi, double eta, int i) {
        const double n[BIQUAD_NODES][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}, 
                                         {0,-1}, {1,0}, {0,1}, {-1,0}, {0,0}};
        if(i < 4) return 0.25*(xi*n[i][0] + 1)*(eta*n[i][1] + 1)*(xi*n[i][0] + eta*n[i][1] - 1);
        if(i < 8) return 0.5*(1 - xi*xi)*(eta*n[i-4][1] + 1) + 0.5*(1 - eta*eta)*(xi*n[i-4][0] + 1);
        if(i == 8) return (1 - xi*xi)*(1 - eta*eta); 
        return (1 - xi*xi)*(1 - eta*eta);
    }

    void biquadratic_grad(double xi, double eta, int i, double& dx, double& dy) {
        const double eps = 1e-6;
        dx = (biquadratic_basis(xi+eps, eta, i) - biquadratic_basis(xi-eps, eta, i))/(2*eps);
        dy = (biquadratic_basis(xi, eta+eps, i) - biquadratic_basis(xi, eta-eps, i))/(2*eps);
    }

    std::vector<int> get_edge_nodes(int n1, int n2) {
        const int edge_indices[4][3] = {{0, 4, 1}, {1, 5, 2}, {2, 6, 3}, {3, 7, 0}};
        
        for(const auto& elem : elements) {
            for(int e = 0; e < 4; ++e) {
                int a = elem.nodes[edge_indices[e][0]];
                int b = elem.nodes[edge_indices[e][2]];
                if((a == n1 && b == n2) || (a == n2 && b == n1)) {
                    return {elem.nodes[edge_indices[e][0]],
                            elem.nodes[edge_indices[e][1]],
                            elem.nodes[edge_indices[e][2]]};
                }
            }
        }
        return {};
    }

    double calculate_edge_length(int n1, int n2) {
        return std::hypot(nodes[n2].x - nodes[n1].x, nodes[n2].y - nodes[n1].y);
    }

    double get_dirichlet_value(int func_id, double x, double y) {
        return x*x + y*y; // Пример аналитического решения
    }

    double get_neumann_value(int func_id) {
        return -4.0; // Для теста с u = x² + y²
    }

public:

    void init_sparsity() {
        std::vector<std::unordered_map<int, double>> temp(nodes.size());
        
        for(const auto& elem : elements) {
            for(int i = 0; i < BIQUAD_NODES; ++i) {
                int row = elem.nodes[i];
                //if(!nodes[row].is_active) continue;
                
                for(int j = 0; j < BIQUAD_NODES; ++j) {
                    int col = elem.nodes[j];
                    if(!nodes[col].is_active || col > row) continue;
                    temp[row][col] += 0.0;
                }
            }
        }
        
        matrix.size = nodes.size();
        matrix.row_ptr.push_back(0);
        for(int i = 0; i < matrix.size; ++i) {
            matrix.diag.push_back(0.0);
            std::vector<int> cols;
            for(const auto& pair : temp[i]) {
                if(pair.first != i) cols.push_back(pair.first);
            }
            std::sort(cols.begin(), cols.end());
            matrix.col_idx.insert(matrix.col_idx.end(), cols.begin(), cols.end());
            matrix.row_ptr.push_back(matrix.col_idx.size());
        }
        matrix.lower.resize(matrix.col_idx.size(), 0.0);
    }

    void load_nodes(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error opening nodes file\n";
            return;
        }

        double dummy;
        int n_nodes;
        file >> dummy >> dummy >> n_nodes;
        
        nodes.resize(n_nodes);
        for(int i = 0; i < n_nodes; ++i) {
            file >> nodes[i].x >> nodes[i].y;
            std::cout << "Загружено узлов: " << nodes.size() << "\n"; // Должно быть 9

        }
    }

    void load_elems(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Error opening elements file\n";
            return;
        }

        int n_elems;
        file >> n_elems;
        elements.resize(n_elems);
        
        for(auto& elem : elements) {
            for(int i = 0; i < BIQUAD_NODES; ++i) file >> elem.nodes[i];
            for(int i = 0; i < BILIN_NODES; ++i) file >> elem.gamma[i];
        }
    }

    void load_bc(const std::string& filename, int bc_type) {
        std::ifstream file(filename);
        if (!file) return;

        int n_edges;
        file >> n_edges;
        
        for(int i = 0; i < n_edges; ++i) {
            EdgeBC bc;
            file >> bc.node1 >> bc.node2 >> bc.func_id;
            if(bc_type == 3) file >> bc.beta;
            
            auto& target = (bc_type == 1) ? bc_first : 
                          (bc_type == 2) ? bc_second : bc_third;
            target.push_back(bc);
        }
    }

    void apply_boundary_conditions() {
        // Инициализация вектора solution перед использованием
        solution.resize(nodes.size(), 0.0);
        rhs.resize(nodes.size(), 0.0); 

        // Обработка условий Дирихле
        for(const auto& bc : bc_first) {
            auto edge_nodes = get_edge_nodes(bc.node1, bc.node2);
            
            // Проверка на пустые ребра
            if(edge_nodes.empty()) {
                std::cerr << "Предупреждение: ребро " << bc.node1 << "-" << bc.node2 
                          << " не найдено, пропуск условий Дирихле\n";
                continue;
            }

            for(int node : edge_nodes) {
                // Проверка валидности индекса
                if(node < 0 || node >= nodes.size()) {
                    std::cerr << "Ошибка: некорректный индекс узла " << node 
                              << " при установке Дирихле\n";
                    continue;
                }

                // Установка значений
                nodes[node].is_active = false;
                solution[node] = get_dirichlet_value(bc.func_id, 
                    nodes[node].x, nodes[node].y);
            }
        }

        // Нейман
        for (const auto& bc : bc_second) {
            auto edge_nodes = get_edge_nodes(bc.node1, bc.node2);
            if (edge_nodes.empty()) continue;

            // Проверка всех узлов ребра
            for (int node : edge_nodes) {
                if (node < 0 || node >= nodes.size()) {
                    std::cerr << "Ошибка: некорректный индекс узла " << node << "\n";
                    return;
                }
            }

        // Вычисление параметров
        double theta = get_neumann_value(bc.func_id);
        double length = calculate_edge_length(edge_nodes[0], edge_nodes[2]);

        // Гауссово интегрирование
        const int num_gauss = 3;
        for(int gi = 0; gi < num_gauss; ++gi) {
            double xi = gauss_points[gi];
            double weight = gauss_weights[gi];
            double N[3] = {0.5*xi*(xi-1), 1.0-xi*xi, 0.5*xi*(xi+1)};

            for(int i = 0; i < 3; ++i) {
                int node = edge_nodes[i];
                
                // Проверка 3: Активен ли узел и корректен ли индекс?
                if(node < 0 || node >= nodes.size()) {
                    std::cerr << "Ошибка: индекс " << node << " вне диапазона!\n";
                    continue;
                }

                if(nodes[node].is_active) {
                    // Проверка 4: Не вышли ли за пределы rhs?
                    if(node >= rhs.size()) {
                        std::cerr << "Ошибка: rhs.size() = " << rhs.size()
                                  << ", node = " << node << "\n";
                        continue;
                    }
                    
                    rhs[node] += theta * N[i] * weight * (length/2.0);
                }
            }
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

            for(int gi = 0; gi < GAUSS_POINTS; ++gi) {
                for(int gj = 0; gj < GAUSS_POINTS; ++gj) {
                    double xi = gauss_points[gi];
                    double eta = gauss_points[gj];
                    double w = gauss_weights[gi] * gauss_weights[gj];

                    // Якобиан
                    double J[2][2] = {0};
                    for(int k = 0; k < BIQUAD_NODES; ++k) {
                        double dx, dy;
                        biquadratic_grad(xi, eta, k, dx, dy);
                        J[0][0] += dx * nodes[elem.nodes[k]].x;
                        J[0][1] += dx * nodes[elem.nodes[k]].y;
                        J[1][0] += dy * nodes[elem.nodes[k]].x;
                        J[1][1] += dy * nodes[elem.nodes[k]].y;
                    }
                    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

                    // Интерполяция γ
                    double gamma = 0.0;
                    for(int k = 0; k < BILIN_NODES; ++k) {
                        gamma += elem.gamma[k] * bilinear_basis(xi, eta, k);
                    }

                    // Градиенты
                    double dN[BIQUAD_NODES][2];
                    for(int k = 0; k < BIQUAD_NODES; ++k) {
                        double dxi, deta;
                        biquadratic_grad(xi, eta, k, dxi, deta);
                        dN[k][0] = (J[1][1]*dxi - J[0][1]*deta)/detJ;
                        dN[k][1] = (-J[1][0]*dxi + J[0][0]*deta)/detJ;
                    }

                    // Вклад в матрицу и вектор
                    for(int i = 0; i < BIQUAD_NODES; ++i) {
                        for(int j = 0; j < BIQUAD_NODES; ++j) {
                            Ke[i][j] += (dN[i][0]*dN[j][0] + dN[i][1]*dN[j][1] +
                                       gamma * biquadratic_basis(xi,eta,i) * 
                                       biquadratic_basis(xi,eta,j)) * w * detJ;
                        }
                        Fe[i] += (-4.0) * biquadratic_basis(xi, eta, i) * w * std::abs(detJ);  // тут важное
                    }
                }
            }

            // Сборка
            for(int i = 0; i < BIQUAD_NODES; ++i) {
                int row = elem.nodes[i];
                if(!nodes[row].is_active) continue; // Пропуск неактивных узлов
                
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
        solution.resize(nodes.size(), 0.0);
        std::vector<double> r = rhs;
        std::vector<double> p = r;
        std::vector<double> Ap(nodes.size(), 0.0);

        for(int iter = 0; iter < params.max_iter; ++iter) {
            // Ap = A*p
            std::fill(Ap.begin(), Ap.end(), 0.0);
            for(int i = 0; i < matrix.size; ++i) {
                Ap[i] += matrix.diag[i] * p[i];
                for(int j = matrix.row_ptr[i]; j < matrix.row_ptr[i+1]; ++j) {
                    Ap[i] += matrix.lower[j] * p[matrix.col_idx[j]];
                    Ap[matrix.col_idx[j]] += matrix.lower[j] * p[i];
                }
            }

            double alpha = 0.0, pAp = 0.0;
            for(size_t i = 0; i < r.size(); ++i) {
                alpha += r[i] * r[i];
                pAp += p[i] * Ap[i];
            }
            if(std::abs(pAp) < EPS) break;
            alpha /= pAp;

            for(size_t i = 0; i < solution.size(); ++i) {
                solution[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
            }

            double rr = 0.0;
            for(double val : r) rr += val*val;
            if(std::sqrt(rr) < params.tolerance) {
                std::cout << "Converged in " << iter << " iterations\n";
                break;
            }

            double beta = rr / (alpha * pAp);
            for(size_t i = 0; i < p.size(); ++i) {
                p[i] = r[i] + beta * p[i];
            }
        }
    }

    void save_results(const std::string& filename) {
        std::ofstream file(filename);
        file << "x,y,u\n";
        for(size_t i = 0; i < nodes.size(); ++i) {
            file << nodes[i].x << "," << nodes[i].y << "," << solution[i] << "\n";
        }
    }
};

int main() {
    FEMSolver solver;
    
    solver.load_nodes("nodes.txt");
    solver.load_elems("elems.txt");
    solver.load_bc("first.txt", 1);
    solver.load_bc("second.txt", 2);
    solver.load_bc("third.txt", 3);
    
    solver.init_sparsity();
    solver.apply_boundary_conditions();
    solver.assemble_system();
    solver.solve();
    solver.save_results("results.csv");

    return 0;
}