#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <locale.h>

// Константы
#define NUM_NODES_PER_ELEMENT 16
#define MAX_ITER 100000
#define EPSILON 1e-12

/////////////////////////////
////        Data        /////
/////////////////////////////

// Глобальная матрица и вектор правой части в формате CSR
double* ggl; // Нижняя треугольная часть матрицы (разреженная)
double* ggu; // Верхняя треугольная часть матрицы (разреженная)
double* di;  // Диагональные элементы матрицы
int* ig;     // Столбцовые индексы ненулевых элементов
int* jg;     // Столбцовые индексы ненулевых элементов
int gn;      // Размерность матрицы
int gk;      // Общее количество ненулевых элементов

// Коэффициенты для LU разложения
double* l_coeff; // Нижние коэффициенты L
double* u_coeff; // Верхние коэффициенты U

// Векторы правой части и решение
double* b; // Вектор правой части
double* q; // Вектор решения

// Для сетки
int xwn, ywn; // Количество узлов по X и Y
double* Xw;    // Координаты узлов по X
double* Yw;    // Координаты узлов по Y
int L;         // Число подобластей
int* W[4];     // Координаты подобластей (левый, правый, нижний, верхний)

int** nvtr;    // Глобальные номера узлов для каждого КЭ [k][NUM_NODES_PER_ELEMENT]
int* nvkat;    // Номера подобластей для каждого КЭ [k]
double* xy[2]; // Координаты узлов [2][n]
bool* fict;    // Флаг фиктивных узлов [n]
int n;         // Общее количество узлов
int k;         // Количество конечных элементов

// Краевые условия
int* kt1;      // Номера краевых условий 1-го рода [kt1_n]
int** l1;      // Границы краевых условий 1-го рода [2][kt1_n]
int kt1_n;     // Количество краевых условий 1-го рода

int* nvk2;     // Номера краевых условий 2-го рода [nvk2_n]
int** nvr2;    // Границы краевых условий 2-го рода [2][nvk2_n]
int nvk2_n;    // Количество краевых условий 2-го рода

int* nvk3;     // Номера краевых условий 3-го рода [nvk3_n]
int** nvr3;    // Границы краевых условий 3-го рода [2][nvk3_n]
int nvk3_n;    // Количество краевых условий 3-го рода

// Локальные матрицы и векторы
double G[NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT];
double C[NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT];
double Ck[NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT];
double Ak[NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT];
double Bk[NUM_NODES_PER_ELEMENT];

/////////////////////////////
////     Prototypes     /////
/////////////////////////////

// Параметры уравнения
double lambda_func(double x, double y);
double gamma_func(double x, double y);
double f_func(double x, double y);

// Краевые условия и их учет
double ug(double x, double y, int index);
void bc1cons();
double teta(double x, double y, int index);
void bc2cons();
double beta_func(double x, double y, int index);
double u_beta_func(double x, double y, int index);
void bc3cons();

// Загрузка и обработка сетки
int getnumW(int i, int j);
void loadnet();
void loadbc();

// Формирование локальных матриц и правых частей
void compile_G_C_matrix();
void compile_local_matrix(int ii);
void compile_local_F(int ii);

// Сборка глобальной СЛАУ
void global_matrix_portrait_comp();
void add_local_to_global(int ii);
void compile_global_a_f();

// Решатель СЛАУ
void lu_decomposition();
void solve_forward(double* y, double* b_vec);
void solve_backward(double* x_vec, double* y);
void mmul(double* res, double* v);
double vmul(double* l_vec, double* r_vec);
void los_lu();

// Вывод решения
void unload();

// Главная функция
int main();

/////////////////////////////
////       Implementation    ///
/////////////////////////////

// Параметры уравнения
double lambda_func(double x, double y) // Постоянная на элементе
{
    return 1.0;
}

double gamma_func(double x, double y)
{
    // Разложение по билинейным базисным функциям
    // Для простоты считаем усреднение по узлам
    // Более точное разложение требует интеграции
    return 0.0;
}

double f_func(double x, double y)
{
    // Для теста 1.1: u = x^3 + y^3
    // Δu = 6x + 6y => f = -Δu = -6x - 6y
    // Для теста 1.2: u = sin(pi*x)*sin(pi*y)
    // Δu = -2*pi^2 * sin(pi*x)*sin(pi*y) => f = 2*pi^2 * sin(pi*x)*sin(pi*y)
    // Здесь используем тест 1.1
    return -6.0 * x - 6.0 * y;
}

// Краевые условия 1-го рода (Дирихле)
double ug(double x, double y, int index) // Каталог
{
    // Для теста 1.1: u = x^3 + y^3
    return pow(x, 3) + pow(y, 3);
}

// Краевые условия 2-го рода (Неймана)
double teta(double x, double y, int index) // Каталог
{
    // Пример: Нулевая производная
    return 0.0;
}

// Краевые условия 3-го рода (Робин)
double beta_func(double x, double y, int index) // Каталог
{
    // Пример: Постоянный коэффициент
    return 1.0;
}

double u_beta_func(double x, double y, int index) // Каталог
{
    // Пример: Постоянная функция
    return 0.0;
}

// Загрузка сетки
int getnumW(int i, int j)
{
    // Получаем номер области, в которой находится узел,
    // не учитываются узлы, принадлежащие верхней и правой границам подобластей, а также фиктивные узлы.
    bool on_boundary = false;
    int region = -1;
    for (int p = 0; p < L; p++)
    {
        if (j >= W[0][p] && j < W[1][p] && i >= W[2][p] && i < W[3][p])
        {
            region = p;
            break;
        }
        else if (j == W[1][p] || i == W[3][p])
        {
            on_boundary = true;
        }
    }
    if (region >= 0)
        return region;
    else if (on_boundary)
        return -2;
    else
        return -1;
}

void loadnet() // Загрузка области и формирование сетки
{
    // Открытие файла st.txt
    FILE* fp;
    if (fopen_s(&fp, "st.txt", "r") != 0 || fp == NULL) {
        printf("Ошибка открытия файла st.txt\n");
        exit(1);
    }

    // Чтение xwn и Xw
    if (fscanf_s(fp, "%d", &xwn) != 1) {
        printf("Ошибка чтения xwn из st.txt\n");
        exit(1);
    }
    Xw = (double*)malloc(xwn * sizeof(double));
    if (Xw == NULL) {
        printf("Ошибка выделения памяти для Xw.\n");
        exit(1);
    }
    for (int i = 0; i < xwn; i++) {
        if (fscanf_s(fp, "%lf", &Xw[i]) != 1) {
            printf("Ошибка чтения Xw[%d] из st.txt\n", i);
            exit(1);
        }
    }

    // Чтение ywn и Yw
    if (fscanf_s(fp, "%d", &ywn) != 1) {
        printf("Ошибка чтения ywn из st.txt\n");
        exit(1);
    }
    Yw = (double*)malloc(ywn * sizeof(double));
    if (Yw == NULL) {
        printf("Ошибка выделения памяти для Yw.\n");
        exit(1);
    }
    for (int i = 0; i < ywn; i++) {
        if (fscanf_s(fp, "%lf", &Yw[i]) != 1) {
            printf("Ошибка чтения Yw[%d] из st.txt\n", i);
            exit(1);
        }
    }

    // Чтение L и W
    if (fscanf_s(fp, "%d", &L) != 1) {
        printf("Ошибка чтения L из st.txt\n");
        exit(1);
    }
    for (int p = 0; p < 4; p++) {
        W[p] = (int*)malloc(L * sizeof(int));
        if (W[p] == NULL) {
            printf("Ошибка выделения памяти для W[%d].\n", p);
            exit(1);
        }
    }
    for (int i = 0; i < L; i++) {
        if (fscanf_s(fp, "%d %d %d %d", &W[0][i], &W[1][i], &W[2][i], &W[3][i]) != 4) {
            printf("Ошибка чтения W[%d] из st.txt\n", i);
            exit(1);
        }
        // Корректировка индексации (0-индексация)
        W[0][i]--;
        W[1][i]--;
        W[2][i]--;
        W[3][i]--;
    }
    fclose(fp);

    // Построение сетки
    n = (3 * xwn - 2) * (3 * ywn - 2);
    fict = (bool*)malloc(n * sizeof(bool));
    if (fict == NULL) {
        printf("Ошибка выделения памяти для fict.\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) fict[i] = true;

    // Выделение памяти для координат узлов
    xy[0] = (double*)malloc(n * sizeof(double));
    if (xy[0] == NULL) {
        printf("Ошибка выделения памяти для xy[0].\n");
        exit(1);
    }
    xy[1] = (double*)malloc(n * sizeof(double));
    if (xy[1] == NULL) {
        printf("Ошибка выделения памяти для xy[1].\n");
        exit(1);
    }

    // Формирование координат узлов
    int z_idx = 0;
    for (int i = 0; i < ywn; i++) {
        for (int j = 0; j < xwn; j++) {
            // Внутри каждого прямоугольника 3x3 узла для биквадратичных функций
            double hx = (j < xwn - 1) ? (Xw[j + 1] - Xw[j]) : (Xw[j] - Xw[j - 1]);
            double hy = (i < ywn - 1) ? (Yw[i + 1] - Yw[i]) : (Yw[i] - Yw[i - 1]);

            // Центральные узлы
            for (int mi = 0; mi < 3; mi++) {
                for (int mj = 0; mj < 3; mj++) {
                    if (z_idx >= n) break;
                    xy[0][z_idx] = Xw[j] + (mj / 2.0) * hx;
                    xy[1][z_idx] = Yw[i] + (mi / 2.0) * hy;
                    z_idx++;
                }
            }
        }
    }

    // Вычисление количества конечных элементов
    k = 0;
    for (int i = 0; i < L; i++) {
        int p = W[1][i] - W[0][i];
        int q = W[3][i] - W[2][i];
        k += p * q;
    }

    // Выделение памяти для nvkat и nvtr
    nvkat = (int*)malloc(k * sizeof(int));
    if (nvkat == NULL) {
        printf("Ошибка выделения памяти для nvkat.\n");
        exit(1);
    }

    nvtr = (int**)malloc(k * sizeof(int*));
    if (nvtr == NULL) {
        printf("Ошибка выделения памяти для nvtr.\n");
        exit(1);
    }
    for (int elem = 0; elem < k; elem++) {
        nvtr[elem] = (int*)malloc(NUM_NODES_PER_ELEMENT * sizeof(int));
        if (nvtr[elem] == NULL) {
            printf("Ошибка выделения памяти для nvtr[%d].\n", elem);
            exit(1);
        }
    }

    // Заполнение массивов nvkat и nvtr
    int elem_idx = 0;
    for (int p = 0; p < L; p++) {
        for (int i = W[2][p]; i < W[3][p]; i++) { // По строкам
            for (int j = W[0][p]; j < W[1][p]; j++) { // По столбцам
                if (elem_idx >= k) break;
                nvkat[elem_idx] = p;
                // Предполагаем, что каждый элемент имеет 16 узлов (4x4 сетка)
                // Реализуем простое заполнение, требуется адаптация под конкретную сетку
                for (int m = 0; m < NUM_NODES_PER_ELEMENT; m++) {
                    // Пример заполнения: линейное расположение узлов
                    // Необходимо адаптировать под вашу конкретную локальную нумерацию
                    // Здесь предполагается, что узлы идут слева направо и снизу вверх
                    int row = i * 3 + (m / 4);
                    int col = j * 3 + (m % 4);
                    if (row >= (3 * ywn - 2) || col >= (3 * xwn - 2)) {
                        nvtr[elem_idx][m] = -1; // Фиктивный узел
                        fict[row * (3 * xwn - 2) + col] = false;
                    }
                    else {
                        nvtr[elem_idx][m] = row * (3 * xwn - 2) + col;
                    }
                }
                elem_idx++;
            }
        }
    }

    // Освобождение памяти для W
    for (int p = 0; p < 4; p++) {
        free(W[p]);
    }

    // Инициализация остальных массивов
    // Вектор правой части
    b = (double*)calloc(n, sizeof(double));
    if (b == NULL) {
        printf("Ошибка выделения памяти для вектора b.\n");
        exit(1);
    }

    // Вектор решения
    q = (double*)calloc(n, sizeof(double));
    if (q == NULL) {
        printf("Ошибка выделения памяти для вектора q.\n");
        exit(1);
    }

    // Инициализация коэффициентов L и U
    l_coeff = (double*)calloc(gk, sizeof(double));
    u_coeff = (double*)calloc(gk, sizeof(double));
    if (l_coeff == NULL || u_coeff == NULL) {
        printf("Ошибка выделения памяти для l_coeff или u_coeff.\n");
        exit(1);
    }

    // Отладочный вывод
    printf("Сетка загружена: n=%d, k=%d\n", n, k);
}

// Загрузка краевых условий
void loadbc()
{
    FILE* fp;

    // Загрузка краевых условий 1-го рода
    if (fopen_s(&fp, "ku1.txt", "r") != 0 || fp == NULL) {
        printf("Ошибка открытия файла ku1.txt\n");
        exit(1);
    }
    if (fscanf_s(fp, "%d", &kt1_n) != 1) {
        printf("Ошибка чтения kt1_n из ku1.txt\n");
        exit(1);
    }
    kt1 = (int*)malloc(kt1_n * sizeof(int));
    if (kt1 == NULL) {
        printf("Ошибка выделения памяти для kt1.\n");
        exit(1);
    }
    l1 = (int**)malloc(2 * sizeof(int*));
    if (l1 == NULL) {
        printf("Ошибка выделения памяти для l1.\n");
        exit(1);
    }
    l1[0] = (int*)malloc(kt1_n * sizeof(int));
    l1[1] = (int*)malloc(kt1_n * sizeof(int));
    if (l1[0] == NULL || l1[1] == NULL) {
        printf("Ошибка выделения памяти для l1[0] или l1[1].\n");
        exit(1);
    }
    for (int i = 0; i < kt1_n; i++) {
        if (fscanf_s(fp, "%d %d %d", &kt1[i], &l1[0][i], &l1[1][i]) != 3) {
            printf("Ошибка чтения краевых условий 1-го рода из ku1.txt\n");
            exit(1);
        }
        // Корректировка индексации (0-индексация)
        l1[0][i]--;
        l1[1][i]--;
    }
    fclose(fp);

    // Загрузка краевых условий 2-го рода
    if (fopen_s(&fp, "ku2.txt", "r") != 0 || fp == NULL) {
        printf("Ошибка открытия файла ku2.txt\n");
        exit(1);
    }
    if (fscanf_s(fp, "%d", &nvk2_n) != 1) {
        printf("Ошибка чтения nvk2_n из ku2.txt\n");
        exit(1);
    }
    nvk2 = (int*)malloc(nvk2_n * sizeof(int));
    if (nvk2 == NULL) {
        printf("Ошибка выделения памяти для nvk2.\n");
        exit(1);
    }
    nvr2 = (int**)malloc(2 * sizeof(int*));
    if (nvr2 == NULL) {
        printf("Ошибка выделения памяти для nvr2.\n");
        exit(1);
    }
    nvr2[0] = (int*)malloc(nvk2_n * sizeof(int));
    nvr2[1] = (int*)malloc(nvk2_n * sizeof(int));
    if (nvr2[0] == NULL || nvr2[1] == NULL) {
        printf("Ошибка выделения памяти для nvr2[0] или nvr2[1].\n");
        exit(1);
    }
    for (int i = 0; i < nvk2_n; i++) {
        if (fscanf_s(fp, "%d %d %d", &nvk2[i], &nvr2[0][i], &nvr2[1][i]) != 3) {
            printf("Ошибка чтения краевых условий 2-го рода из ku2.txt\n");
            exit(1);
        }
        // Корректировка индексации (0-индексация)
        nvr2[0][i]--;
        nvr2[1][i]--;
    }
    fclose(fp);

    // Загрузка краевых условий 3-го рода
    if (fopen_s(&fp, "ku3.txt", "r") != 0 || fp == NULL) {
        printf("Ошибка открытия файла ku3.txt\n");
        exit(1);
    }
    if (fscanf_s(fp, "%d", &nvk3_n) != 1) {
        printf("Ошибка чтения nvk3_n из ku3.txt\n");
        exit(1);
    }
    nvk3 = (int*)malloc(nvk3_n * sizeof(int));
    if (nvk3 == NULL) {
        printf("Ошибка выделения памяти для nvk3.\n");
        exit(1);
    }
    nvr3 = (int**)malloc(2 * sizeof(int*));
    if (nvr3 == NULL) {
        printf("Ошибка выделения памяти для nvr3.\n");
        exit(1);
    }
    nvr3[0] = (int*)malloc(nvk3_n * sizeof(int));
    nvr3[1] = (int*)malloc(nvk3_n * sizeof(int));
    if (nvr3[0] == NULL || nvr3[1] == NULL) {
        printf("Ошибка выделения памяти для nvr3[0] или nvr3[1].\n");
        exit(1);
    }
    for (int i = 0; i < nvk3_n; i++) {
        if (fscanf_s(fp, "%d %d %d", &nvk3[i], &nvr3[0][i], &nvr3[1][i]) != 3) {
            printf("Ошибка чтения краевых условий 3-го рода из ku3.txt\n");
            exit(1);
        }
        // Корректировка индексации (0-индексация)
        nvr3[0][i]--;
        nvr3[1][i]--;
    }
    fclose(fp);
}

// Формирование локальных матриц G и C
void compile_G_C_matrix()
{
    // Для простоты заполняем G и C как единичные матрицы
    // В реальной задаче необходимо вычислить интегралы произведений базисных функций и их производных
    for (int i = 0; i < NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT; i++) {
        G[i] = (i % (NUM_NODES_PER_ELEMENT + 1) == 0) ? 1.0 : 0.0; // Диагональные элементы
        C[i] = (i % (NUM_NODES_PER_ELEMENT + 1) == 0) ? 1.0 : 0.0; // Диагональные элементы
        Ck[i] = (i % (NUM_NODES_PER_ELEMENT + 1) == 0) ? 1.0 : 0.0; // Диагональные элементы
    }
}

// Формирование локальной матрицы жесткости для элемента ii
void compile_local_matrix(int ii)
{
    double lambdak = lambda_func(xy[0][nvtr[ii][8]], xy[1][nvtr[ii][8]]);
    double gammak = 0.0;
    for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
        gammak += gamma_func(xy[0][nvtr[ii][i]], xy[1][nvtr[ii][i]]);
    }
    gammak /= NUM_NODES_PER_ELEMENT;

    double hx = xy[0][nvtr[ii][15]] - xy[0][nvtr[ii][0]];
    double hy = xy[1][nvtr[ii][15]] - xy[1][nvtr[ii][0]];
    double k1 = lambdak * hy / hx;
    double k2 = lambdak * hx / hy;
    double k3 = gammak * hx * hy;

    // Формирование локальной матрицы Ak (16x16)
    for (int i = 0; i < NUM_NODES_PER_ELEMENT * NUM_NODES_PER_ELEMENT; i++) {
        Ak[i] = k1 * G[i] + k2 * C[i] + k3 * Ck[i];
    }
}

// Формирование локального вектора правой части для элемента ii
void compile_local_F(int ii)
{
    double fk[NUM_NODES_PER_ELEMENT];
    double hx = xy[0][nvtr[ii][15]] - xy[0][nvtr[ii][0]];
    double hy = xy[1][nvtr[ii][15]] - xy[1][nvtr[ii][0]];
    double k_area = hx * hy;

    for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
        fk[i] = f_func(xy[0][nvtr[ii][i]], xy[1][nvtr[ii][i]]);
    }

    // Bk = Ck * fk * k_area
    for (int m = 0; m < NUM_NODES_PER_ELEMENT; m++) {
        Bk[m] = 0.0;
        for (int n = 0; n < NUM_NODES_PER_ELEMENT; n++) {
            Bk[m] += Ck[m * NUM_NODES_PER_ELEMENT + n] * fk[n];
        }
        Bk[m] *= k_area;
    }
}

// Формирование портрета матрицы в формате CSR
void global_matrix_portrait_comp()
{
    /*
    Формируем массив ig и jg для CSR формата
    ig - указатели на начало каждой строки в массиве jg
    jg - столбцовые индексы ненулевых элементов
    */
    // Инициализация ig
    ig = (int*)malloc((gn + 1) * sizeof(int));
    if (ig == NULL) {
        printf("Ошибка выделения памяти для ig.\n");
        exit(1);
    }

    // Используем временный массив для подсчета количества ненулевых элементов
    int* temp = (int*)calloc(gn, sizeof(int));
    if (temp == NULL) {
        printf("Ошибка выделения памяти для temp.\n");
        exit(1);
    }

    // Подсчет количества ненулевых элементов в каждой строке
    for (int ielem = 0; ielem < k; ielem++) {
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
            int row = nvtr[ielem][i];
            if (!fict[row]) continue;
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; j++) {
                int col = nvtr[ielem][j];
                if (!fict[col]) continue;
                if (col <= row) { // Нижняя треугольная часть
                    temp[row]++;
                }
            }
        }
    }

    // Формирование массива ig
    ig[0] = 0;
    for (int i = 0; i < gn; i++) {
        ig[i + 1] = ig[i] + temp[i];
    }

    // Общее количество ненулевых элементов
    gk = ig[gn];

    // Выделение памяти для jg, ggl и ggu
    jg = (int*)malloc(gk * sizeof(int));
    if (jg == NULL) {
        printf("Ошибка выделения памяти для jg.\n");
        exit(1);
    }
    ggl = (double*)calloc(gk, sizeof(double));
    ggu = (double*)calloc(gk, sizeof(double));
    if (ggl == NULL || ggu == NULL) {
        printf("Ошибка выделения памяти для ggl или ggu.\n");
        exit(1);
    }

    // Используем временный массив для заполнения jg
    int* current = (int*)calloc(gn, sizeof(int));
    if (current == NULL) {
        printf("Ошибка выделения памяти для current.\n");
        exit(1);
    }

    for (int ielem = 0; ielem < k; ielem++) {
        for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
            int row = nvtr[ielem][i];
            if (!fict[row]) continue;
            for (int j = 0; j < NUM_NODES_PER_ELEMENT; j++) {
                int col = nvtr[ielem][j];
                if (!fict[col]) continue;
                if (col <= row) { // Нижняя треугольная часть
                    int pos = ig[row] + current[row];
                    jg[pos] = col;
                    current[row]++;
                }
            }
        }
    }

    // Освобождение временных массивов
    free(temp);
    free(current);
}

// Добавление локальной матрицы и вектора в глобальную матрицу и вектор
void add_local_to_global(int ii)
{
    // Добавление локальной матрицы Ak в глобальную матрицу
    for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
        int row = nvtr[ii][i];
        if (!fict[row]) continue;
        for (int j = 0; j < NUM_NODES_PER_ELEMENT; j++) {
            int col = nvtr[ii][j];
            if (!fict[col]) continue;
            if (col > row) continue; // Нижняя треугольная часть

            // Поиск позиции (row, col) в CSR формате
            int start = ig[row];
            int end = ig[row + 1];
            int pos = -1;
            for (int p = start; p < end; p++) {
                if (jg[p] == col) {
                    pos = p;
                    break;
                }
            }
            if (pos != -1) {
                ggl[pos] += Ak[i * NUM_NODES_PER_ELEMENT + j];
                ggu[pos] += Ak[i * NUM_NODES_PER_ELEMENT + j];
            }
            else {
                printf("Ошибка: Элемент (%d, %d) не найден в матрице.\n", row, col);
            }
        }
    }

    // Добавление локального вектора Bk в глобальный вектор b
    for (int i = 0; i < NUM_NODES_PER_ELEMENT; i++) {
        int row = nvtr[ii][i];
        if (!fict[row]) continue;
        b[row] += Bk[i];
    }
}

// Сборка глобальной матрицы и вектора правой части
void compile_global_a_f()
{
    // Формирование локальных матриц G, C, Ck
    compile_G_C_matrix();

    // Формирование портрета глобальной матрицы
    global_matrix_portrait_comp();

    // Инициализация диагональных элементов и вектора правой части
    di = (double*)calloc(gn, sizeof(double));
    if (di == NULL) {
        printf("Ошибка выделения памяти для di.\n");
        exit(1);
    }

    // Сборка матрицы и вектора
    for (int ielem = 0; ielem < k; ielem++) {
        compile_local_matrix(ielem);
        compile_local_F(ielem);
        add_local_to_global(ielem);
    }

    // Учет краевых условий
    bc3cons(); // Краевые условия 3-го рода
    bc2cons(); // Краевые условия 2-го рода
    bc1cons(); // Краевые условия 1-го рода
}

// Учет краевых условий 1-го рода (Дирихле)
void bc1cons()
{
    for (int i = 0; i < kt1_n; i++) {
        int ind = l1[0][i];
        if (ind < 0 || ind >= gn) {
            printf("Ошибка: Индекс узла %d вне диапазона [0, %d)\n", ind, gn);
            continue;
        }
        di[ind] = 1.0;
        // Обнуление нижней треугольной части строки
        for (int j = ig[ind]; j < ig[ind + 1]; j++) {
            ggl[j] = 0.0;
        }
        // Обнуление верхней треугольной части столбца
        for (int j = 0; j < gk; j++) {
            if (jg[j] == ind) {
                ggu[j] = 0.0;
            }
        }
        // Установка значения вектора правой части
        b[ind] = ug(xy[0][ind], xy[1][ind], kt1[i]);
    }
}

// Учет краевых условий 2-го рода (Неймана)
void bc2cons()
{
    for (int i = 0; i < nvk2_n; i++) {
        int ind_start = nvr2[0][i];
        int ind_end = nvr2[1][i];
        double theta = teta(xy[0][ind_start], xy[1][ind_start], nvk2[i]);

        // Добавление вклада вектора правой части
        for (int ind = ind_start; ind <= ind_end; ind++) {
            if (ind < 0 || ind >= gn) {
                printf("Ошибка: Индекс узла %d вне диапазона [0, %d)\n", ind, gn);
                continue;
            }
            b[ind] += theta;
        }
    }
}

// Учет краевых условий 3-го рода (Робин)
void bc3cons()
{
    for (int i = 0; i < nvk3_n; i++) {
        int ind_start = nvr3[0][i];
        int ind_end = nvr3[1][i];
        double beta_val = beta_func(xy[0][ind_start], xy[1][ind_start], nvk3[i]);
        double u_beta_val = u_beta_func(xy[0][ind_start], xy[1][ind_start], nvk3[i]);

        // Добавление вклада в диагональные элементы и вектор правой части
        for (int ind = ind_start; ind <= ind_end; ind++) {
            if (ind < 0 || ind >= gn) {
                printf("Ошибка: Индекс узла %d вне диапазона [0, %d)\n", ind, gn);
                continue;
            }
            di[ind] += beta_val;
            b[ind] += beta_val * u_beta_val;
        }
    }
}

// LU разложение (разделение матрицы на L и U)
void lu_decomposition()
{
    // Выделение памяти для L и U коэффициентов
    l_coeff = (double*)calloc(gk, sizeof(double));
    u_coeff = (double*)calloc(gk, sizeof(double));
    if (l_coeff == NULL || u_coeff == NULL) {
        printf("Ошибка выделения памяти для l_coeff или u_coeff.\n");
        exit(1);
    }

    // Инициализация диагональных элементов
    for (int i = 0; i < gn; i++) {
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            if (jg[j] < i) {
                l_coeff[j] = ggl[j] / di[jg[j]];
                ggl[j] = l_coeff[j];
            }
            else if (jg[j] == i) {
                di[i] -= l_coeff[j] * ggu[j];
                if (di[i] == 0.0) {
                    printf("Ошибка: Нулевой диагональный элемент при LU разложении.\n");
                    exit(1);
                }
                u_coeff[j] = ggu[j] / di[i];
                ggu[j] = u_coeff[j];
            }
            else {
                u_coeff[j] = ggu[j];
            }
        }
    }
}

// Прямой ход: L * y = b
void solve_forward(double* y, double* b_vec)
{
    for (int i = 0; i < gn; i++) {
        y[i] = b_vec[i];
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            if (jg[j] < i) {
                y[i] -= l_coeff[j] * y[jg[j]];
            }
        }
        y[i] /= di[i];
    }
}

// Обратный ход: U * x = y
void solve_backward(double* x_vec, double* y)
{
    for (int i = gn - 1; i >= 0; i--) {
        x_vec[i] = y[i];
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            if (jg[j] > i) {
                x_vec[i] -= u_coeff[j] * x_vec[jg[j]];
            }
        }
    }
}

// Умножение матрицы на вектор (A * v)
void mmul(double* res, double* v)
{
    // Инициализация результата нулями
    for (int i = 0; i < gn; i++) {
        res[i] = 0.0;
    }

    // Умножение нижней треугольной части
    for (int i = 0; i < gn; i++) {
        for (int j = ig[i]; j < ig[i + 1]; j++) {
            int col = jg[j];
            res[i] += ggl[j] * v[col];
            if (col != i) { // Верхняя треугольная часть
                res[col] += ggu[j] * v[i];
            }
        }
    }
}

// Скалярное произведение
double vmul(double* l_vec, double* r_vec)
{
    double res = 0.0;
    for (int i = 0; i < gn; i++) {
        res += l_vec[i] * r_vec[i];
    }
    return res;
}

// Метод LOS с LU разложением
void los_lu()
{
    // Выделение памяти для необходимых векторов
    double* x_sol = (double*)calloc(gn, sizeof(double));
    double* r = (double*)calloc(gn, sizeof(double));
    double* z = (double*)calloc(gn, sizeof(double));
    double* p = (double*)calloc(gn, sizeof(double));
    double* y = (double*)calloc(gn, sizeof(double));
    double* tmp = (double*)calloc(gn, sizeof(double));

    if (x_sol == NULL || r == NULL || z == NULL || p == NULL || y == NULL || tmp == NULL) {
        printf("Ошибка выделения памяти для векторов решения.\n");
        exit(1);
    }

    // Инициализация начального приближения
    for (int i = 0; i < gn; i++) {
        x_sol[i] = 0.0;
    }

    // Вычисление начальной невязки r = b - A*x
    mmul(tmp, x_sol);
    for (int i = 0; i < gn; i++) {
        r[i] = b[i] - tmp[i];
    }

    // Решение L * y = r
    solve_forward(y, r);

    // Решение U * z = y
    solve_backward(z, y);

    // Инициализация p = z
    for (int i = 0; i < gn; i++) {
        p[i] = z[i];
    }

    double rho = vmul(r, z);
    double rho_old = rho;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        // Вычисление A*p
        double Ap[10];
        mmul(Ap, p);

        // Вычисление альфа = rho / (p, Ap)
        double alpha = rho / vmul(p, Ap);

        // Обновление решения x = x + alpha * p
        for (int i = 0; i < gn; i++) {
            x_sol[i] += alpha * p[i];
        }

        // Обновление невязки r = r - alpha * Ap
        for (int i = 0; i < gn; i++) {
            r[i] -= alpha * Ap[i];
        }

        // Решение L * y = r
        solve_forward(y, r);

        // Решение U * z = y
        solve_backward(z, y);

        // Вычисление rho_new = (r, z)
        double rho_new = vmul(r, z);

        // Проверка сходимости
        if (sqrt(rho_new) < EPSILON) {
            printf("Сходимость достигнута за %d итераций.\n", iter + 1);
            break;
        }

        // Вычисление бета = rho_new / rho
        double beta_val = rho_new / rho;

        // Обновление направления поиска p = z + beta * p
        for (int i = 0; i < gn; i++) {
            p[i] = z[i] + beta_val * p[i];
        }

        // Обновление rho
        rho = rho_new;

        // Если достигнуто максимальное количество итераций
        if (iter == MAX_ITER - 1) {
            printf("Максимальное количество итераций (%d) достигнуто.\n", MAX_ITER);
        }
    }

    // Сохранение решения
    q = x_sol;

    // Освобождение памяти
    free(r);
    free(z);
    free(p);
    free(y);
    free(tmp);
}

// Вывод решения в файл
void unload()
{
    FILE* fp;
    if (fopen_s(&fp, "q.txt", "w") != 0 || fp == NULL) {
        printf("Ошибка открытия файла q.txt для записи\n");
        exit(1);
    }
    for (int i = 0; i < gn; i++) {
        fprintf(fp, "%.7lf\t", q[i]);
    }
    fclose(fp);
}

// Главная функция
int main()
{
    setlocale(LC_ALL, "RU"); // RU локализация консоли
    // Загрузка сетки
    loadnet();

    // Загрузка краевых условий
    loadbc();

    // Сборка глобальной матрицы и вектора правой части
    compile_global_a_f();

    // Решение СЛАУ методом LOS с LU разложением
    lu_decomposition();
    los_lu();

    // Вывод решения
    unload();

    printf("Решение записано в файл q.txt\n");
    system("notepad q.txt"); // Открыть файл q.txt после выполнения

    // Освобождение выделенной памяти
    free(Xw);
    free(Yw);
    for (int elem = 0; elem < k; elem++) {
        free(nvtr[elem]);
    }
    free(nvtr);
    free(nvkat);
    free(fict);
    free(xy[0]);
    free(xy[1]);
    free(b);
    free(q);
    free(ggl);
    free(ggu);
    free(jg);
    free(di);
    free(l_coeff);
    free(u_coeff);
    free(kt1);
    free(l1[0]);
    free(l1[1]);
    free(l1);
    free(nvk2);
    free(nvr2[0]);
    free(nvr2[1]);
    free(nvr2);
    free(nvk3);
    free(nvr3[0]);
    free(nvr3[1]);
    free(nvr3);

    return 0;
}
