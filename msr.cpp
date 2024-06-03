#include <algorithm>
#include <memory>
#include <string.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include "msr.h"

double f_0(double, double) {
    return 1;
}

double f_1(double x, double) {
    return x;
}

double f_2(double, double y) {
    return y;
}

double f_3(double x, double y) {
    return x + y;
}

double f_4(double x, double y) {
    return sqrt(x*x + y*y);
}

double f_5(double x, double y) {
    return x*x + y*y;
}

double f_6(double x, double y) {
    return exp(x*x - y*y);
}

double f_7(double x, double y) {
    return 1. / (25*(x*x + y*y) + 1);
}

#define FUNC(I, J) do { ij2l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++; } \
                  while (0)

#define F(I, J) (f(a + (I)*hx, с + (J)*hy))

void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k) {
    int i, i1, i2, l, J; double s;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        s = A[i] * x[i];
        l = I[i+1] - I[i];
        J = I[i];
        for (int j = 0; j < l; ++j) {
            s += A[J + j] * x[I[J + j]];
        }

        y[i] = s;
    }
}

void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v1, double* v2, int flag, int p, int k) {
    double w = 1; // гиперпараметр (метод Зейделя при w = 1)
    if (flag) {
        solve_lsystem(n, I, A, v2, v1, w, p, k);
    } else {
        solve_rsystem(n, I, A, v2, v1, w, p, k);
    }

    reduce_sum<int>(p);
}

bool step(int n, double* A, int* I, double* x, double* r, double* u, double* v, double prec, int p, int k) {
    matrix_mult_vector_msr(n, A, I, v, u, p, k); // u = Av
    double c1 = scalar_product(n, u, r, p, k); // C1 = (u, r)
    double c2 = scalar_product(n, u, u, p, k); // C2 = (u, u)
    if (c1 < prec || c2 < prec) {
        return true;
    }

    double tau = c1 / c2;
    mult_sub_vector(n, x, v, tau, p, k); // x -= tau * v
    mult_sub_vector(n, r, u, tau, p, k); // r -= tau * u
    return false;
}

int minimal_residual_msr_matrix(int n, double* A, int* I, double* b, double* x,
    double* r, double* u, double* v, double eps, int maxit, int p, int k) {

    double prec, b_norm2;
    int it;
    b_norm2 = scalar_product(n, b, b, p, k); // (b, b)
    prec = b_norm2 * eps * eps;
    matrix_mult_vector_msr(n, A, I, x, r, p, k); // r = Ax
    mult_sub_vector(n, r, b, 1., p, k); // r -= 1*b
    
    for (it = 0; it < maxit; ++it) {
        apply_preconditioner_msr_matrix(n, A, I, v, r, 0, p, k); // v = np.linalg.solve(L, r)
        if (step(n, A, I, x, r, u, v, prec, p, k)) {
            break;
        }

        matrix_mult_vector_msr(n, A, I, x, u, p, k); // u = Ax
        mult_sub_vector(n, u, b, 1., p, k); // u -= b 
        apply_preconditioner_msr_matrix(n, A, I, v, u, 1, p, k); // v = np.linalg.solve(R, a@x-b)
        if (step(n, A, I, x, r, u, v, prec, p, k)) {
            break;
        }
    }

    if (it >= maxit) {
        return -1;
    }

    return it;
}

int minimal_residual_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, 
    double eps, int maxit, int maxsteps, int p, int k) {

    int step, ret, its = 0;
    for (step = 0; step < maxsteps; ++step) {
        ret = minimal_residual_msr_matrix(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        if (ret >= 0) { its += ret; break; }
        its += maxit;
    }

    if (step >= maxsteps) {
        return -1;
    } 

    return its;
}

void thread_rows(int n, int p, int k, int& i1, int& i2) {
    i1 = n*k; i1 /= p; i2 = n*(k + 1); i2 /= p; 
}

double scalar_product(int n, double* x, double* y, int p, int k) {
    int i1, i2, i; double s = 0;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        s += x[i] * y[i];
    }

    s = reduce_sum_det(p, k, s);
    return s;
}

void mult_sub_vector(int n, double* x, double* y, double tau, int p, int k) {
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        x[i] -= tau * y[i];
    }

    reduce_sum<int>(p);
}

void ij2l(int nx, int, int i, int j, int& l) {
    l = i + j * (nx + 1);
}

void l2ij(int nx, int, int& i, int& j, int l) {
    j = l / (nx + 1);
    i = l - j * (nx + 1);
}

int get_len_msr(int nx, int ny) {
    return (nx + 1)*(ny + 1)
            + 6*(nx - 1)*(ny - 1)
            + 4*(2*(nx-1) + 2*(ny-1))
            + 2*3 + 2*2; 
}

#define FUNC(I, J) do { ij2l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++; } \
                  while (0)

// странно считает 0 1 0 1 2 2 1 1e-10 1000 1
int get_off_diag(int nx, int ny, int i, int j, int* I_ij) {
    int m = 0; int k;
    if (i < nx) { FUNC(i+1, j); }
    if (j > 0) { FUNC(i, j-1); }
    if (i > 0 && j > 0) { FUNC(i-1, j-1); }
    if (i > 0) { FUNC(i - 1, j); }
    if (j < ny) { FUNC(i, j+1); }
    if (i < nx && j < ny) { FUNC(i+1, j+1); }

    return m;
}

int get_len_msr_off_diag(int nx, int ny) {
    int m = 0; int i, j;
    for (i = 0; i <= nx; ++i) {
        for (j = 0; j <= ny; ++j) {
            m += get_off_diag(nx, ny, i, j);
        }
    }    

    return m;
}

int allocate_msr_matrix(int nx, int ny, double** p_A, int** p_I) {
    int diag_len = (nx+1)*(ny+1);
    int off_diag = get_len_msr_off_diag(nx, ny);
    int len = diag_len + off_diag + 1;

    double* A = nullptr; int* I = nullptr;
    A = new double[len]; if (A == nullptr) { return 1; }
    I = new int[len]; if (I == nullptr) { return 2; }
    *p_A = A; *p_I = I;
    return 0;  
}

void fill_I(int nx, int ny, int* I) { 
    int N = (nx+1) * (ny+1);
    int l; int i, j;
    int m; int r = N + 1;
    for (l = 0; l < N; ++l) {
        l2ij(nx, ny, i, j, l);
        I[l] = r;
        m = get_off_diag(nx, ny, i, j, I + r);
        r += m;
    }

    I[l] = r;
}

void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_off_diag) {

    double s = hx*hy;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        *A_diag = s / 2;
        for (int l = 0; l < 6; ++l) {
            A_off_diag[l] = s / 12; 
        }
    } else if (j == 0 && i > 0 && i < nx) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 2*s/24;
    } else if (j == ny && i > 0 && i < nx) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 1*s/24;
    } else if (i == 0 && j > 0 && j < ny) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 2*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 1*s/24;
        A_off_diag[3] = 2*s/24;
    } else if (i == nx && j > 0 && j < ny) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 1*s/24;
    } else if (i == 0 && j == 0) {
        *A_diag = 2*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 2*s/24;
    } else if (i == nx && j == ny) {
        *A_diag = 2*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 1*s/24;
    } else if (i == 0 && j == ny) {
        *A_diag = 1*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
    } else if (i == nx && j == 0) {
        *A_diag = 1*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
    }

    // abort() если сюда дошли
}

void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k) {
    int l1, l2;
    int i, j;
    int N = (nx + 1) * (ny + 1);    
    l1 = N * k; l1 /= p;
    l2 = N * (k + 1); l2 /= p;

    for (int l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        double* A_diag = A + l;
        double* A_off_diag = A + I[l];

        fill_A_ij(nx, ny, hx, hy, i, j, A_diag, A_off_diag);
    }

    reduce_sum<int>(p);
}

// обязательно проверить матрицу построенную на симметричность иначе потом долго ошибку искать
int check_symm(int nx, int ny, int* I, double* A, double eps, int p, int k) {
    int l1, l2, l;
    int N = (nx+1)*(ny+1);
    l1 = k * N; l1 /= p;
    l2 = (k + 1) * N; l2 /= p;

    int err = 0;

    for (l = l1; l < l2; ++l) {
        int m = I[l + 1] - I[l];
        double* A_off_diag = A + I[l];
        for (int q = 0; q < m; ++q) {
            double a_ij = A_off_diag[q];
            int j = I[I[l] + q];
            // найти в j-ой строке внедиагональный элемент с номером l
            int m2 = I[j + 1] - I[j];
            int q2;
            for (q2 = 0; q2 < m2; ++q2) {
                if (I[I[j] + q2] == l) {
                    break;
                }
            }

            if (q2 >= m2) {
                err++;
            } else if (fabs(A[I[j] + q2] - a_ij) > eps) {
                err++;
            } 
        }
    }

    reduce_sum<int>(p, &err, 1);
    return err;
}

double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double)) {
    double w = hx*hy/192;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        return w * (36*F(i, j) + 20*(F(i+0.5, j) + F(i, j-0.5) + F(i-0.5,j-0.5) + F(i-0.5, j) + F(i, j+0.5) + F(i+0.5,j+0.5))
            + 4*(F(i+0.5, j-0.5) + F(i-0.5, j-1) + F(i-1,j-0.5) + F(i-0.5,j+0.5) + F(i+0.5,j+1) + F(i+1,j+0.5))
            + 2*(F(i+1,j) + F(i, j-1) + F(i-1,j-1) + F(i-1,j) + F(i,j+1) + F(i+1,j+1))
        );
    }
    
    if (i > 0 && i < nx && j == 0) {
        return w * ( 
            18*F(i,j) + 10*(F(i+0.5,j) + F(i-0.5,j)) + 20*(F(i,j+0.5) + F(i+0.5,j+0.5))
            + 4*(F(i-0.5,j+0.5) + F(i+0.5,j+1) + F(i+1,j+0.5)) + 1*(F(i-1,j) + F(i+1,j)) + 2*(F(i,j+1) + F(i+1,j+1))
        );
    }

    if (i > 0 && i < nx && j == ny) {
        return w * (
            18*F(i, j) + 10*(F(i+0.5,j) + F(i-0.5,j)) + 20*(F(i,j-0.5) + F(i-0.5,j-0.5))
            + 4*(F(i+0.5,j-0.5) + F(i-0.5,j-1) + F(i-1,j-0.5)) + 1*(F(i-1,j) + F(i+1,j)) + 2*(F(i,j-1) + F(i-1,j-1))
        );
    }

    if (i == 0 && j > 0 && j < ny) {
        return w*(
            18*F(i,j) + 10*(F(i,j-0.5) + F(i,j+0.5)) + 20*(F(i+0.5,j) + F(i+0.5,j+0.5))
            + 4*(F(i+0.5,j-0.5) + F(i+0.5,j+1) + F(i+1,j+0.5)) + 1*(F(i,j-1) + F(i,j+1)) + 2*(F(i+1,j) + F(i+1,j+1))
        );     
    }

    if (i == nx && j > 0 && j < ny) {
        return w * (
            18*F(i, j) + 10*(F(i,j-0.5) + F(i,j+0.5)) + 20*(F(i-0.5,j) + F(i-0.5,j-0.5))
            + 4*(F(i-0.5,j-1) + F(i-1,j-0.5) + F(i-0.5,j+0.5)) + 1*(F(i,j-1) + F(i,j+1)) + 2*(F(i-1,j) + F(i-1,j-1))
        );
    }

    if (i == 0 && j == 0) {
        return w*(
            12*F(i,j) + 10*(F(i+0.5,j) + F(i,j+0.5)) + 20*F(i+0.5,j+0.5)
            + 4*(F(i+1,j+0.5) + F(i+0.5,j+1)) + 1*(F(i+1,j) + F(i,j+1)) + 2*F(i+1,j+1)
        );     
    }    

    if (i == nx && j == ny) {
        return w*(
            12*F(i,j) + 10*(F(i-0.5,j) + F(i,j-0.5)) + 20*F(i-0.5,j-0.5)
            + 4*(F(i-0.5,j-1) + F(i-1,j-0.5)) + 1*(F(i,j-1) + F(i-1,j)) + 2*F(i-1,j-1)
        );   
    }

    if (i == 0 && j == ny) {
        return w*(
            6*F(i, j) + 10*(F(i+0.5,j) + F(i,j-0.5)) + 4*F(i+0.5,j-0.5) + F(i+1,j) + F(i,j-1) 
        );
    }

    if (i == nx && j == 0) {
        return w*(
            6*F(i, j) + 10*(F(i-0.5,j) + F(i,j+0.5)) + 4*F(i-0.5,j+0.5) + F(i-1,j) + F(i,j+1) 
        );      
    }

    return 1e308;
}

void fill_B(int nx, int ny, double hx, double hy, double a, double c, double* B, double (*f)(double, double), int p, int k) {
    int l1, l2;
    int i, j;
    int N = (nx + 1) * (ny + 1);    
    l1 = N * k; l1 /= p;
    l2 = N * (k + 1); l2 /= p;

    for (int l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        B[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f);
    }

    reduce_sum<int>(p);
}

void solve_rsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k) {

    int i1, i2;
    thread_rows(n, p, k, i1, i2);

    for (int i = i2 - 1; i >= i1; --i) {
        int m = I[i + 1] - I[i];
        double s = 0;
        for (int q = 0; q < m; ++q) {
            int j = I[I[i] + q];
            if (j >= i1 && j < i2 && j > i) {
                s += w * x[j] * U[I[i] + q];
            }
        }

        x[i] = (b[i] - s) / U[i];
    }    
}

void solve_lsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k) {

    int i1, i2;
    thread_rows(n, p, k, i1, i2);

    for (int i = i2 - 1; i >= i1; --i) {
        int m = I[i + 1] - I[i];
        double s = 0;
        for (int q = 0; q < m; ++q) {
            int j = I[I[i] + q];
            if (j >= i1 && j < i2 && j < i) {
                s += w * x[j] * U[I[i] + q];
            }
        }

        x[i] = (b[i] - s) / U[i];
    }    
}

double get_cpu_time() {
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1e-6;
}

static std::unique_ptr<double[]> results = nullptr;

int init_reduce_sum(int p) {
    results.reset(new double[p]{});
    memset(results.get(), 0, p * sizeof (double));
    if (results == nullptr) {
        return -1;
    }
    return 0;
}

double reduce_sum_det(int p, int k, double s) {
    double sum = 0; int l;
    results[k] = s;
    reduce_sum<int>(p);
    for (l = 0; l < p; ++l) {
        sum += results[l];
    }

    reduce_sum<int>(p);
    return sum;
}

double r1(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    int N = (nx+1) * (ny+1);
    int i1, i2;
    int i, j;
    thread_rows(N, p, k, i1, i2);
    double res = -1;
    for (int l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        if (i == nx || j == ny) {
            continue;
        }

        double res1 = fabs(f(a + (i + 2./3)*hx, c + (j + 1./3)*hy) - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3);
        double res2 = fabs(f(a + (i + 1./3)*hx, c + (j + 2./3)*hy) - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3);
        res = std::max(res, std::max(res1, res2));
    }

    reduce_sum(p, &res, 1, &max);    
    return res;
} 

double r2(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    int N = (nx+1) * (ny+1);
    int i1, i2;
    int i, j;
    thread_rows(N, p, k, i1, i2);
    double sum = 0;
    for (int l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        if (i == nx || j == ny) {
            continue;
        }

        double res1 = fabs(f(a + (i + 2./3)*hx, c + (j + 1./3)*hy) - (x[l] + x[l + 1] + x[l + 1 + nx + 1]) / 3);
        double res2 = fabs(f(a + (i + 1./3)*hx, c + (j + 2./3)*hy) - (x[l] + x[l + nx + 1] + x[l + 1 + nx + 1]) / 3);
        sum += res1 + res2;
    }

    sum = reduce_sum_det(p, k, sum);    
    return hx*hy*sum/2;
} 

double r3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    int N = (nx+1) * (ny+1);
    int i1, i2;
    int i, j;
    thread_rows(N, p, k, i1, i2);
    double res = -1;
    for (int l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        res = std::max(res, fabs(f(a + i*hx, c + j*hy) - x[l]));
    }

    reduce_sum(p, &res, 1, &max);    
    return res;
}

double r4(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    int N = (nx+1) * (ny+1);
    int i1, i2;
    int i, j;
    thread_rows(N, p, k, i1, i2);
    double sum = 0;
    for (int l = i1; l < i2; ++l) {
        l2ij(nx, ny, i, j, l);
        sum += fabs(f(a + i*hx, c + j*hy) - x[l]);
    }

    sum = reduce_sum_det(p, k, sum);    
    return hx*hy*sum;
}

void* solution(void *ptr) {
    Args *args = (Args*) ptr;
    pthread_cond_t *cond = args->cond;
    pthread_mutex_t *mutex = args->mutex;
    int kk = args->k;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (kk % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    while (!(args->close_window)) {
        if (!args->ready) {
            double *A = args->A;
            int *I = args->I;
            double *B = args->B;
            double *x = args->x;
            double *u = args->u;
            double *v = args->v;
            double *r = args->r;
            double a = args->a, b = args->b;
            double c = args->c, d = args->d;
            double eps = args->eps;
            int nx = args->nx, ny = args->ny;
            int k = args->k;
            int p = args->p;
            int maxit = args->maxit;
            int maxsteps = 300; // гиперпараметр
            double (*f) (double, double) = args->f;

            if (k == 0) init_reduce_sum(p);
            reduce_sum<int>(p);

            double hx = (b - a) / nx;
            double hy = (d - c) / ny;
            int N = (nx + 1) * (ny + 1);

            fill_A(nx, ny, hx, hy, I, A, p, k);
            fill_B(nx, ny, hx, hy, a, c, B, f, p, k);

            args->t1 = get_cpu_time();
            int its = minimal_residual_msr_matrix_full(N, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, k);
            args->t1 = get_cpu_time() - args->t1;
            args->its = its;

            args->t2 = get_cpu_time();
            double res_1 = r1(nx, ny, a, c, hx, hy, x, f, p, k);
            double res_2 = r2(nx, ny, a, c, hx, hy, x, f, p, k);
            double res_3 = r3(nx, ny, a, c, hx, hy, x, f, p, k);
            double res_4 = r4(nx, ny, a, c, hx, hy, x, f, p, k);
            args->t2 = get_cpu_time() - args->t2;

            args->res_1 = res_1;
            args->res_2 = res_2;
            args->res_3 = res_3;
            args->res_4 = res_4;

            args->ready = true;
            reduce_sum<int>(p);

            pthread_mutex_lock(mutex);
            while (!(args->ready)) {
                pthread_cond_wait(cond, mutex);
            }
            pthread_mutex_unlock(mutex);
        }
    }

    args->ready = true;
    return nullptr;
}
