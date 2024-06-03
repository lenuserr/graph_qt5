#include <string.h>
#include "gui_data.h"
#include "msr.h"

void gui_data::read_data(char* argv[]) {
    sscanf (argv[1], "%lf", &a);
    sscanf (argv[2], "%lf", &b);
    sscanf (argv[3], "%lf", &c);
    sscanf (argv[4], "%lf", &d);
    sscanf (argv[5], "%d", &nx);
    sscanf (argv[6], "%d", &ny);
    sscanf (argv[7], "%d", &mx);
    sscanf (argv[8], "%d", &my);
    sscanf (argv[9], "%d", &m);
    sscanf (argv[10], "%lf", &eps);
    sscanf (argv[11], "%d", &maxit);
    sscanf (argv[12], "%d", &p);
}

void gui_data::realloc_data() {
    if (I) delete[] I;
    if (A) delete[] A;
    if (B) delete[] B;
    if (x) delete[] x;
    if (r) delete[] r;
    if (u) delete[] u;
    if (v) delete[] v;

    I = nullptr;
    A = nullptr;
    allocate_msr_matrix(nx, ny, &A, &I);
    init_reduce_sum(p);
    int N = (nx + 1) * (ny + 1);
    B = new double[N];
    x = new double[N];
    r = new double[N];
    u = new double[N];
    v = new double[N];

    fill_I(nx, ny, I);

    memset(x, 0, N*sizeof(double));
}

gui_data::~gui_data() {
    if (A) delete[] A;
    if (I) delete[] I;
    if (B) delete[] B;
    if (x) delete[] x;
    if (u) delete[] u;
    if (v) delete[] v;
    if (r) delete[] r;
}
