#pragma once

struct point2d {
    double x;
    double y;

    point2d() : x(0), y(0)
        {}
    point2d(double _x, double _y) : x(_x), y(_y)
        {}
};


struct gui_data {
    double *A = nullptr;
    int *I = nullptr;
    double *B = nullptr;
    double *x = nullptr;
    double *u = nullptr;
    double *v = nullptr;
    double *r = nullptr;

    double a;
    double b;
    double c;
    double d;
    double eps;
    int nx;
    int ny;
    int mx;
    int my;
    int m;
    int maxit;
    int p;
    double (*f)(double, double);

    void read_data(char* argv[]);
    void realloc_data();
    double f_min_max(double& abs_min, double& abs_max);
    void pf_min_max(double& abs_min, double& abs_max);
    void residual_min_max(double& abs_min, double& abs_max);
    double pf(point2d p);
    ~gui_data();
};
