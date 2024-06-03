#include <QPainter>
#include <stdio.h>
#include "window.h"

#define L2G(X,Y) (l2g ((X), (Y), min_y, max_y))
#define EPS 1e-14
#define BUF 1024

Window::Window (QWidget *parent)
  : QWidget (parent) {

    widget = parent;
}

QSize Window::minimumSizeHint () const {
  return QSize (100, 100);
}

QSize Window::sizeHint () const {
  return QSize (1000, 1000);
}

bool Window::msr_ready() {
  for (int i = 0; i < data.p; i++) {
      if (!args[i].ready) {
          return false;
      }
  }

  return true;
}

void Window::waiting_msr() {
    if (msr_ready()) {
        printf (
          "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\n\
                It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
                program_name, task, args[0].res_1, args[0].res_2, args[0].res_3, args[0].res_4,
                args[0].t1, args[0].t2, args[0].its, data.eps, func_id, data.nx, data.ny, data.p);
        update();
    }
    else {
        QTimer::singleShot(200, this, &Window::waiting_msr);
    }
}

void Window::close_window() {
  if (msr_ready()) {
      widget->close();
  } else {
      QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
  }
}

int Window::parse_command_line (int argc, char* argv[]) {
    if (argc != 13) {
        return -1;
    }

    data.read_data(argv);
    data.realloc_data();

    args = new Args[data.p];
    tid = new pthread_t[data.p];
    program_name = argv[0];
    func_id = data.m - 1;

    for (int i = 0; i < data.p; i++) {
        args[i].copy_data(data);

        args[i].k = i;
        args[i].mutex = &mutex;
        args[i].cond = &cond;
    }

    change_func();

    return 0;
}

void Window::select_f() {
  double (*f)(double, double) = data.f;
  switch (func_id) {
      case 0:
        f_name = "f (x, y) = 1";
        f = f_0;
        break;
      case 1:
        f_name = "f (x, y) = x";
        f = f_1;
        break;
      case 2:
        f_name = "f (x, y) = y";
        f = f_2;
        break;
      case 3:
        f_name = "f (x, y) = x + y";
        f = f_3;
        break;
      case 4:
        f_name = "f (x, y) = sqrt(x*x + y*y)";
        f = f_4;
        break;
      case 5:
        f_name = "f (x, y) = x*x + y*y";
        f = f_5;
        break;
      case 6:
        f_name = "f (x, y) = exp(x*x - y*y)";
        f = f_6;
        break;
      case 7:
        f_name = "f (x, y) = 1 / (25(x*x + y*y) + 1)";
        f = f_7;
        break;
  }

  data.f = f;
}

void Window::change_func() {
    if (msr_ready()) {
        func_id = (func_id + 1) % 8;
        select_f();

        double abs_min, abs_max;
        f_max = data.f_min_max(abs_min, abs_max);

        data.realloc_data();

        for (int i = 0; i < data.p; i++) {
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        if (threads_created) {
            pthread_cond_broadcast(&cond);
        } else {
            threads_created = true;
            for (int i = 0; i < data.p; i++) {
                pthread_create(&args[i].tid, 0, solution, args + i);
            }
        }

        QTimer::singleShot(200, this, &Window::waiting_msr);

    } else {
        QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
    }
}

void Window::previous_func() {
    if (msr_ready()) {
        if (func_id) {
            func_id = (func_id - 1) % 8;
        } else {
            func_id = 7;
        }
        select_f();

        double abs_min, abs_max;
        f_max = data.f_min_max(abs_min, abs_max);

        data.realloc_data();

        for (int i = 0; i < data.p; i++) {
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond);
        QTimer::singleShot(200, this, &Window::waiting_msr);

    } else {
        QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
    }
}

void Window::increase_n() {
    if (msr_ready()) {
        data.nx *= 2;
        data.ny *= 2;
        data.realloc_data();

        for (int i = 0; i < data.p; i++) {
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond); // wake up
        QTimer::singleShot(200, this, &Window::waiting_msr);

    } else {
        QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
    }
}

void Window::decrease_n() {
    if (msr_ready()) {
        data.nx /= 2;
        data.ny /= 2;
        data.realloc_data();

        for (int i = 0; i < data.p; i++) {
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond); // wake up
        QTimer::singleShot(200, this, &Window::waiting_msr);

    } else {
        QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
    }
}

void Window::change_graph() {
    n_graph = (n_graph + 1) % 3;
    update();
}

void Window::increase_p() {
    if (msr_ready()) {
        p_ += 1;
        data.realloc_data();

        for (int i = 0; i < data.p; i++) {
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond); // wake up
        QTimer::singleShot(200, this, &Window::waiting_msr);

    } else {
        QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
    }
}

void Window::decrease_p() {
    if (msr_ready()) {
        p_ -= 1;
        data.realloc_data();

        for (int i = 0; i < data.p; i++) {
            args[i].copy_data(data);
            args[i].k = i;
            args[i].mutex = &mutex;
            args[i].cond = &cond;
            args[i].ready = false;
        }

        pthread_cond_broadcast(&cond); // wake up
        QTimer::singleShot(200, this, &Window::waiting_msr);

    } else {
        QMessageBox::warning(0, "Waiting for MSR!", "Calculations are not completed!");
    }
}

void Window::increase_s() {
    double len_x = data.b - data.a;
    double len_y = data.d - data.c;

    data.a += len_x / 4;
    data.b -= len_x / 4;
    data.c += len_y / 4;
    data.d -= len_y / 4;

    update();
}

void Window::decrease_s() {
    double len_x = data.b - data.a;
    double len_y = data.d - data.c;

    data.a -= len_x / 2;
    data.b += len_x / 2;
    data.c -= len_y / 2;
    data.d += len_y / 2;

    update();
}

void Window::increase_m() {
    data.mx *= 2;
    data.my *= 2;

    double abs_min, abs_max;
    f_max = data.f_min_max(abs_min, abs_max);

    update();
}

void Window::decrease_m() {
    data.mx /= 2;
    data.my /= 2;

    double abs_min, abs_max;
    f_max = data.f_min_max(abs_min, abs_max);

    update();
}

double get_normalize(double value, double abs_min, double abs_max) {
    if (value - abs_min < 0) return 0;
    if (abs_max - value < 0) return 1;
    if (std::fabs(abs_max - abs_min) < EPS) return 1;
    return std::fabs((value - abs_min) / (abs_max - abs_min));
}

void Window::draw_triangle(QPainter* painter, point2d p1, point2d p2, point2d p3, QColor color) {
    QPainterPath path;
    double min_y = data.c, max_y = data.d;

    path.moveTo(L2G(p1.x, p1.y));
    path.lineTo(L2G(p2.x, p2.y));
    path.lineTo(L2G(p3.x, p3.y));
    path.lineTo(L2G(p1.x, p1.y));

    painter->setPen (Qt::NoPen);
    painter->fillPath (path, QBrush(color));
}

void f_rgb(int& R, int& G, int& B, double coef) {
    R = coef * 255;
    G = coef * 0;
    B = coef * 0;
}

void pf_rgb(int& R, int& G, int& B, double coef) {
    R = coef * 0;
    G = coef * 255;
    B = coef * 0;
}

void residual_rgb(int& R, int& G, int& B, double coef) {
    R = coef * 0;
    G = coef * 0;
    B = coef * 255;
}

double gui_data::f_min_max(double& abs_min, double& abs_max) {
    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value1, value2;
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            value1 = f(a + (i + 1./3)*hx, c + (j + 2./3) * hy);
            value2 = f(a + (i + 2./3)*hx, c + (j + 1./3) * hy);
            abs_min = std::min(abs_min, value1);
            abs_max = std::max(abs_max, value1);
            abs_min = std::min(abs_min, value2);
            abs_max = std::max(abs_max, value2);
        }
    }

    return std::max(abs_min, abs_max);
}

void gui_data::pf_min_max(double& abs_min, double& abs_max) {
    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value1, value2;
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            value1 = pf({a + (i + 1./3)*hx, c + (j + 2./3) * hy});
            value2 = pf({a + (i + 2./3)*hx, c + (j + 1./3) * hy});
            abs_min = std::min(abs_min, value1);
            abs_max = std::max(abs_max, value1);
            abs_min = std::min(abs_min, value2);
            abs_max = std::max(abs_max, value2);
        }
    }
}

void gui_data::residual_min_max(double& abs_min, double& abs_max) {
    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double f_val, pf_val, value1, value2;
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {

            f_val = f(a + (i + 1./3)*hx, c + (j + 2./3) * hy);
            pf_val = pf({a + (i + 1./3)*hx, c + (j + 2./3) * hy});
            value1 = std::fabs(f_val - pf_val);

            f_val = f(a + (i + 2./3)*hx, c + (j + 1./3) * hy);
            pf_val = pf({a + (i + 2./3)*hx, c + (j + 1./3) * hy});
            value2 = std::fabs(f_val - pf_val);

            abs_min = std::min(abs_min, value1);
            abs_max = std::max(abs_max, value1);
            abs_min = std::min(abs_min, value2);
            abs_max = std::max(abs_max, value2);
        }
    }
}

void Window::draw_f(QPainter* painter) {
    double a = data.a; double b = data.b; double c = data.c; double d = data.d;
    double mx = data.mx; double my = data.my; double (*f)(double, double) = data.f;

    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value, coef;
    int R, G, B;
    point2d p1, p2, p3, p;
    double abs_min = f(a + hx/3, c + (2./3) * hy);
    double abs_max = abs_min;
    data.f_min_max(abs_min, abs_max);
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            // верхний треугольник в квадрате
            p1 = {a + i*hx, c + j*hy};
            p2 = {p1.x, c + (j + 1)*hy};
            p3 = {a + (i + 1)*hx, p2.y};

            value = f(a + (i + 1./3)*hx, c + (j + 2./3) * hy);
            coef = get_normalize(value, abs_min, abs_max);
            f_rgb(R, G, B, coef);
            draw_triangle(painter, p1, p2, p3, QColor(R, G, B));

            // нижний треугольник в квадрате
            p2 = {a + (i + 1)*hx, c + j*hy};
            value = f(a + (i + 2./3)*hx, c + (j + 1./3) * hy);
            coef = get_normalize(value, abs_min, abs_max);
            f_rgb(R, G, B, coef);
            draw_triangle(painter, p1, p2, p3, QColor(R, G, B));
        }
    }

    char max_[BUF];
    double f_max = std::max(std::fabs(abs_min), std::fabs(abs_max));
    sprintf(max_, "max|f| = %.2e", f_max);
    painter->setPen("white");
    painter->drawText(5, 40, max_);
    draw_text(painter);
}

void Window::draw_Pf(QPainter* painter) {
    double a = data.a; double b = data.b; double c = data.c; double d = data.d;
    double mx = data.mx; double my = data.my;

    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value, coef;
    int R, G, B;
    point2d p1, p2, p3, p;
    double abs_min = data.pf({a + hx/3, c + (2./3) * hy});
    double abs_max = abs_min;

    data.pf_min_max(abs_min, abs_max);
    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            // верхний треугольник в квадрате
            p1 = {a + i*hx, c + j*hy};
            p2 = {p1.x, c + (j + 1)*hy};
            p3 = {a + (i + 1)*hx, p2.y};

            value = data.pf({a + (i + 1./3)*hx, c + (j + 2./3) * hy});
            coef = get_normalize(value, abs_min, abs_max);
            pf_rgb(R, G, B, coef);
            draw_triangle(painter, p1, p2, p3, QColor(R, G, B));

            // нижний треугольник в квадрате
            p2 = {a + (i + 1)*hx, c + j*hy};

            value = data.pf({a + (i + 2./3)*hx, c + (j + 1./3) * hy});
            coef = get_normalize(value, abs_min, abs_max);
            pf_rgb(R, G, B, coef);
            draw_triangle(painter, p1, p2, p3, QColor(R, G, B));
        }
    }

    char max_[BUF];
    double pf_max = std::max(std::fabs(abs_min), std::fabs(abs_max));
    sprintf(max_, "max|Pf| = %.2e", pf_max);
    painter->setPen("blue");
    painter->drawText(5, 40, max_);
    draw_text(painter);
}

double length(double x, double y) {
    return sqrt(x*x + y*y);
}

double gui_data::pf(point2d p) {
    int i, j, l;
    double hnx = (b - a) / nx;
    double hny = (d - c) / ny;
    double w1, w2;

    i = (p.x - a) / hnx;
    j = (p.y - c) / hny;
    ij2l(nx, ny, i, j, l);

    if (std::fabs((p.x - a) - i * hnx) < EPS && fabs((p.y - c) - j * hny) < EPS) {
        return x[l];
    }

    double d2 = length((p.x - a) - (i + 1) * hnx, (p.y - c) - j*hny);
    double d3 = length((p.x - a) - i * hnx, (p.y - c) - (j + 1)*hny);

    if (d2 < d3) { // d2 vs d3
        w1 = ((i + 1) * hnx - (p.x - a)) / hnx;
        w2 = (p.y - c - j * hny) / hny;
        return  w1 * x[l] + w2 * x[l + 1 + (nx + 1)] + (1 - w1 - w2) * x[l + 1];
    }

    w1 = (p.x - a - i * hnx) / hnx;
    w2 = ((j + 1) * hny - (p.y - c)) / hny;
    return w1 * x[l + 1 + (nx + 1)] + w2 * x[l] + (1 - w1 - w2) * x[l + (nx + 1)];
}

void Window::draw_residual(QPainter *painter) {
    double a = data.a; double b = data.b; double c = data.c; double d = data.d;
    double mx = data.mx; double my = data.my; double (*f)(double, double) = data.f;

    double hx = (b - a) / mx;
    double hy = (d - c) / my;
    double value, f_val, pf_val, coef;
    int R, G, B;
    point2d p1, p2, p3, p;
    double abs_min = std::fabs(f(a + hx/3, c + (2./3) * hy) - data.pf({a + hx/3, c + (2./3) * hy}));
    double abs_max = abs_min;
    data.residual_min_max(abs_min, abs_max);

    for (int i = 0; i < mx; ++i) {
        for (int j = 0; j < my; ++j) {
            // верхний треугольник в квадрате
            p1 = {a + i*hx, c + j*hy};
            p2 = {p1.x, c + (j + 1)*hy};
            p3 = {a + (i + 1)*hx, p2.y};

            f_val = f(a + (i + 1./3)*hx, c + (j + 2./3) * hy);
            pf_val = data.pf({a + (i + 1./3)*hx, c + (j + 2./3) * hy});

            value = std::fabs(f_val - pf_val);
            coef = get_normalize(value, abs_min, abs_max);
            residual_rgb(R, G, B, coef);
            draw_triangle(painter, p1, p2, p3, QColor(R, G, B));

            // нижний треугольник в квадрате
            p2 = {a + (i + 1)*hx, c + j*hy};
            f_val = f(a + (i + 2./3)*hx, c + (j + 1./3) * hy);
            pf_val = data.pf({a + (i + 2./3)*hx, c + (j + 1./3) * hy});

            value = std::fabs(f_val - pf_val);
            coef = get_normalize(value, abs_min, abs_max);
            residual_rgb(R, G, B, coef);
            draw_triangle(painter, p1, p2, p3, QColor(R, G, B));
        }
    }

    char max_[BUF];
    double res_max = std::max(std::fabs(abs_min), std::fabs(abs_max));
    sprintf(max_, "max|Pf - f| = %.2e", res_max);
    painter->setPen("yellow");
    painter->drawText(5, 40, max_);
    draw_text(painter);
}

QPointF Window::l2g (double x_loc, double y_loc, double y_min, double y_max) {
  double x_gl = (x_loc - data.a) / (data.b - data.a) * width ();
  double y_gl = (y_max - y_loc) / (y_max - y_min) * height ();
  return QPointF (x_gl, y_gl);
}

void Window::draw_text(QPainter* painter) {
    painter->drawText (5, 20, f_name);
    std::string text1 = "a = " + std::to_string(data.a) + " b = " + std::to_string(data.b);
    std::string text2 = "c = " + std::to_string(data.c) + " d = " + std::to_string(data.d);
    std::string text3 = "nx = " + std::to_string(data.nx) + " ny = " + std::to_string(data.ny);
    std::string text4 = "mx = " + std::to_string(data.mx) + " my = " + std::to_string(data.my);
    std::string text5 = "p = " + std::to_string(p_);
    painter->drawText (5, 60, text1.c_str());
    painter->drawText (5, 80, text2.c_str());
    painter->drawText (5, 100, text3.c_str());
    painter->drawText (5, 120, text4.c_str());
    painter->drawText(5, 140, text5.c_str());
}

void Window::paintEvent (QPaintEvent * /* event */) {
  QPainter painter (this);

  if (n_graph == 0) {
      draw_f(&painter);
  } else if (n_graph == 1) {
      draw_Pf(&painter);
  } else {
      draw_residual(&painter);
  }
}
