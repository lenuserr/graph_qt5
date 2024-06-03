#ifndef WINDOW_H
#define WINDOW_H
#include <QtWidgets/QtWidgets>
#include "msr.h"
#include "gui_data.h"

class Window : public QWidget {
  Q_OBJECT

private:
  QWidget *widget;
  int func_id;
  const char *f_name;
  const char* program_name;
  int n_graph = 0;
  double f_max;
  int p_ = 0;

  Args* args = nullptr;
  pthread_t* tid = nullptr;
  gui_data data;
  const int task = 3;
  bool threads_created = false;
  pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

public:
  Window (QWidget *parent);

  ~Window() {
    if (args) delete[] args;
    if (tid) delete[] tid;
  }

  QSize minimumSizeHint () const;
  QSize sizeHint () const;

  int parse_command_line (int argc, char *argv[]);
  void select_f();
  void draw_triangle(QPainter* painter, point2d p1, point2d p2, point2d p3, QColor color);
  void draw_f(QPainter* painter);
  void draw_Pf(QPainter* painter);
  void draw_residual(QPainter* painter);
  void draw_text(QPainter* painter);
  bool msr_ready();
  void waiting_msr();
  QPointF l2g (double x_loc, double y_loc, double y_min, double y_max);
public slots:
  void change_func();
  void previous_func();
  void change_graph();
  void increase_s();
  void decrease_s();
  void increase_n();
  void decrease_n();
  void increase_p();
  void decrease_p();
  void increase_m();
  void decrease_m();
  void close_window();


protected:
  void paintEvent (QPaintEvent *event);
};

#endif
