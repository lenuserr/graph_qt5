#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>

#include "window.h"

int main (int argc, char *argv[])
{
  QApplication app (argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window);
  QAction *action;

  if (graph_area->parse_command_line (argc, argv)) {
      QMessageBox::warning (0, "Wrong input arguments!", 
                            "Wrong input arguments!");
      return -1;
  }

  action = tool_bar->addAction ("&Change function", graph_area, SLOT (change_func()));
  action->setShortcut (QString ("0"));

  action = tool_bar->addAction ("&Change graph", graph_area, SLOT (change_graph()));
  action->setShortcut (QString ("1"));

  action = tool_bar->addAction ("&Increase scale", graph_area, SLOT (increase_s()));
  action->setShortcut (QString ("2"));

  action = tool_bar->addAction ("&Decrease scale", graph_area, SLOT (decrease_s()));
  action->setShortcut (QString ("3"));

  action = tool_bar->addAction ("&Increase n", graph_area, SLOT (increase_n()));
  action->setShortcut (QString ("4"));

  action = tool_bar->addAction ("&Decrease n", graph_area, SLOT (decrease_n()));
  action->setShortcut (QString ("5"));

  action = tool_bar->addAction ("&Increase p", graph_area, SLOT (increase_p()));
  action->setShortcut (QString ("6"));

  action = tool_bar->addAction ("&Decrease p", graph_area, SLOT (decrease_p()));
  action->setShortcut (QString ("7"));

  action = tool_bar->addAction ("&Increase m", graph_area, SLOT (increase_m()));
  action->setShortcut (QString ("8"));

  action = tool_bar->addAction ("&Decrease m", graph_area, SLOT (decrease_m()));
  action->setShortcut (QString ("9"));

  action = tool_bar->addAction ("&Previous function", graph_area, SLOT (previous_func()));
  action->setShortcut (QString ("-"));

  action = tool_bar->addAction("E&xit", graph_area, SLOT(close_window()));
  action->setShortcut(QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");

  window->show ();
  app.exec ();

  graph_area->close_window();

  delete window;
  return 0;
}
