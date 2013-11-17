#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "mesh.h"
#include <Eigen/Core>
#include <QObject>
#include <QTimer>

class MainWindow;

class Controller : public QObject
{
    Q_OBJECT

public:
    Controller(MainWindow &mw);

    void renderMesh();

public slots:
    void exportOBJ(std::string filename);
    void importOBJ(std::string filename);
    void updateParameters(ProblemParameters params);
    void findMode();
    void quit();
    void centerCamera();
    void updateGL();
    void tick();

private:    
    MainWindow &mw_;
    Mesh m_;
    QTimer *timer_;
    double time_;
};

#endif // CONTROLLER_H
