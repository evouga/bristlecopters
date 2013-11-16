#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "mesh.h"
#include <Eigen/Core>
#include <QObject>

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
    void simulate();
    void quit();
    void centerCamera();
    void updateGL();

private:    
    MainWindow &mw_;
    Mesh m_;
};

#endif // CONTROLLER_H
