#include "controller.h"
#include "mainwindow.h"
#include "mesh.h"
#include <string>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <QMessageBox>
#include <fstream>

using namespace std;
using namespace Eigen;

Controller::Controller(MainWindow &mw) : mw_(mw), m_(), time_(0)
{
    ProblemParameters params = m_.getParameters();
    mw_.setParameters(params);
    timer_ = new QTimer(this);
    connect(timer_, SIGNAL(timeout()), this, SLOT(tick()));
    timer_->start(100);
}

void Controller::tick()
{
    time_ += 0.1;
    updateGL();
}

void Controller::quit()
{
    QMetaObject::invokeMethod(&mw_, "close");
}

void Controller::renderMesh()
{
    m_.render(time_);
}

void Controller::exportOBJ(string filename)
{
    if(!m_.exportOBJ(filename.c_str()))
    {
        QString msg = "Couldn't write file " + QString::fromStdString(filename) + ". Save failed.";
        QMessageBox::warning(&mw_, "Couldn't Write File", msg, QMessageBox::Ok);
        return;
    }
}

void Controller::importOBJ(string filename)
{
    if(!m_.importOBJ(filename.c_str()))
    {
        string err = string("Couldn't load OBJ: ") + filename;
        QMetaObject::invokeMethod(&mw_, "showError", Q_ARG(std::string, err));
    }
    else
    {
        centerCamera();
    }
    double freq = m_.getModeFrequency();
    QMetaObject::invokeMethod(&mw_, "setModeFrequency", Q_ARG(double, freq));
    updateGL();
}

void Controller::updateParameters(ProblemParameters params)
{    
    m_.setParameters(params);
    double freq = m_.getModeFrequency();
    QMetaObject::invokeMethod(&mw_, "setModeFrequency", Q_ARG(double, freq));
    updateGL();
}

void Controller::findMode()
{
    m_.findMode();
    double freq = m_.getModeFrequency();
    QMetaObject::invokeMethod(&mw_, "setModeFrequency", Q_ARG(double, freq));
    updateGL();
}

void Controller::centerCamera()
{
    Vector3d centroid = m_.centroid();
    double radius = m_.radius();
    QMetaObject::invokeMethod(&mw_, "centerCamera", Q_ARG(Eigen::Vector3d, centroid), Q_ARG(double, radius));
}

void Controller::updateGL()
{
    QMetaObject::invokeMethod(&mw_, "repaintMesh");
}
