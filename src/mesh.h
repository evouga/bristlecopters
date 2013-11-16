#ifndef MESH_H
#define MESH_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "elasticenergy.h"
#include <QMutex>

class Controller;

struct MyTraits : public OpenMesh::DefaultTraits
{
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> OMMesh;

struct ProblemParameters : public ElasticParameters
{
    // simulation
    int maxiters;
    int maxlinesearchiters;
    double tol;
    double rho;

    double eulerTimestep;
    double dampingCoeff;
    int numEulerIters;

    // rendering
    bool showWireframe;
    bool smoothShade;
};

class Mesh
{
public:
    Mesh();

    bool simulate(Controller &cont);

    int numedges() const;
    int numverts() const;
    const ProblemParameters &getParameters() const;
    void setParameters(ProblemParameters params);

    // Rendering methods. These run concurrently and must all lock the meshLock before reading from the mesh.
    void render();
    Eigen::Vector3d centroid();
    double radius();
    // End rendering methods

    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

    friend class EmbeddingMinimizer;
    friend class MetricFit;
    friend class ImplicitEulerStep;

private:
    void dofsFromGeometry(Eigen::VectorXd &h) const;
    void dofsToGeometry(const Eigen::VectorXd &h);

    double cotanWeight(int eidx, const Eigen::VectorXd &q);
    void dirichletLaplacian(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &L);
    void neumannLaplacian(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &L);
    void setBump();

    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    void buildMassMatrix(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &M) const;
    void buildInvMassMatrix(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &M) const;
    double barycentricDualArea(const Eigen::VectorXd &q, int vidx) const;
    double faceArea(const Eigen::VectorXd &q, int fidx) const;

    Eigen::Vector3d colormap(double val) const;
    Eigen::Vector3d colormap(double val, double max) const;
    Eigen::Vector3d HSLtoRGB(const Eigen::Vector3d &hsl) const;

    OMMesh *mesh_;
    ProblemParameters params_;

    // The rendering thread reads the mesh and its edge data. Any function must lock this before writing to
    // to the mesh. (The rendering thread does not write to the mesh so reads from the worker thread do not
    // need to lock.)
    QMutex meshLock_;
};

#endif // MESH_H
