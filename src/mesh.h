#ifndef MESH_H
#define MESH_H

#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <QMutex>

class Controller;

struct MyTraits : public OpenMesh::DefaultTraits
{
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> OMMesh;

struct ProblemParameters
{
    // rendering
    bool showWireframe;
    bool smoothShade;
    int curMode;

    // material
    double h;
    double PoissonRatio;
    double YoungsModulus;
    double rho;
};

class Mesh
{
public:
    Mesh();

    bool findMode();

    int numedges() const;
    int numverts() const;
    const ProblemParameters &getParameters() const;
    void setParameters(ProblemParameters params);
    double getModeFrequency() const;

    // Rendering methods. These run concurrently and must all lock the meshLock before reading from the mesh.
    void render(double t);
    Eigen::Vector3d centroid();
    double radius();
    // End rendering methods

    bool exportOBJ(const char *filename);
    bool importOBJ(const char *filename);

private:
    void dofsFromGeometry(Eigen::VectorXd &q) const;
    void dofsToGeometry(const Eigen::VectorXd &q);

    double cotanWeight(int eidx, const Eigen::VectorXd &q);
    void buildExtendedFreeBoundaryLaplacian(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &L);

    void edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2);
    void edgeEndpointsWithMode(OMMesh::EdgeHandle edge, OMMesh::Point &p1, OMMesh::Point &p2, int mode, double amp);
    void buildExtendedMassMatrix(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &M) const;
    void buildExtendedInvMassMatrix(const Eigen::VectorXd &q, Eigen::SparseMatrix<double> &M) const;
    double barycentricDualArea(const Eigen::VectorXd &q, int vidx) const;
    double faceArea(const Eigen::VectorXd &q, int fidx) const;

    double pointModeValue(OMMesh::VertexHandle vert, int mode);
    double modeAmp(int mode, double time);
    void pointWithMode(OMMesh::VertexHandle vert, OMMesh::Point &pt, int mode, double amp);

    Eigen::Vector3d colormap(double val) const;
    Eigen::Vector3d colormap(double val, double min, double max) const;
    Eigen::Vector3d HSLtoRGB(const Eigen::Vector3d &hsl) const;

    OMMesh *mesh_;
    ProblemParameters params_;

    Eigen::MatrixXd modes_;
    Eigen::VectorXd modeFrequencies_;

    // The rendering thread reads the mesh and its edge data. Any function must lock this before writing to
    // to the mesh. (The rendering thread does not write to the mesh so reads from the worker thread do not
    // need to lock.)
    QMutex meshLock_;
};

#endif // MESH_H
