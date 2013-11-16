#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"

using namespace Eigen;
using namespace OpenMesh;
using namespace std;

Mesh::Mesh() : meshLock_(QMutex::Recursive)
{
    params_.h = 1;
    params_.YoungsModulus = 1;
    params_.PoissonRatio = 0.5;
    params_.maxiters = 50;
    params_.maxlinesearchiters = 100;
    params_.tol = 1e-6;
    params_.rho = 1.0;
    params_.dampingCoeff = 0.9;
    params_.eulerTimestep = 0.1;
    params_.numEulerIters = 100;

    params_.smoothShade = true;
    params_.showWireframe = true;

    mesh_ = new OMMesh();
}

int Mesh::numverts() const
{
    return mesh_->n_vertices();
}

int Mesh::numedges() const
{
    return mesh_->n_edges();
}

void Mesh::dofsFromGeometry(Eigen::VectorXd &h) const
{
    h.resize(numverts());

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        h[i] = pt[2];
    }
}

void Mesh::dofsToGeometry(const VectorXd &h)
{    
    meshLock_.lock();
    {
        assert(h.size() == numverts());

        for(int i=0; i<(int)mesh_->n_vertices(); i++)
        {
            OMMesh::Point &pt = mesh_->point(mesh_->vertex_handle(i));
            pt[2] = h[i];
        }
    }
    meshLock_.unlock();
}

void Mesh::edgeEndpoints(OMMesh::EdgeHandle eh, OMMesh::Point &pt1, OMMesh::Point &pt2)
{
    OMMesh::HalfedgeHandle heh1 = mesh_->halfedge_handle(eh, 0);
    pt1 = mesh_->point(mesh_->from_vertex_handle(heh1));
    pt2 = mesh_->point(mesh_->to_vertex_handle(heh1));
}

bool Mesh::exportOBJ(const char *filename)
{
    OpenMesh::IO::Options opt;
    mesh_->request_face_normals();
    mesh_->request_vertex_normals();
    mesh_->update_normals();
    opt.set(OpenMesh::IO::Options::VertexNormal);
    return OpenMesh::IO::write_mesh(*mesh_, filename, opt);
}

bool Mesh::importOBJ(const char *filename)
{
    bool success = true;
    meshLock_.lock();
    {
        OpenMesh::IO::Options opt;
        mesh_->request_face_normals();
        mesh_->request_vertex_normals();
        opt.set(OpenMesh::IO::Options::VertexNormal);
        success = OpenMesh::IO::read_mesh(*mesh_, filename, opt);
        mesh_->update_normals();
    }
    meshLock_.unlock();
    return success;
}

const ProblemParameters &Mesh::getParameters() const
{
    return params_;
}

void Mesh::setParameters(ProblemParameters params)
{
    params_ = params;
}
