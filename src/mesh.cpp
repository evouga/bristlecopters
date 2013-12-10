#include <OpenMesh/Core/IO/MeshIO.hh>
#include "mesh.h"

using namespace Eigen;
using namespace OpenMesh;
using namespace std;

Mesh::Mesh() : meshLock_(QMutex::Recursive)
{
    params_.h = 0.001;
    params_.YoungsModulus = 1e7;
    params_.PoissonRatio = 0.5;
    params_.rho = 100;
    params_.scale = 1.0;
    params_.animAmplitude = 0.1;
    params_.animSpeed = 1.0;

    params_.smoothShade = true;
    params_.showWireframe = true;

    params_.curMode = 0;

    mesh_ = new OMMesh();
}

int Mesh::numverts() const
{
    return mesh_->n_vertices();
}

int Mesh::numdofs() const
{
    return 3*mesh_->n_vertices();
}

int Mesh::numedges() const
{
    return mesh_->n_edges();
}

void Mesh::dofsFromGeometry(Eigen::VectorXd &q, Eigen::VectorXd &g) const
{
    q.resize(numdofs());
    g.resize(numedges());

    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        for(int j=0; j<3; j++)
            q[3*i+j] = pt[j];
    }

    for(int i=0; i<(int)mesh_->n_edges(); i++)
    {
        g[i] = mesh_->data(mesh_->edge_handle(i)).restlen();
    }
}

void Mesh::dofsToGeometry(const VectorXd &q, const VectorXd &g)
{
    meshLock_.lock();
    {
        assert(q.size() == numdofs());
        assert(g.size() == numedges());

        for(int i=0; i<(int)mesh_->n_vertices(); i++)
        {
            OMMesh::Point &pt = mesh_->point(mesh_->vertex_handle(i));
            for(int j=0; j<3; j++)
                pt[j] = q[3*i+j];
        }

        for(int i=0; i<(int)mesh_->n_edges(); i++)
        {
            mesh_->data(mesh_->edge_handle(i)).setRestlen(g[i]);
        }
    }
    meshLock_.unlock();
}

void Mesh::setIntrinsicLengthsToCurrentLengths()
{
    meshLock_.lock();
    {
        for(OMMesh::EdgeIter ei = mesh_->edges_begin(); ei != mesh_->edges_end(); ++ei)
        {
            double length = mesh_->calc_edge_length(ei.handle());
            mesh_->data(ei.handle()).setRestlen(length);
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

void Mesh::edgeEndpointsWithMode(OMMesh::EdgeHandle edge, OMMesh::Point &p1, OMMesh::Point &p2, int mode, double amp)
{
    OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(edge,0);

    OMMesh::VertexHandle to = mesh_->to_vertex_handle(heh);
    OMMesh::VertexHandle from = mesh_->from_vertex_handle(heh);

    p1 = mesh_->point(from);
    p2 = mesh_->point(to);

    if(modes_.rows() == 3*mesh_->n_vertices() && mode >= 0 && modes_.cols() > mode)
    {
        for(int j=0; j<3; j++)
        {
            p1[j] += amp*modes_.col(mode)[3*from.idx()+j];
            p2[j] += amp*modes_.col(mode)[3*to.idx()+j];
        }
    }
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
        modeFrequencies_.resize(0);
        modes_.resize(0,0);
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

double Mesh::pointModeValue(OMMesh::VertexHandle vert, int mode)
{
    double result = 0;

    if(modes_.rows() == 3*mesh_->n_vertices() && mode >= 0 && mode < modes_.cols())
    {
        OMMesh::Normal n;
        mesh_->calc_vertex_normal_correct(vert, n);
        double norm=0;
        for(int i=0; i<3; i++)
        {
            norm += n[i]*n[i];
            result += n[i]*modes_.col(mode)[3*vert.idx()+i];
        }
        result /= sqrt(norm);
    }
    return result;
}

double Mesh::modeAmp(int mode, double time)
{
    double result = 0;

    if(modes_.rows() == 3*mesh_->n_vertices() && mode >= 0 && mode < modes_.cols())
    {
        result = sin(modeFrequencies_[mode]*time);
    }
    return result;
}

void Mesh::pointWithMode(OMMesh::VertexHandle vert, OMMesh::Point &pt, int mode, double amp)
{
    pt = mesh_->point(vert);

    if(modes_.rows() == 3*mesh_->n_vertices() && mode >= 0 && mode < modes_.cols())
    {
        for(int j=0; j<3; j++)
            pt[j] += amp*modes_.col(mode)[3*vert.idx()+j];
    }
}

double Mesh::getModeFrequency() const
{
    if(modes_.rows() == 3*mesh_->n_vertices() && params_.curMode >= 0 && params_.curMode < modes_.cols())
    {
        return modeFrequencies_[params_.curMode];
    }
    return 0;
}
