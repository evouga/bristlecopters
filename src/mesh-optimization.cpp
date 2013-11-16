#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "newton.h"

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

double Mesh::cotanWeight(int edgeid, const VectorXd &q)
{
    OMMesh::EdgeHandle eh = mesh_->edge_handle(edgeid);
    double weight = 0;
    for(int i=0; i<2; i++)
    {
        OMMesh::HalfedgeHandle heh = mesh_->halfedge_handle(eh,i);

        if(mesh_->is_boundary(heh))
            continue;
        OMMesh::HalfedgeHandle next = mesh_->next_halfedge_handle(heh);

        Vector3d e1, e2;
        OMMesh::VertexHandle oppv = mesh_->to_vertex_handle(next);
        OMMesh::VertexHandle v1 = mesh_->to_vertex_handle(heh);
        OMMesh::VertexHandle v2 = mesh_->from_vertex_handle(heh);
        e1 = q.segment<3>(3*v1.idx())-q.segment<3>(3*oppv.idx());
        e2 = q.segment<3>(3*v2.idx())-q.segment<3>(3*oppv.idx());
        double cosang = e1.dot(e2);
        double sinang = (e1.cross(e2)).norm();
        weight += 0.5*cosang/sinang;
    }
    return weight;
}

void Mesh::dirichletLaplacian(const VectorXd &q, Eigen::SparseMatrix<double> &L)
{
    int numverts = mesh_->n_vertices();
    L.resize(numverts, numverts);
    L.setZero();
    vector<Tr> Lcoeffs;

    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        double totweight = 0;

        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::EdgeHandle eh = mesh_->edge_handle(voh.handle());
            double weight = cotanWeight(eh.idx(), q);
            OMMesh::VertexHandle nb = mesh_->to_vertex_handle(voh.handle());
            Lcoeffs.push_back(Tr(vi.handle().idx(), nb.idx(), weight));
            totweight += weight;
        }
        Lcoeffs.push_back(Tr(vi.handle().idx(), vi.handle().idx(), -totweight));
    }
    L.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());
    SparseMatrix<double> Minv;
    buildInvMassMatrix(q, Minv);
    L = Minv*L;
}

void Mesh::neumannLaplacian(const VectorXd &q, Eigen::SparseMatrix<double> &L)
{
    int numverts = mesh_->n_vertices();
    L.resize(numverts, numverts);
    L.setZero();
    vector<Tr> Lcoeffs;

    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        double totweight = 0;

        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::EdgeHandle eh = mesh_->edge_handle(voh.handle());
            double weight = cotanWeight(eh.idx(), q);
            OMMesh::VertexHandle nb = mesh_->to_vertex_handle(voh.handle());
            Lcoeffs.push_back(Tr(vi.handle().idx(), nb.idx(), weight));
            totweight += weight;
        }
        Lcoeffs.push_back(Tr(vi.handle().idx(), vi.handle().idx(), -totweight));
    }

    // PITA boundary

    for(OMMesh::HalfedgeIter hi = mesh_->halfedges_begin(); hi != mesh_->halfedges_end(); ++hi)
    {
        if(!mesh_->is_boundary(hi.handle()))
            continue;

        OMMesh::HalfedgeHandle opp = mesh_->opposite_halfedge_handle(hi.handle());
        OMMesh::HalfedgeHandle next = mesh_->next_halfedge_handle(opp);

        OMMesh::VertexHandle intvert = mesh_->to_vertex_handle(next);
        OMMesh::VertexHandle v1 = mesh_->to_vertex_handle(opp);
        OMMesh::VertexHandle v2 = mesh_->from_vertex_handle(opp);

        Vector3d ivp, v1p, v2p;
        for(int j=0; j<3; j++)
        {
            ivp[j] = mesh_->point(intvert)[j];
            v1p[j] = mesh_->point(v1)[j];
            v2p[j] = mesh_->point(v2)[j];
        }

        double sintheta = ((v1p-ivp).cross(v2p-ivp)).norm();
        double costheta = ((v1p-ivp).dot(v2p-ivp));
        double magprod = (v1p-ivp).norm()*(v2p-ivp).norm();
        double weight = sintheta/(magprod+costheta);
        Lcoeffs.push_back(Tr(v1.idx(), v1.idx(), 0.5*weight));
        Lcoeffs.push_back(Tr(v1.idx(), v2.idx(), 0.5*weight));
        Lcoeffs.push_back(Tr(v1.idx(), intvert.idx(), -weight));
        Lcoeffs.push_back(Tr(v2.idx(), v1.idx(), 0.5*weight));
        Lcoeffs.push_back(Tr(v2.idx(), v2.idx(), 0.5*weight));
        Lcoeffs.push_back(Tr(v2.idx(), intvert.idx(), -weight));
    }

    L.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());
    SparseMatrix<double> Minv;
    buildInvMassMatrix(q, Minv);
    L = Minv*L;
}

class ImplicitEulerStep : public NewtonObjective
{
public:
    ImplicitEulerStep(Mesh &m,
                      Controller &cont,
                      const VectorXd &q0,
                      const SparseMatrix<double> &Ld,
                      const SparseMatrix<double> &Ln,
                      const SparseMatrix<double> &M,
                      const SparseMatrix<double> &Bdry,
                      double h, double kb, double ks, double A, double B) :
        m_(m), cont_(cont), q0_(q0), Ld_(Ld), Ln_(Ln), M_(M), Bdry_(Bdry), h_(h), kb_(kb), ks_(ks), A_(A), B_(B) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        int nverts = m_.numverts();
        SparseMatrix<double> hcube(nverts, nverts);
        vector<Tr> hcubecoeffs;

        VectorXd LJ(nverts);
        VectorXd LJ0(nverts);
        for(int j=0; j<nverts; j++)
        {
            hcubecoeffs.push_back(Tr(j,j,q[j]*q[j]*q[j]));
            LJ[j] = -3.0*A_/q[j]/q[j]/q[j]/q[j] + 9.0*B_/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j];
            LJ0[j] = -3.0*A_/q0_[j]/q0_[j]/q0_[j]/q0_[j] + 9.0*B_/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j];
        }
        hcube.setFromTriplets(hcubecoeffs.begin(),hcubecoeffs.end());

        SparseMatrix<double> bendingTerm = -Ln_.transpose()*M_*Ln_;

        VectorXd fsurface = kb_*bendingTerm*q + ks_*M_*Ln_*q + LJ;
        VectorXd fsurface0 = kb_*bendingTerm*q0_ + ks_* M_*Ln_*q0_ + LJ;

        SparseMatrix<double> id(nverts,nverts);
        id.setIdentity();

        VectorXd f = (q-q0_)/h_ + hcube*Ld_*fsurface;
        return 0.5*f.dot(f);
    }

    virtual double getEnergyAndDerivatives(const VectorXd &q, VectorXd &grad, SparseMatrix<double> &hess) const
    {
        int nverts = m_.numverts();
        SparseMatrix<double> id(nverts,nverts);
        id.setIdentity();

        SparseMatrix<double> hcube(nverts, nverts);
        vector<Tr> hcubecoeffs;

        VectorXd LJ(nverts);
        VectorXd LJ0(nverts);
        for(int j=0; j<nverts; j++)
        {
            hcubecoeffs.push_back(Tr(j,j,q[j]*q[j]*q[j]));
            LJ[j] = -3.0*A_/q[j]/q[j]/q[j]/q[j] + 9.0*B_/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j];
            LJ0[j] = -3.0*A_/q0_[j]/q0_[j]/q0_[j]/q0_[j] + 9.0*B_/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j]/q0_[j];
        }
        hcube.setFromTriplets(hcubecoeffs.begin(),hcubecoeffs.end());

        SparseMatrix<double> dLJ(nverts,nverts);
        vector<Tr> dLJcoeffs;
        for(int j=0; j<nverts; j++)
        {
            dLJcoeffs.push_back(Tr(j,j,12.0*A_/q[j]/q[j]/q[j]/q[j]/q[j] - 90.0*B_/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]/q[j]));
        }
        dLJ.setFromTriplets(dLJcoeffs.begin(), dLJcoeffs.end());

        SparseMatrix<double> bendingTerm = -Ln_.transpose()*M_*Ln_;

        VectorXd fsurface = kb_*bendingTerm*q + ks_* M_*Ln_*q + LJ;
        VectorXd fsurface0 = kb_*bendingTerm*q0_ + ks_* M_*Ln_*q0_ + LJ;
        VectorXd f = (q-q0_)/h_ + hcube*Ld_*fsurface;

        SparseMatrix<double> dfsurface(nverts,nverts);
        dfsurface = kb_*bendingTerm + ks_*M_*Ln_ + dLJ;
        SparseMatrix<double> df(nverts,nverts);
        df.setIdentity();
        df /= h_;
        df += hcube*Ld_*dfsurface;


        SparseMatrix<double> dhcube(nverts,nverts);
        vector<Tr> dhcubecoeffs;

        for(int j=0; j<nverts; j++)
        {
            double dhcfac = 3.0*q[j]*q[j];
            dhcubecoeffs.push_back(Tr(j,j,dhcfac));
        }
        dhcube.setFromTriplets(dhcubecoeffs.begin(), dhcubecoeffs.end());

        SparseMatrix<double> fsurfacediag(nverts,nverts);
        vector<Tr> fsurfacediagcoeffs;
        for(int j=0; j<nverts; j++)
        {
            fsurfacediagcoeffs.push_back(Tr(j,j,fsurface[j]));
        }
        fsurfacediag.setFromTriplets(fsurfacediagcoeffs.begin(), fsurfacediagcoeffs.end());

        df += dhcube*Ld_*fsurfacediag;

        grad = df.transpose()*f;
        hess = df.transpose()*df;

        return 0.5*f.dot(f);
    }

    virtual void showCurrentIteration(const VectorXd &q) const
    {
        m_.dofsToGeometry(q);
        cont_.updateGL();
    }

private:
    Mesh &m_;
    Controller &cont_;
    const VectorXd &q0_;
    const SparseMatrix<double> &Ld_;
    const SparseMatrix<double> &Ln_;
    const SparseMatrix<double> &M_;
    const SparseMatrix<double> &Bdry_;
    double h_;
    double kb_;
    double ks_;
    double A_;
    double B_;
};

bool Mesh::simulate(Controller &cont)
{
    double ts = params_.eulerTimestep;
    int numsteps = params_.numEulerIters;

    VectorXd fullq(3*numverts());
    for(int i=0; i<numverts(); i++)
    {
        OMMesh::Point pt = mesh_->point(mesh_->vertex_handle(i));
        for(int j=0; j<3; j++)
            fullq[3*i+j] = pt[j];
    }
    SparseMatrix<double> Ld;
    SparseMatrix<double> Ln;
    SparseMatrix<double> M;
    SparseMatrix<double> Bdry(numverts(), numverts());
    dirichletLaplacian(fullq,Ld);
    neumannLaplacian(fullq, Ln);
    buildMassMatrix(fullq, M);

    vector<Tr> bdrycoeffs;
    for(int i=0; i<numverts(); i++)
    {
        if(mesh_->is_boundary(mesh_->vertex_handle(i)))
            bdrycoeffs.push_back(Tr(i,i,1));
    }
    Bdry.setFromTriplets(bdrycoeffs.begin(), bdrycoeffs.end());

    setBump();
    cont.updateGL();
    VectorXd h(numverts());
    dofsFromGeometry(h);

    double kb = 1e-5;
    double ks = 0*1;
    double A = .01;
    double B = A/30.0;//*1e-5;

    for(int iter = 0; iter < numsteps; iter++)
    {
        ImplicitEulerStep *obj = new ImplicitEulerStep(*this, cont, h, Ld, Ln, M, Bdry, ts, kb, ks, A, B);
        Newton n(*obj);
        NewtonParameters params;
        params.tol = params_.tol;
        params.maxiters = params_.maxiters;
        params.lsmaxiters = params_.maxlinesearchiters;
        //params.descentcheck = false;
        VectorXd result;
        Newton::SolverStatus ss = n.solve(params, h, result);
        std::cout << n.solverStatusMessage(ss) << std::endl;

        delete obj;

        if(ss != Newton::CONVERGED)
            return false;

        h = result;
        dofsToGeometry(h);
        cont.updateGL();
    }
    return true;
}

void Mesh::setBump()
{
    for(int i=0; i<(int)mesh_->n_vertices(); i++)
    {
        OMMesh::VertexHandle vh = mesh_->vertex_handle(i);
        OMMesh::Point pos = mesh_->point(vh);
        double r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
        double z = 1.1 - 0.1*tanh(r)*tanh(r);
        mesh_->point(vh)[2] = z;
    }
}

void Mesh::buildMassMatrix(const VectorXd &q, Eigen::SparseMatrix<double> &M) const
{
    int numverts = mesh_->n_vertices();
    M.resize(numverts, numverts);
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(q, vidx);
        entries.push_back(Tr(vidx, vidx, area));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildInvMassMatrix(const VectorXd &q, Eigen::SparseMatrix<double> &M) const
{
    int numverts = mesh_->n_vertices();
    M.resize(numverts, numverts);
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(q, vidx);
        entries.push_back(Tr(vidx, vidx, 1.0/area));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

double Mesh::barycentricDualArea(const VectorXd &q, int vidx) const
{
    double result = 0;
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    for(OMMesh::VertexFaceIter vfi = mesh_->vf_iter(vh); vfi; ++vfi)
    {
        result += faceArea(q, vfi.handle().idx());
    }
    return result/3.0;
}

double Mesh::faceArea(const VectorXd &q, int fidx) const
{
    FaceHandle fh = mesh_->face_handle(fidx);
    int verts[3];
    int idx=0;
    for(OMMesh::FaceVertexIter fvi = mesh_->fv_iter(fh); fvi; ++fvi)
    {
        verts[idx++] = fvi.handle().idx();
    }

    Vector3d q0 = q.segment<3>(3*verts[0]);
    Vector3d q1 = q.segment<3>(3*verts[1]);
    Vector3d q2 = q.segment<3>(3*verts[2]);

    double A = ((q1-q0).cross(q2-q0)).norm();
    return 0.5*A;
}
