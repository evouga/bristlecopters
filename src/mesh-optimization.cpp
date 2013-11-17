#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <Eigen/Core>

typedef Eigen::Triplet<double> Tr;

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

void Mesh::buildExtendedFreeBoundaryLaplacian(const VectorXd &q, Eigen::SparseMatrix<double> &L)
{
    int numverts = mesh_->n_vertices();
    L.resize(3*numverts, 3*numverts);
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
            for(int j=0; j<3; j++)
                Lcoeffs.push_back(Tr(3*vi.handle().idx()+j, 3*nb.idx()+j, weight));
            totweight += weight;
        }
        for(int j=0; j<3; j++)
            Lcoeffs.push_back(Tr(3*vi.handle().idx()+j, 3*vi.handle().idx()+j, -totweight));
    }

    // PITA boundary

    // a1 = s1 cos theta1
    // a2 = s2 cos theta2
    // alt = s1 sin theta1

    // diff = (v2*a1+v1*a2)/s3 - v0

    // tot = diff/alt * s3

    // tot = v2*a1/alt + v1*a2/alt - v0*s3/alt
    // = v2*(s1 cos theta1)/(s1 sin theta1) + v1 * (s2 cos theta2)/(s2 sin theta2) - v0*s3/alt
    // = v2 cot theta1 + v1 cot theta2 - v0*(a1/alt + a2/alt)
    // = v2 cot theta1 + v1 cot theta2 - v0 (cot theta1+ cot theta2)
    // = (v2-v0) cot theta1 + (v1-v0) cot theta2

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
        ivp = q.segment<3>(3*intvert.idx());
        v1p = q.segment<3>(3*v1.idx());
        v2p = q.segment<3>(3*v2.idx());


        double costheta1 = ((ivp-v1p).dot(v2p-v1p));
        double sintheta1 = ((ivp-v1p).cross(v2p-v1p)).norm();
        double cottheta1 = costheta1/sintheta1;

        double costheta2 = ((ivp-v2p).dot(v1p-v2p));
        double sintheta2 = ((ivp-v2p).cross(v1p-v2p)).norm();
        double cottheta2 = costheta2/sintheta2;

        for(int j=0; j<3; j++)
        {
            Lcoeffs.push_back(Tr(3*v1.idx()+j, 3*v1.idx()+j, 0.5*cottheta2));
            Lcoeffs.push_back(Tr(3*v1.idx()+j, 3*v2.idx()+j, 0.5*cottheta1));
            Lcoeffs.push_back(Tr(3*v1.idx()+j, 3*intvert.idx()+j, -0.5*cottheta1-0.5*cottheta2));
            Lcoeffs.push_back(Tr(3*v2.idx()+j, 3*v1.idx()+j, 0.5*cottheta2));
            Lcoeffs.push_back(Tr(3*v2.idx()+j, 3*v2.idx()+j, 0.5*cottheta1));
            Lcoeffs.push_back(Tr(3*v2.idx()+j, 3*intvert.idx()+j, -0.5*cottheta1-0.5*cottheta2));
        }
    }

    L.setFromTriplets(Lcoeffs.begin(), Lcoeffs.end());
    SparseMatrix<double> Minv;
    buildExtendedInvMassMatrix(q, Minv);
    L = Minv*L;
}

bool Mesh::findMode(void)
{
    double coeff = params_.YoungsModulus*params_.h*params_.h*params_.h/(24.0*(1.0-params_.PoissonRatio*params_.PoissonRatio));
    VectorXd q(3*numverts());
    dofsFromGeometry(q);
    SparseMatrix<double> M(3*numverts(),3*numverts());
    buildExtendedMassMatrix(q, M);

    SparseMatrix<double> L(3*numverts(), 3*numverts());
    buildExtendedFreeBoundaryLaplacian(q, L);

    SparseMatrix<double> N(3*numverts(), numverts());
    vector<Tr> Ncoeffs;
    for(int i=0; i<numverts(); i++)
    {
        OMMesh::VertexHandle vh = mesh_->vertex_handle(i);
        OMMesh::Normal n;
        mesh_->calc_vertex_normal_correct(vh, n);
        Vector3d normal;
        for(int j=0; j<3; j++)
            normal[j] = n[j];
        double norm =normal.norm();
        if(fabs(norm) > 1e-8)
            normal /= norm;
        else
            normal.setZero();

        for(int j=0; j<3; j++)
            Ncoeffs.push_back(Tr(3*i+j, i, normal[j]));
    }
    N.setFromTriplets(Ncoeffs.begin(), Ncoeffs.end());

    SparseMatrix<double> left = params_.rho*params_.h*N.transpose()*M*N;
    SparseMatrix<double> right = coeff*N.transpose()*L.transpose()*M*L*N;

    std::cout << "solving" << std::endl;
    GeneralizedSelfAdjointEigenSolver<MatrixXd > solver(right, left);

    modes_ = N*solver.eigenvectors();
    modeFrequencies_ = solver.eigenvalues();
    for(int i=0; i<modeFrequencies_.size(); i++)
    {
        if(modeFrequencies_[i] < 0)
            modeFrequencies_[i] = 0;
        else
            modeFrequencies_[i] = sqrt(modeFrequencies_[i]);
    }

    for(int i=0; i<modes_.cols(); i++)
    {
        double maxnorm = 0;
        for(int j=0; j<numverts(); j++)
        {
            Vector3d vec = modes_.col(i).segment<3>(3*j);
            maxnorm = std::max(maxnorm, vec.norm());
        }
        if(maxnorm > 0)
            modes_.col(i) /= maxnorm;
    }
    modes_.col(0).setZero();

    return true;
}

void Mesh::buildExtendedMassMatrix(const VectorXd &q, Eigen::SparseMatrix<double> &M) const
{
    int numverts = mesh_->n_vertices();
    M.resize(3*numverts, 3*numverts);
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(q, vidx);
        for(int j=0; j<3; j++)
            entries.push_back(Tr(3*vidx+j, 3*vidx+j, area));
    }

    M.setFromTriplets(entries.begin(), entries.end());
}

void Mesh::buildExtendedInvMassMatrix(const VectorXd &q, Eigen::SparseMatrix<double> &M) const
{
    int numverts = mesh_->n_vertices();
    M.resize(3*numverts, 3*numverts);
    vector<Tr> entries;
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        int vidx = vi.handle().idx();
        double area = barycentricDualArea(q, vidx);
        for(int j=0; j<3; j++)
            entries.push_back(Tr(3*vidx+j, 3*vidx+j, 1.0/area));
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
