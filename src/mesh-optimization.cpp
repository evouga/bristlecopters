#include "mesh.h"
#include <iomanip>
#include "controller.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "elasticenergy.h"

typedef Eigen::Triplet<double> Tr;

using namespace std;
using namespace Eigen;
using namespace OpenMesh;

void Mesh::elasticEnergy(const VectorXd &q,
                         const VectorXd &g,
                         double &energyB,
                         double &energyS,
                         VectorXd &gradq,
                         Eigen::SparseMatrix<double> &hessq,
                         Eigen::SparseMatrix<double> &gradggradq,
                         int derivativesRequested) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energyB = energyS = 0;

    if(derivativesRequested & ElasticEnergy::DR_DQ)
    {
        gradq.resize(numdofs());
        gradq.setZero();
    }

    if(derivativesRequested & ElasticEnergy::DR_HQ)
        hessq.resize(numdofs(), numdofs());

    if(derivativesRequested & ElasticEnergy::DR_DGDQ)
    {
        gradggradq.resize(numedges(), numdofs());
    }


    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;
    vector<Tr> dgdqcoeffs;

    // bending energy
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        vector<int> spokeidx;
        vector<int> rightoppidx;
        vector<int> nbidx;
        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::HalfedgeHandle heh = voh.handle();
            int eidx = mesh_->edge_handle(heh).idx();
            spokeidx.push_back(eidx);

            OMMesh::VertexOHalfedgeIter nextoh = voh;
            ++nextoh;
            if(!nextoh)
                nextoh = mesh_->voh_iter(vi.handle());

            OMMesh::VertexHandle nextvert = mesh_->to_vertex_handle(nextoh.handle());

            OMMesh::HalfedgeHandle opp = mesh_->next_halfedge_handle(heh);;
            if(mesh_->to_vertex_handle(opp) != nextvert)
            {
                opp = mesh_->prev_halfedge_handle(mesh_->opposite_halfedge_handle(heh));
                assert(mesh_->from_vertex_handle(opp) == nextvert);
            }

            int oidx = mesh_->edge_handle(opp).idx();
            rightoppidx.push_back(oidx);

            OMMesh::VertexHandle vh = mesh_->to_vertex_handle(heh);
            nbidx.push_back(vh.idx());
        }

        int centidx = vi.handle().idx();

        energyB += ElasticEnergy::bendingEnergy(q, g, centidx, nbidx, spokeidx, rightoppidx, gradq, Hqcoeffs, dgdqcoeffs, params_, derivativesRequested);
    }

    // Stretching energy
    for(OMMesh::FaceIter it = mesh_->faces_begin(); it != mesh_->faces_end(); ++it)
    {
        int qidx[3];
        int gidx[3];

        int idx=0;
        for(OMMesh::FaceHalfedgeIter fhi = mesh_->fh_iter(it.handle()); fhi; ++fhi)
        {
            assert(idx < 3);
            OMMesh::HalfedgeHandle heh = fhi.handle();
            OMMesh::EdgeHandle eh = mesh_->edge_handle(heh);
            OMMesh::VertexHandle from = mesh_->from_vertex_handle(heh);
            gidx[idx] = eh.idx();
            qidx[(idx+1)%3] = from.idx();
            idx++;
        }
        assert(idx == 3);

        energyS += ElasticEnergy::stretchingEnergy(q, g, qidx, gidx, gradq, Hqcoeffs, dgdqcoeffs, params_, derivativesRequested);
    }

    if(derivativesRequested & ElasticEnergy::DR_HQ)
        hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    if(derivativesRequested & ElasticEnergy::DR_DGDQ)
        gradggradq.setFromTriplets(dgdqcoeffs.begin(), dgdqcoeffs.end());
}


bool Mesh::findMode(void)
{
    // Mq'' = HV
    // Bq = 0
    setIntrinsicLengthsToCurrentLengths();
    VectorXd q(3*numverts());
    VectorXd g(numedges());
    dofsFromGeometry(q,g);
    SparseMatrix<double> M(3*numverts(),3*numverts());
    buildExtendedMassMatrix(q, M);

    double energyB, energyS;

    int derivs = ElasticEnergy::DR_HQ;

    std::cout << "Calculating energy hessian" << std::endl;

    SparseMatrix<double> hessq, gradggradq;
    VectorXd gradq;

    elasticEnergy(q, g, energyB, energyS, gradq, hessq, gradggradq, derivs);

    std::cout << "Building matrices" << std::endl;

    MatrixXd left = hessq;
    MatrixXd right = params_.rho*params_.h*M;

    double leftscale = left.trace()/numdofs();
    double rightscale = right.trace()/numdofs();

    std::cout << "Scales: " << leftscale << " " << rightscale << std::endl;

    left /= leftscale;
    right /= rightscale;

    std::cout << "Solving for spectrum" << std::endl;
    GeneralizedSelfAdjointEigenSolver<MatrixXd > solver(left, right);

    MatrixXd rawmodes = solver.eigenvectors();
    modeFrequencies_ = solver.eigenvalues()*leftscale/rightscale;

    std::cout << "Done" << std::endl;
    for(int i=0; i<modeFrequencies_.size(); i++)
    {
        if(modeFrequencies_[i] < 0)
            modeFrequencies_[i] = 0;
        else
            modeFrequencies_[i] = sqrt(modeFrequencies_[i]);
    }

    modes_ = rawmodes;

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
        double area = circumcentricDualArea(q, vidx);
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
        double area = circumcentricDualArea(q, vidx);
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

double Mesh::circumcentricDualArea(const VectorXd &q, int vidx) const
{
    double result = 0;
    OMMesh::VertexHandle vh = mesh_->vertex_handle(vidx);
    for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vh); voh; ++voh)
    {
        if(mesh_->is_boundary(voh.handle()))
            continue;

        OMMesh::VertexHandle v1 = mesh_->to_vertex_handle(voh.handle());
        OMMesh::HalfedgeHandle nextheh = mesh_->next_halfedge_handle(voh.handle());
        OMMesh::VertexHandle v2 = mesh_->to_vertex_handle(nextheh);

        Vector3d cp = q.segment<3>(3*vh.idx());
        Vector3d v1p = q.segment<3>(3*v1.idx());
        Vector3d v2p = q.segment<3>(3*v2.idx());

        double sinalpha = ((cp-v1p).cross(v2p-v1p)).norm();
        double cosalpha = (cp-v1p).dot(v2p-v1p);
        double cotalpha = cosalpha/sinalpha;

        double sinbeta = ((cp-v2p).cross(v1p-v2p)).norm();
        double cosbeta = (cp-v2p).dot(v1p-v2p);
        double cotbeta = cosbeta/sinbeta;
        double v = (cp-v1p).norm();
        double u = (cp-v2p).norm();
        result += (u*u*cotalpha+v*v*cotbeta)/8.0;
    }
    return result;
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
