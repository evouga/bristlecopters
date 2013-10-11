#include "mesh.h"
#include <fadbad.h>
#include <fadiff.h>
#include <badiff.h>
#include <iomanip>
#include "autodifftemplates.h"
#include "controller.h"
#include <Eigen/Dense>
#include "newton.h"

using namespace std;
using namespace Eigen;
using namespace fadbad;
using namespace OpenMesh;

double Mesh::stretchOne(const VectorXd &qs, const VectorXd &gs, int *qidx, int *gidx, VectorXd &dq, std::vector<Tr> &hq, std::vector<Tr> &dgdq, bool derivs) const
{
    // sqrt(det g) tr(g^-1 a - I)^2

    double g[3];
    Vector3d q[3];
    for(int i=0; i<3; i++)
    {
        g[i] = gs[gidx[i]];
        q[i] = qs.segment<3>(3*qidx[i]);
    }

    double detg = g[0]*g[0]*g[1]*g[1] - (g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])/4.0;

    double A = (q[2]-q[1]).dot(2.0*g[1]*g[1]*(q[2]-q[1])+(g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[0]))
            +(q[2]-q[0]).dot(2.0*g[0]*g[0]*(q[2]-q[0])+(g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[1]))
            -4.0*detg;

    if(derivs)
    {
        Vector3d dAdq[3];
        dAdq[0] = 2.0*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(q[2]-q[1])-4.0*g[0]*g[0]*(q[2]-q[0]);
        dAdq[1] = 2.0*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(q[2]-q[0])-4.0*g[1]*g[1]*(q[2]-q[1]);
        dAdq[2] = 2.0*(g[1]*g[1]+g[2]*g[2]-g[0]*g[0])*(q[2]-q[1])+2.0*(g[0]*g[0]+g[2]*g[2]-g[1]*g[1])*(q[2]-q[0]);

        double dqcoeff = A/(4.0*detg*sqrt(detg));

        for(int i=0; i<3; i++)
            dq.segment<3>(3*qidx[i]) += dqcoeff*dAdq[i];

        for(int i=0; i<3; i++)
        {
            Matrix3d dqidqi;
            dqidqi.setIdentity();
            dqidqi *= A*g[i]*g[i];
            dqidqi += 0.25*(dAdq[i]*dAdq[i].transpose());
            dqidqi *= 1.0/(detg*sqrt(detg));

            for(int j=0; j<3; j++)
                for(int k=0; k<3; k++)
                {
                    hq.push_back(Tr(3*qidx[i]+j, 3*qidx[i]+k, dqidqi(k,j)));
                }
        }

        Matrix3d dq0dq1;
        dq0dq1.setIdentity();
        dq0dq1 *= 0.5*(g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*A;
        dq0dq1 += 0.25*dAdq[1]*dAdq[0].transpose();
        dq0dq1 *= 1.0/(detg*sqrt(detg));
        for(int j=0; j<3; j++)
        {
            for(int k=0; k<3; k++)
            {
                hq.push_back(Tr(3*qidx[0]+j, 3*qidx[1]+k, dq0dq1(k,j)));
                hq.push_back(Tr(3*qidx[1]+k, 3*qidx[0]+j, dq0dq1(k,j)));
            }
        }

        for(int i=0; i<2; i++)
        {
            Matrix3d dq2dqi;
            dq2dqi.setIdentity();
            dq2dqi *= 0.5*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2]-2.0*g[i]*g[i])*A;
            dq2dqi += 0.25*dAdq[i]*dAdq[2].transpose();
            dq2dqi *= 1.0/(detg*sqrt(detg));
            for(int j=0; j<3; j++)
                for(int k=0; k<3; k++)
                {
                    hq.push_back(Tr(3*qidx[2]+j, 3*qidx[i]+k, dq2dqi(k,j)));
                    hq.push_back(Tr(3*qidx[i]+k, 3*qidx[2]+j, dq2dqi(k,j)));
                }
        }

        Vector3d dgdetg(g[0]*(g[1]*g[1]+g[2]*g[2]-g[0]*g[0]),
                g[1]*(g[0]*g[0]+g[2]*g[2]-g[1]*g[1]),
                g[2]*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2]));

        Vector3d dgA = -4.0*dgdetg;
        dgA[0] += 4.0*(g[0]*(q[0]-q[2]).dot(q[0]-q[1]));
        dgA[1] += 4.0*(g[1]*(q[1]-q[2]).dot(q[1]-q[0]));
        dgA[2] += 4.0*(g[2]*(q[2]-q[1]).dot(q[2]-q[0]));

        Matrix3d dgdqA[3];
        dgdqA[0].col(0) = 4.0*g[0]*(q[0]-q[2])+4.0*g[0]*(q[0]-q[1]);
        dgdqA[0].col(1) = 4.0*g[1]*(q[2]-q[1]);
        dgdqA[0].col(2) = 4.0*g[2]*(q[1]-q[2]);

        dgdqA[1].col(0) = 4.0*g[0]*(q[2]-q[0]);
        dgdqA[1].col(1) = 4.0*g[1]*(q[1]-q[0])+4.0*g[1]*(q[1]-q[2]);
        dgdqA[1].col(2) = 4.0*g[2]*(q[0]-q[2]);

        dgdqA[2].col(0) = 4.0*g[0]*(q[1]-q[0]);
        dgdqA[2].col(1) = 4.0*g[1]*(q[0]-q[1]);
        dgdqA[2].col(2) = 4.0*g[2]*(q[2]-q[0]) + 4.0*g[2]*(q[2]-q[1]);

        for(int i=0; i<3; i++)
        {
            Matrix3d dgdqi = A/(4.0*detg*sqrt(detg)) * dgdqA[i];
            dgdqi += 1.0/(4.0*detg*sqrt(detg)) * dAdq[i]*dgA.transpose();
            dgdqi -= 3.0*A/(8.0 * detg * detg * sqrt(detg)) * dAdq[i]*dgdetg.transpose();

            for(int j=0; j<3; j++)
                for(int k=0; k<3; k++)
                {
                    dgdq.push_back(Tr(gidx[j], 3*qidx[i]+k, dgdqi(k,j)));
                }
        }
    }

    return 0.5*A*A/(4.0*detg*sqrt(detg));
}

double Mesh::stretchTwo(const VectorXd &qs, const VectorXd &gs, int *qidx, int *gidx, VectorXd &dq, std::vector<Tr> &hq, std::vector<Tr> &dgdq, bool derivs) const
{
    // 0.5 sqrt(det g) det(g^-1 a - I)

    double g[3];
    Vector3d q[3];
    for(int i=0; i<3; i++)
    {
        g[i] = gs[gidx[i]];
        q[i] = qs.segment<3>(3*qidx[i]);
    }

    double detg = g[0]*g[0]*g[1]*g[1] - (g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2])/4.0;

    double A = (q[2]-q[1]).dot(q[2]-q[1])*(q[2]-q[0]).dot(q[2]-q[0])-(q[2]-q[1]).dot(q[2]-q[0])*(q[2]-q[1]).dot(q[2]-q[0]);
    double B = (g[0]*g[0]+g[1]*g[1]-g[2]*g[2])*(q[2]-q[1]).dot(q[2]-q[0]) - g[0]*g[0]*(q[2]-q[0]).dot(q[2]-q[0]) - g[1]*g[1]*(q[2]-q[1]).dot(q[2]-q[1]);

    if(derivs)
    {
        Vector3d dqA[3];
        dqA[0] = 2.0*(q[2]-q[1]).dot(q[2]-q[0])*(q[2]-q[1]) - 2.0*(q[2]-q[1]).dot(q[2]-q[1])*(q[2]-q[0]);
        dqA[1] = 2.0*(q[2]-q[1]).dot(q[2]-q[0])*(q[2]-q[0]) - 2.0*(q[2]-q[0]).dot(q[2]-q[0])*(q[2]-q[1]);
        dqA[2] = 2.0*(q[2]-q[0]).dot(q[1]-q[0])*(q[2]-q[1]) - 2.0*(q[2]-q[1]).dot(q[1]-q[0])*(q[2]-q[0]);

        Vector3d dqB[3];
        dqB[0] = (g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[1]) + 2.0*g[0]*g[0]*(q[2]-q[0]);
        dqB[1] = (g[2]*g[2]-g[0]*g[0]-g[1]*g[1])*(q[2]-q[0]) + 2.0*g[1]*g[1]*(q[2]-q[1]);
        dqB[2] = (g[1]*g[1]-g[2]*g[2]-g[0]*g[0])*(q[2]-q[0]) + (g[0]*g[0]-g[1]*g[1]-g[2]*g[2])*(q[2]-q[1]);

        for(int i=0; i<3; i++)
            dq.segment<3>(3*qidx[i]) += (dqA[i]+dqB[i])/(2.0*sqrt(detg));

        Matrix3d dqdq0A[3];
        dqdq0A[0].setIdentity();
        dqdq0A[0] *= 2.0*(q[2]-q[1]).dot(q[2]-q[1]);
        dqdq0A[0] -= 2.0*(q[2]-q[1])*(q[2]-q[1]).transpose();

        dqdq0A[1].setIdentity();
        dqdq0A[1] *= -2.0*(q[2]-q[1]).dot(q[2]-q[0]);
        dqdq0A[1] += 4.0*(q[2]-q[0])*(q[2]-q[1]).transpose()-2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();

        dqdq0A[2].setIdentity();
        dqdq0A[2] *= 2.0*(q[2]-q[1]).dot(q[1]-q[0]);
        dqdq0A[2] += 2.0*(q[2]-q[1])*(q[2]-q[1]).transpose();
        dqdq0A[2] += 2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
        dqdq0A[2] -= 4.0*(q[2]-q[0])*(q[2]-q[1]).transpose();

        Matrix3d dqdq1A[3];
        dqdq1A[0].setIdentity();
        dqdq1A[0] *= -2.0*(q[2]-q[1]).dot(q[2]-q[0]);
        dqdq1A[0] += 4.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
        dqdq1A[0] += -2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();

        dqdq1A[1].setIdentity();
        dqdq1A[1] *= 2.0*(q[2]-q[0]).dot(q[2]-q[0]);
        dqdq1A[1] += -2.0*(q[2]-q[0])*(q[2]-q[0]).transpose();

        dqdq1A[2].setIdentity();
        dqdq1A[2] *= 2.0*(q[2]-q[0]).dot(q[0]-q[1]);
        dqdq1A[2] += 2.0*(q[2]-q[0])*(q[2]-q[0]).transpose();
        dqdq1A[2] += 2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();
        dqdq1A[2] += -4.0*(q[2]-q[1])*(q[2]-q[0]).transpose();

        Matrix3d dqdq2A[3];
        dqdq2A[0].setIdentity();
        dqdq2A[0] *= 2.0*(q[2]-q[1]).dot(q[1]-q[0]);
        dqdq2A[0] += -2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
        dqdq2A[0] += -2.0*(q[2]-q[1])*(q[1]-q[0]).transpose();
        dqdq2A[0] += 2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();

        dqdq2A[1].setIdentity();
        dqdq2A[1] *= -2.0*(q[2]-q[0]).dot(q[1]-q[0]);
        dqdq2A[1] += 2.0*(q[2]-q[1])*(q[2]-q[0]).transpose();
        dqdq2A[1] += -2.0*(q[2]-q[0])*(q[2]-q[1]).transpose();
        dqdq2A[1] += 2.0*(q[2]-q[0])*(q[1]-q[0]).transpose();

        dqdq2A[2].setIdentity();
        dqdq2A[2] *= 2.0*(q[1]-q[0]).dot(q[1]-q[0]);
        dqdq2A[2] += -2.0*(q[1]-q[0])*(q[1]-q[0]).transpose();

        Matrix3d dqdq0B[3];
        dqdq0B[0].setIdentity();
        dqdq0B[0] *= -2.0*g[0]*g[0];

        dqdq0B[1].setIdentity();
        dqdq0B[1] *= (g[0]*g[0]+g[1]*g[1]-g[2]*g[2]);

        dqdq0B[2].setIdentity();
        dqdq0B[2] *= (g[0]*g[0]+g[2]*g[2]-g[1]*g[1]);

        Matrix3d dqdq1B[3];
        dqdq1B[0].setIdentity();
        dqdq1B[0] *= (g[0]*g[0]+g[1]*g[1]-g[2]*g[2]);

        dqdq1B[1].setIdentity();
        dqdq1B[1] *= -2.0*g[1]*g[1];

        dqdq1B[2].setIdentity();
        dqdq1B[2] *= (g[1]*g[1]+g[2]*g[2]-g[0]*g[0]);

        Matrix3d dqdq2B[3];
        dqdq2B[0].setIdentity();
        dqdq2B[0] *= (g[0]*g[0]+g[2]*g[2]-g[1]*g[1]);

        dqdq2B[1].setIdentity();
        dqdq2B[1] *= (g[1]*g[1]+g[2]*g[2]-g[0]*g[0]);

        dqdq2B[2].setIdentity();
        dqdq2B[2] *= -2.0*g[2]*g[2];

        for(int i=0; i<3; i++)
        {
            Matrix3d dqdq0E = (dqdq0A[i]+dqdq0B[i])/(2.0*sqrt(detg));
            Matrix3d dqdq1E = (dqdq1A[i]+dqdq1B[i])/(2.0*sqrt(detg));
            Matrix3d dqdq2E = (dqdq2A[i]+dqdq2B[i])/(2.0*sqrt(detg));
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    hq.push_back(Tr(3*qidx[i]+j, 3*qidx[0]+k, dqdq0E(k,j)));
                    hq.push_back(Tr(3*qidx[i]+j, 3*qidx[1]+k, dqdq1E(k,j)));
                    hq.push_back(Tr(3*qidx[i]+j, 3*qidx[2]+k, dqdq2E(k,j)));
                }
            }
        }

        Matrix3d dgdqB[3];

        dgdqB[0].col(0) = 2.0*g[0]*(q[2]-q[0]) + 2.0*g[0]*(q[1]-q[0]);
        dgdqB[0].col(1) = -2.0*g[1]*(q[2]-q[1]);
        dgdqB[0].col(2) = 2.0*g[2]*(q[2]-q[1]);

        dgdqB[1].col(0) = -2.0*g[0]*(q[2]-q[0]);
        dgdqB[1].col(1) = 2.0*g[1]*(q[2]-q[1]) + 2.0*g[1]*(q[0]-q[1]);
        dgdqB[1].col(2) = 2.0*g[2]*(q[2]-q[0]);

        dgdqB[2].col(0) = -2.0*g[0]*(q[1]-q[0]);
        dgdqB[2].col(1) = 2.0*g[1]*(q[1]-q[0]);
        dgdqB[2].col(2) = 2.0*g[2]*(q[0]-q[2]) + 2.0*g[2]*(q[1]-q[2]);

        Vector3d dgdetg(g[0]*(g[1]*g[1]+g[2]*g[2]-g[0]*g[0]),
                g[1]*(g[0]*g[0]+g[2]*g[2]-g[1]*g[1]),
                g[2]*(g[0]*g[0]+g[1]*g[1]-g[2]*g[2]));


        for(int i=0; i<3; i++)
        {
            Matrix3d dgdqi = -1.0/(4.0*detg*sqrt(detg))*(dqA[i]+dqB[i])*dgdetg.transpose();
            dgdqi += 1.0/(2.0*sqrt(detg)) * dgdqB[i];
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    dgdq.push_back(Tr(gidx[j], 3*qidx[i]+k, dgdqi(k,j)));
                }
            }
        }
    }

    return 0.5*sqrt(detg) + 0.5*(A+B)/sqrt(detg);
}


void Mesh::elasticEnergy(const VectorXd &q, const VectorXd &g, double &energyB, double &energyS) const
{
    EnergyDerivatives derivs = NONE;
    VectorXd gradq, gradg;
    SparseMatrix<double> hessq, hessg, gradggradq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergyG(const VectorXd &q,
                          const VectorXd &g,
                          double &energyB,
                          double &energyS,
                          VectorXd &gradg,
                          Eigen::SparseMatrix<double> &hessg) const
{
    EnergyDerivatives derivs = G;
    VectorXd gradq;
    SparseMatrix<double> hessq, gradggradq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergyQ(const VectorXd &q,
                          const VectorXd &g,
                          double &energyB,
                          double &energyS,
                          VectorXd &gradq,
                          Eigen::SparseMatrix<double> &hessq) const
{
    EnergyDerivatives derivs = Q;
    VectorXd gradg;
    SparseMatrix<double> hessg, gradggradq;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergyGQ(const VectorXd &q, const VectorXd &g, VectorXd &gradq, Eigen::SparseMatrix<double> &gradggradq)
{
    double energyB, energyS;
    VectorXd gradg;
    SparseMatrix<double> hessg, hessq;
    EnergyDerivatives derivs = BOTH;
    elasticEnergy(q, g, energyB, energyS, gradq, gradg, hessq, hessg, gradggradq, derivs);
}

void Mesh::elasticEnergy(const VectorXd &q,
                         const VectorXd &g,
                         double &energyB,
                         double &energyS,
                         VectorXd &gradq,
                         VectorXd &gradg,
                         Eigen::SparseMatrix<double> &hessq,
                         Eigen::SparseMatrix<double> &hessg,
                         Eigen::SparseMatrix<double> &gradggradq,
                         EnergyDerivatives derivs) const
{
    assert(q.size() == numdofs());
    assert(g.size() == numedges());
    energyB = energyS = 0;

    if(derivs & Q)
    {
        gradq.resize(numdofs());
        gradq.setZero();
        hessq.resize(numdofs(), numdofs());
    }

    if(derivs & G)
    {
        gradg.resize(numedges());
        gradg.setZero();
        hessg.resize(numedges(), numedges());
    }

    if(derivs & G && derivs & Q)
    {
        gradggradq.resize(numedges(), numdofs());
    }


    vector<Tr> Hqcoeffs;
    vector<Tr> Hgcoeffs;
    vector<Tr> dgdqcoeffs;

    // bending energy
    double bendcoeff = params_.h*params_.h*params_.YoungsModulus/24.0/(1.0+params_.PoissonRatio);
    for(OMMesh::VertexIter vi = mesh_->vertices_begin(); vi != mesh_->vertices_end(); ++vi)
    {
        if(mesh_->is_boundary(vi.handle()))
            continue;

        B<F<double> > stencilenergy = 0;
        // bending energy Laplacian term
        double Lcoeff = 1.0/(1.0-params_.PoissonRatio);
        vector<double> spokelens;
        vector<int> spokeidx;
        vector<double> rightopplens;
        vector<int> rightoppidx;

        vector<Vector3d> nbq;
        vector<int> nbidx;
        for(OMMesh::VertexOHalfedgeIter voh = mesh_->voh_iter(vi.handle()); voh; ++voh)
        {
            OMMesh::HalfedgeHandle heh = voh.handle();
            int eidx = mesh_->edge_handle(heh).idx();
            spokelens.push_back(g[eidx]);
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
            rightopplens.push_back(g[oidx]);
            rightoppidx.push_back(oidx);

            OMMesh::VertexHandle vh = mesh_->to_vertex_handle(heh);
            nbq.push_back(q.segment<3>(3*vh.idx()));
            nbidx.push_back(vh.idx());
        }

        int numnbs = (int)spokelens.size();
        int diffqvars;
        if(derivs & Q)
            diffqvars = 3*(numnbs+1);
        else
            diffqvars = 0;

        int diffvars;
        if(derivs & G)
            diffvars = diffqvars + 2*numnbs;
        else
            diffvars = diffqvars;

        F<double> *dspokelens = new F<double>[numnbs];
        B<F<double> > *ddspokelens = new B<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dspokelens[i] = spokelens[i];
            if(derivs & G)
                dspokelens[i].diff(diffqvars+i, diffvars);
            ddspokelens[i] = dspokelens[i];
            //if(derivs & G)
            //    ddspokelens[i].diff(diffqvars+i, diffvars);
        }
        F<double> *dopplens = new F<double>[numnbs];
        B<F<double> > *ddopplens = new B<F<double> >[numnbs];
        for(int i=0; i<numnbs; i++)
        {
            dopplens[i] = rightopplens[i];
            if(derivs & G)
                dopplens[i].diff(diffqvars+numnbs+i, diffvars);
            ddopplens[i] = dopplens[i];
            //if(derivs & G)
            //    ddopplens[i].diff(diffqvars+numnbs+i, diffvars);
        }

        F<double> dq[3];
        B<F<double> > ddq[3];
        int vidx = vi.handle().idx();
        Vector3d thisq = q.segment<3>(3*vi.handle().idx());
        for(int i=0; i<3; i++)
        {
            dq[i] = thisq[i];
            if(derivs & Q)
                dq[i].diff(i,diffvars);
            ddq[i] = dq[i];
        }

        F<double> *dnbq[3];
        B<F<double> > *ddnbq[3];
        for(int i=0; i<3; i++)
        {
            dnbq[i] = new F<double>[numnbs];
            ddnbq[i] = new B<F<double> >[numnbs];
            for(int j=0; j<numnbs; j++)
            {
                dnbq[i][j] = nbq[j][i];
                if(derivs & Q)
                    dnbq[i][j].diff(3+3*j+i,diffvars);
                ddnbq[i][j] = dnbq[i][j];
            }
        }

        {
            B<F<double> > Lr[3];
            for(int j=0; j<3; j++)
                Lr[j] = L(ddq[j], numnbs, ddnbq[j], ddspokelens, ddopplens);
            B<F<double> > materialArea = dualbarycentricarea(numnbs, ddspokelens, ddopplens);

            stencilenergy += bendcoeff*Lcoeff*normSquared(Lr)/materialArea;

            // bending energy det term
            vector<B<F<double> > *> nbddq;
            for(int i=0; i<numnbs; i++)
            {
                B<F<double> > *ddq = new B<F<double> >[3];
                nbddq.push_back(ddq);
                for(int j=0; j<3; j++)
                {
                    ddq[j] = ddnbq[j][i];
                }
            }

            B<F<double> > Kcurv = K(ddq, nbddq);
            B<F<double> > embeddedArea = dualbarycentricarea(ddq, nbddq);
            stencilenergy += -2.0*bendcoeff*Kcurv*embeddedArea*embeddedArea/materialArea;
            for(int i=0; i<numnbs; i++)
                delete[] nbddq[i];
        }

        stencilenergy.diff(0,1);

        energyB += stencilenergy.val().val();

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                gradq[3*vidx+j] += ddq[j].d(0).val();
            }
        }

        for(int i=0; i<numnbs; i++)
        {
            if(derivs & Q)
            {
                for(int j=0; j<3; j++)
                {
                    gradq[3*nbidx[i]+j] += ddnbq[j][i].d(0).val();
                }
            }
            if(derivs & G)
            {
                gradg[spokeidx[i]] += ddspokelens[i].d(0).val();
                gradg[rightoppidx[i]] += ddopplens[i].d(0).val();
            }
        }

        if(derivs & Q)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    double hess = ddq[j].d(0).d(k);
                    if(hess != 0)
                        Hqcoeffs.push_back(Tr(3*vidx+j,3*vidx+k,hess));
                }
                for(int k=0; k<numnbs; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        double hess = ddq[j].d(0).d(3+3*k+l);
                        if(hess != 0)
                        {
                            Hqcoeffs.push_back(Tr(3*vidx+j,3*nbidx[k]+l,hess));
                            Hqcoeffs.push_back(Tr(3*nbidx[k]+l,3*vidx+j,hess));
                        }
                    }
                }
            }

            for(int i=0; i<numnbs; i++)
            {
                for(int j=0; j<numnbs; j++)
                {
                    for(int k=0; k<3; k++)
                    {
                        for(int l=0; l<3; l++)
                        {
                            double hess = ddnbq[k][i].d(0).d(3+3*j+l);
                            if(hess != 0)
                                Hqcoeffs.push_back(Tr(3*nbidx[i]+k,3*nbidx[j]+l,hess));
                        }
                    }
                }
            }
        }

        if(derivs & G)
        {
            for(int i=0; i<numnbs; i++)
            {
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddspokelens[i].d(0).d(diffqvars+j);
                    if(hess != 0)
                        Hgcoeffs.push_back(Tr(spokeidx[i],spokeidx[j],hess));
                }
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddspokelens[i].d(0).d(diffqvars+numnbs+j);
                    if(hess != 0)
                    {
                        Hgcoeffs.push_back(Tr(spokeidx[i],rightoppidx[j],hess));
                        Hgcoeffs.push_back(Tr(rightoppidx[j],spokeidx[i],hess));
                    }
                }
                for(int j=0; j<numnbs; j++)
                {
                    double hess = ddopplens[i].d(0).d(diffqvars+numnbs+j);
                    if(hess != 0)
                        Hgcoeffs.push_back(Tr(rightoppidx[i],rightoppidx[j],hess));
                }

                if(derivs & Q)
                {
                    for(int j=0; j<3; j++)
                    {
                        double hess = ddspokelens[i].d(0).d(j);
                        if(hess != 0)
                            dgdqcoeffs.push_back(Tr(spokeidx[i],3*vidx+j,hess));
                        hess = ddopplens[i].d(0).d(j);
                        if(hess != 0)
                            dgdqcoeffs.push_back(Tr(rightoppidx[i],3*vidx+j,hess));
                    }

                    for(int j=0; j<numnbs; j++)
                    {
                        for(int k=0; k<3; k++)
                        {
                            double hess = ddspokelens[i].d(0).d(3+3*j+k);
                            if(hess != 0)
                                dgdqcoeffs.push_back(Tr(spokeidx[i], 3*nbidx[j]+k, hess));

                            hess = ddopplens[i].d(0).d(3+3*j+k);
                            if(hess != 0)
                                dgdqcoeffs.push_back(Tr(rightoppidx[i], 3*nbidx[j]+k, hess));
                        }
                    }
                }
            }
        }

        delete[] dspokelens;
        delete[] ddspokelens;
        delete[] dopplens;
        delete[] ddopplens;
        for(int i=0; i<3; i++)
        {
            delete[] dnbq[i];
            delete[] ddnbq[i];
        }

    }

    // Stretching energy
    double stretchcoeff = params_.YoungsModulus/8.0/(1.0+params_.PoissonRatio);

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        B<F<double> > stencilenergy = 0;
        int vidx[3];
        Vector3d pts[3];
        double elens[3];
        int eidx[3];
        int idx=0;
        for(OMMesh::FaceHalfedgeIter fei = mesh_->fh_iter(fi.handle()); fei; ++fei)
        {
            vidx[idx] = mesh_->to_vertex_handle(fei.handle()).idx();
            pts[idx] = q.segment<3>(3*mesh_->to_vertex_handle(fei.handle()).idx());
            int eid = mesh_->edge_handle(fei.handle()).idx();
            eidx[idx] = eid;
            elens[idx] = g[eid];
            idx++;
        }
        assert(idx==3);

        int diffqvars;
        if(derivs & Q)
            diffqvars = 9;
        else
            diffqvars = 0;

        int diffvars;
        if(derivs & G)
            diffvars = diffqvars + 3;
        else
            diffvars = diffqvars;

        F<double> dq[3][3];
        B<F<double> > ddq[3][3];
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                dq[i][j] = pts[i][j];
                if(derivs & Q)
                    dq[i][j].diff(3*i+j, diffvars);
                ddq[i][j] = dq[i][j];
            }
        }

        F<double> delen[3];
        B<F<double> > ddelen[3];
        for(int i=0; i<3; i++)
        {
            delen[i] = elens[i];
            if(derivs & G)
                delen[i].diff(diffqvars+i, diffvars);
            ddelen[i] = delen[i];
        }

        {
            B<F<double> > emblens[3];
            for(int i=0; i<3; i++)
            {
                int previdx = (3+i-1)%3;
                B<F<double> > vec[3];
                diff(ddq[i],ddq[previdx],vec);
                emblens[i] = norm(vec);
            }

            B<F<double> > g[4];
            B<F<double> > a[4];
            metrictri(ddelen[0], ddelen[1], ddelen[2], g);
            metrictri(emblens[0], emblens[1], emblens[2], a);
            B<F<double> > ginva[4];
            invGtimesH(g, a, ginva);
            ginva[0] -= 1.0;
            ginva[3] -= 1.0;
            B<F<double> > matarea = 0.5*sqrt(det(g));
            //stencilenergy += matarea*stretchcoeff/(1.0-params_.PoissonRatio)*tr(ginva)*tr(ginva);
            stencilenergy += matarea*stretchcoeff*-2.0*det(ginva);
        }

        stencilenergy.diff(0,1);

        energyS += stencilenergy.val().val();

        VectorXd sdq(numdofs());
        sdq.setZero();
        vector<Tr> hq;
        vector<Tr> dgdq;
        int qidx[3];
        qidx[0] = vidx[1];
        qidx[1] = vidx[2];
        qidx[2] = vidx[0];
        double checkenergy = stretchTwo(q, g, qidx, eidx, sdq, hq, dgdq, true);
        SparseMatrix<double> Hq(numdofs(), numdofs());
        Hq.setFromTriplets(hq.begin(), hq.end());
        Hq  *= -2.0*stretchcoeff;
        SparseMatrix<double> DgDq(numedges(), numdofs());
        DgDq.setFromTriplets(dgdq.begin(), dgdq.end());
        DgDq *= -2.0*stretchcoeff;

        std::cout << stretchcoeff*-2.0*checkenergy << " " << stencilenergy.val().val() << endl;

        for(int i=0; i<3; i++)
        {
            if(derivs & Q)
            {
                std::cout << -2.0*stretchcoeff*sdq.segment<3>(3*vidx[i]).transpose() << " ";
                for(int j=0; j<3; j++)
                {
                    gradq[3*vidx[i]+j] += ddq[i][j].d(0).val();
                    std::cout << ddq[i][j].d(0).val() << " ";
                }
                std::cout << std::endl;
            }
            if(derivs & G)
                gradg[eidx[i]] += ddelen[i].d(0).val();
        }

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                if(derivs & Q)
                {
                    for(int k=0; k<3; k++)
                    {
                        for(int l=0; l<3; l++)
                        {
                            double hess = ddq[i][k].d(0).d(3*j+l);
                            if(hess != 0.0)
                            {
                                Hqcoeffs.push_back(Tr(3*vidx[i]+k,3*vidx[j]+l,hess));
                                std::cout << "H" << hess << " " << Hq.coeffRef(3*vidx[i]+k, 3*vidx[j]+l) << endl;
                            }
                        }
                    }
                }

                if(derivs & G)
                {
                    double hess = ddelen[i].d(0).d(diffqvars+j);
                    if(hess != 0)
                    {
                        Hgcoeffs.push_back(Tr(eidx[i],eidx[j],hess));
                    }                                       
                }
            }

            if(derivs & G && derivs & Q)
            {
                for(int k=0; k<3; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        double hess = ddelen[i].d(0).d(3*k+l);
                        if(hess != 0.0)
                        {
                            dgdqcoeffs.push_back(Tr(eidx[i], 3*vidx[k]+l, hess));
                            std::cout << "dd " << hess << " " << DgDq.coeffRef(eidx[i], 3*vidx[k]+l) << endl;
                        }
                    }
                }
            }
        }
    }

    if(derivs & Q)
        hessq.setFromTriplets(Hqcoeffs.begin(), Hqcoeffs.end());
    if(derivs & G)
        hessg.setFromTriplets(Hgcoeffs.begin(), Hgcoeffs.end());
    if(derivs & G && derivs & Q)
        gradggradq.setFromTriplets(dgdqcoeffs.begin(), dgdqcoeffs.end());
}

double Mesh::triangleInequalityLineSearch(const VectorXd &g, const VectorXd &dg) const
{
    assert(g.size() == numedges());
    assert(dg.size() == numedges());

    double maxt = std::numeric_limits<double>::infinity();

    for(OMMesh::FaceIter fi = mesh_->faces_begin(); fi != mesh_->faces_end(); ++fi)
    {
        vector<double> gs;
        vector<double> dgs;
        for(OMMesh::FaceEdgeIter fei = mesh_->fe_iter(fi.handle()); fei; ++fei)
        {
            gs.push_back(g[fei.handle().idx()]);
            dgs.push_back(dg[fei.handle().idx()]);
        }
        assert(gs.size() == 3);
        for(int i=0; i<3; i++)
        {
            int idx[3];
            for(int j=0; j<3; j++)
            {
                idx[j] = (i+j)%3;
            }
            double thismax = triangleInequalityLineSearch(gs[idx[0]], gs[idx[1]], gs[idx[2]],
                    dgs[idx[0]], dgs[idx[1]], dgs[idx[2]]);

            maxt = std::min(maxt, thismax);
        }
    }

    return maxt;
}

double Mesh::triangleInequalityLineSearch(double g0, double g1, double g2, double dg0, double dg1, double dg2) const
{
    double num = g0+g1-g2;
    double denom = dg2-dg0-dg1;
    if(denom == 0)
        return std::numeric_limits<double>::infinity();
    double cand = num/denom;
    if(cand < 0)
        return std::numeric_limits<double>::infinity();
    return cand;
}

class EmbeddingMinimizer : public NewtonObjective
{
public:
    EmbeddingMinimizer(Mesh &m, Controller &cont, const VectorXd &g) : m_(m), cont_(cont), g_(g) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        double energyB, energyS;
        VectorXd gradq;
        SparseMatrix<double> hessq;
        m_.elasticEnergyQ(q, g_, energyB, energyS, gradq, hessq);
        return energyB+energyS;
    }

    virtual void getGradient(const VectorXd &q, VectorXd &grad) const
    {
        double energyB, energyS;
        SparseMatrix<double> hessq;
        m_.elasticEnergyQ(q, g_, energyB, energyS, grad, hessq);
    }

    virtual void getHessian(const VectorXd &q, Eigen::SparseMatrix<double> &hess) const
    {
        double energyB, energyS;
        VectorXd grad;
        m_.elasticEnergyQ(q, g_, energyB, energyS, grad, hess);
    }

    virtual void showCurrentIteration(const VectorXd &q) const
    {
        m_.dofsToGeometry(q, g_);
        cont_.updateGL();
    }

private:
    Mesh &m_;
    Controller &cont_;
    const VectorXd &g_;
};

class MetricFit : public NewtonObjective
{
public:
    MetricFit(Mesh &m, Controller &cont, const VectorXd &q) : m_(m), cont_(cont), q_(q) {}

    virtual double getEnergy(const VectorXd &q) const
    {
        VectorXd gradq;
        SparseMatrix<double> dgdq;
        m_.elasticEnergyGQ(q_, q, gradq, dgdq);
        return gradq.norm();
    }

    virtual void getGradient(const VectorXd &q, VectorXd &grad) const
    {
        SparseMatrix<double> dgdq;
        VectorXd gradq;
        m_.elasticEnergyGQ(q_, q, gradq, dgdq);
        grad = dgdq*gradq;
    }

    virtual void getHessian(const VectorXd &q, Eigen::SparseMatrix<double> &hess) const
    {
        SparseMatrix<double> dgdq;
        VectorXd gradq;
        m_.elasticEnergyGQ(q_, q, gradq, dgdq);
        hess = dgdq*dgdq.transpose();
    }

    virtual void showCurrentIteration(const VectorXd &q) const
    {
        m_.dofsToGeometry(q_, q);
        cont_.updateGL();
    }

private:
    Mesh &m_;
    Controller &cont_;
    const VectorXd &q_;
};

bool Mesh::relaxEnergy(Controller &cont, RelaxationType type)
{
    VectorXd q(numdofs());
    VectorXd g(numedges());
    dofsFromGeometry(q, g);

    NewtonObjective *obj = NULL;
    VectorXd guess;

    if(type == RelaxEmbedding)
    {
        obj = new EmbeddingMinimizer(*this, cont, g);
        guess = q;
    }
    else if(type == FitMetric)
    {
        obj = new MetricFit(*this, cont, q);
        guess = g;
    }
    else
        return false;

    Newton n(*obj);
    NewtonParameters params;
    params.tol = params_.tol;
    params.maxiters = params_.maxiters;
    params.lsmaxiters = params_.maxlinesearchiters;
    VectorXd result;
    Newton::SolverStatus ss = n.solve(params, guess, result);
    std::cout << n.solverStatusMessage(ss) << std::endl;

    if(ss != Newton::CONVERGED)
        return false;

    if(type == RelaxEmbedding)
    {
        q = result;
    }
    else if(type == FitMetric)
    {
        g = result;
    }

    dofsToGeometry(q, g);
    cont.updateGL();
    return true;
}

bool Mesh::largestMagnitudeEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue)
{
    int dim = M.cols();
    VectorXd v(dim);
    v.setRandom();
    v.normalize();
    double curestimate=0;
    int iter=0;

    for(iter=0; iter < params_.maxpoweriters; iter++)
    {
        v = M*v;        
        v.normalize();

        curestimate = v.dot(M*v);
        double err = (M*v-curestimate*v).norm();
        if(err < params_.powertol)
            break;
    }
    eigenvalue = curestimate;
    return iter < params_.maxpoweriters;
}

bool Mesh::smallestEigenvalue(const Eigen::SparseMatrix<double> &M, double &eigenvalue)
{
    double largesteval;
    if(!largestMagnitudeEigenvalue(M, largesteval))
        return false;

    if(largesteval < 0)
    {
        eigenvalue = largesteval;
        return true;
    }

    int dim = M.cols();
    SparseMatrix<double> shift(dim,dim);
    shift.setIdentity();
    shift *= -largesteval;
    SparseMatrix<double> newM = M + shift;
    double newlargest;
    if(!largestMagnitudeEigenvalue(newM, newlargest))
        return false;
    eigenvalue = newlargest+largesteval;
    return true;
}
