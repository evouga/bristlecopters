#include "newton.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

using namespace Eigen;

Newton::Newton(const NewtonObjective &no) : no_(no)
{
}

Newton::SolverStatus Newton::solve(NewtonParameters params, const VectorXd &guess, VectorXd &result)
{
    VectorXd q = guess;

    result.resize(q.size());
    result.setZero();

    int iter=0;

    for(; iter<params.maxiters; iter++)
    {
        SparseMatrix<double> hessian;
        VectorXd gradient;
        double energy = no_.getEnergyAndDerivatives(q, gradient, hessian);
        double energycheck = no_.getEnergy(q);
        if(fabs(energy-energycheck) > 1e-10)
            return BAD_INPUT;

        if(isnan(energy))
            return BAD_INPUT;

        std::cout << std::fixed << std::setprecision(8) << "Iter " << iter+1
                  << "   E " << energy
                  << "   |dq| " << gradient.norm();


        if(gradient.norm() < params.tol)
        {
            std::cout << std::endl;
            break;
        }

        SimplicialLDLT<SparseMatrix<double> > solver;
        solver.compute(hessian);
        VectorXd searchdir = -solver.solve(gradient);

        double residual = (hessian*searchdir+gradient).norm();

        std::cout << "   r = " << residual << "   sd = " << searchdir.dot(gradient);

        double initialenergy = energy;

        VectorXd newq = q + searchdir;
        energy = no_.getEnergy(newq);
        double neededeps=0;
        double eps = params.lmfactor;
        int lsiters=0;

        while(energy > initialenergy || isnan(energy))
        {
            neededeps = eps;
            if(++lsiters > params.lsmaxiters)
            {
                std::cout << std::endl;
                return LSITERS_EXCEEDED;
            }
            SparseMatrix<double> shift(hessian.rows(), hessian.cols());
            std::vector<Tr> shiftcoeffs;
            for(int i=0; i<hessian.rows(); i++)
                shiftcoeffs.push_back(Tr(i,i,eps*hessian.coeffRef(i,i)));
            shift.setFromTriplets(shiftcoeffs.begin(), shiftcoeffs.end());
            SparseMatrix<double> newh = hessian + shift;
            SimplicialLDLT<SparseMatrix<double> > solver;
            solver.compute(newh);
            searchdir = -solver.solve(gradient);
            newq = q + searchdir;
            residual = (newh*searchdir+gradient).norm();
            eps *= 10;
            energy = no_.getEnergy(newq);
            //std::cerr << residual << " " << searchdir.norm() << " " << energy << " " << initialenergy << std::endl;
        }

        std::cout << std::fixed << std::setprecision(8) << "   lm " << neededeps << std::endl;

        q = newq;
        no_.showCurrentIteration(q);
    }

    if(iter == params.maxiters)
        return MAXITERS_EXCEEDED;

    result = q;

    return CONVERGED;
}

std::string Newton::solverStatusMessage(SolverStatus status) const
{
    switch(status)
    {
    case CONVERGED:
        return std::string("Converged");
    case MAXITERS_EXCEEDED:
        return std::string("Maximum outer iterations exceeded");
    case LSITERS_EXCEEDED:
        return std::string("Maximum line search iterations exceeded");
    case BAD_INPUT:
        return std::string("Bad function data supplied by objective functor");
    case NONE:
    default:
        return std::string("Invalid Solver Status");
    }
}
