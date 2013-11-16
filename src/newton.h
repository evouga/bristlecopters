#ifndef NEWTON_H
#define NEWTON_H

#include <Eigen/Core>
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> Tr;

class NewtonObjective
{
public:
    virtual ~NewtonObjective() {}
    virtual double getEnergy(const Eigen::VectorXd &q) const = 0;
    virtual double getEnergyAndDerivatives(const Eigen::VectorXd &q, Eigen::VectorXd &grad, Eigen::SparseMatrix<double> &hess) const = 0;
    virtual void showCurrentIteration(const Eigen::VectorXd &) const {}
};

struct NewtonParameters
{
    NewtonParameters() : tol(1e-6), maxiters(100), lsmaxiters(100), lmfactor(1e-9) {}
    double tol;
    int maxiters;
    int lsmaxiters;
    double lmfactor;
};

class Newton
{
public:
    enum SolverStatus {
        NONE,
        CONVERGED,
        MAXITERS_EXCEEDED,
        LSITERS_EXCEEDED,
        BAD_INPUT
    };

    Newton(const NewtonObjective &no);
    SolverStatus solve(NewtonParameters params, const Eigen::VectorXd &guess, Eigen::VectorXd &result);
    std::string solverStatusMessage(SolverStatus status) const;

private:
    const NewtonObjective &no_;
};

#endif // NEWTON_H
