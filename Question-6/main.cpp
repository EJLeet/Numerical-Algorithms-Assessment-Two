#include <iostream>
#include <cmath>
#include <fstream>
#include <array>
#include "Eigen/Core"
#include "LBFGS.h"
#include "LBFGSB.h"

using std::cout;
using std::endl;
using Eigen::VectorXd;
using namespace LBFGSpp;

double f(double ep)
{ return (ep - 1) * (ep - 1) - 1; }

class Energy
{
    const int n;
    int p;

public:
    Energy(int n_, int p_) : n(n_), p(p_)  {}
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double V = 0;

        // initialise gradient vector to zero
        for (int i = 0; i < 3 * n; i++)
            grad[i] = 0;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                // update energy
                double disx = x[3 * i] - x[3 * j];
                double disy = x[3 * i + 1] - x[3 * j + 1];
                double disz = x[3 * i + 2] - x[3 * j + 2];
                double r = sqrt(disx * disx + disy * disy + disz * disz);
                double ep = exp(p * (r - 1));
                double v = f(ep);
                V += v;

                // get gradients
                double deriv = 2 * p * ep * (ep - 1);
                double dx = deriv * disx / r;
                double dy = deriv * disy / r;
                double dz = deriv * disz / r;

                // update positive for one grad negative for other
                grad[i * 3]     += dx;
                grad[i * 3 + 1] += dy;
                grad[i * 3 + 2] += dz;
                grad[j * 3]     -= dx;
                grad[j * 3 + 1] -= dy;
                grad[j * 3 + 2] -= dz;
            
            }
        }
        return V;
    }
};

void find_minima(VectorXd x, const int n, int p);

int main()
{
    double ep;
    const int n = 2;
    std::array<int, 4> p = {3, 6, 10, 14};

    // outfiles for part a graphs
    std::ofstream p3("3.txt");
    std::ofstream p6("6.txt");
    std::ofstream p10("10.txt");
    std::ofstream p14("14.txt");

    // traverse each p value and all r values to create graphs
    for (auto i : p)
    {
        for (double r = 0.6; r <= 2.0; r += 0.001)
        {
            ep = exp(i * (r - 1));
            switch(i)
            {
                case 3 : { p3  << r << " " << f(ep) << endl; break; }
                case 6 : { p6  << r << " " << f(ep) << endl; break; }
                case 10: { p10 << r << " " << f(ep) << endl; break; }
                case 14: { p14 << r << " " << f(ep) << endl; break; }
            }
        }
    }

    /*
        Find energy configuration between two 
        atoms for p = 3, 6, 10
                                            */
    // Initial atoms decleration
    VectorXd x = VectorXd::Zero(3 * n);

    // first particle xyz
    x[0] = 0.8;
    x[1] = -0.1;
    x[2] = 18;

    // second particle xyz
    x[3] = 0.7;
    x[4] = -1.1;
    x[5] = -0.3;

    // find minima across different p
    // exclude p = 14 due to p = 14 being too large
    for (int i = 0; i < p.size() - 1; i++)
        find_minima(x, n, p[i]);

    return 0;
}


void find_minima(VectorXd x, const int n, int p)
{
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 10000000;

    // Create solver and function object
    LBFGSSolver<double> solver(param);
    Energy fun(n, p);

    // fx will be overwritten to minima found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    cout << "p = " << p << endl;
    cout << niter << " iterations" << endl;
    cout << "New xyz for particles = " << x.transpose() << endl;
    cout << "Minima of f(x) = " << fx << endl << endl;
}

