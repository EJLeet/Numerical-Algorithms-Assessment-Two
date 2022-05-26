#include <iostream>
#include <cmath>
#include <fstream>
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
private:
    int n, p;

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
                double posx = x[3 * i] - x[3 * j];
                double posy = x[3 * i + 1] - x[3 * j + 1];
                double posz = x[3 * i + 2] - x[3 * j + 2];
                double r = sqrt(posx * posx + posy * posy + posz * posz);
                double ep = exp(p * (r - 1));
                double v = f(ep);
                V += v;

                // get gradients
                double deriv = 2 * p * ep * (ep - 1);
                double dx = deriv * posx / r;
                double dy = deriv * posy / r;
                double dz = deriv * posz / r;

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

void part_a();
void part_b();


int main()
{
    part_a();
    part_b();
    return 0;
}

void part_a()
{// get function and r values for 5 constant p vals

    double ep;
    std::ofstream p3("3.txt");
    for (double r = 0.6; r <= 2.0; r += 0.001)
    {
        ep = exp(3 * (r - 1));
        p3 << r << " " << f(ep) << endl;
    }

    std::ofstream p6("6.txt");
    for (double r = 0.6; r <= 2.0; r += 0.001)
    {
        ep = exp(6 * (r - 1));
        p6 << r << " " << f(ep) << endl;
    }

    std::ofstream p10("10.txt");
    for (double r = 0.6; r <= 2.0; r += 0.001)
    {
        ep = exp(10 * (r - 1));
        p10 << r << " " << f(ep) << endl;
    }

    std::ofstream p14("14.txt");
    for (double r = 0.6; r <= 2.0; r += 0.001)
    {
        ep = exp(14 * (r - 1));
        p14 << r << " " << f(ep) << endl;
    }
}

void part_b()
{
    const int n = 2;
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 10000000;

    // Create solver and function object
    LBFGSSolver<double> solver(param);
    Energy fun(n, 3);

    // Initial guess
    VectorXd x = VectorXd::Zero(3 * n);

    // first particle xyz
    x[0] = 1.8;
    x[1] = 2.8;
    x[2] = 3.8;

    // second particle xyz
    x[3] = 4;
    x[4] = -1;
    x[5] = 18;

    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "New xyz for particles = " << x.transpose() << std::endl;
    std::cout << "Minima of f(x) = " << fx << std::endl;
}