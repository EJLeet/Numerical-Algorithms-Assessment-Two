#include "lbfgs.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <random>
#include <fstream>

using std::cout;
using std::endl;

double f(double ep)
{ return (ep - 1) * (ep - 1) - 1; }
int N = 10, P = 6;

void set_position(Eigen::VectorXd& x);

class MinimizationExample
{
public:
    std::tuple<Eigen::VectorXd, int> run(Eigen::VectorXd& x)
    {
        double finalCost;

        /* Set the minimization parameters */
        lbfgs::lbfgs_parameter_t params;
        params.g_epsilon = 1.0e-8;
        params.past = 3;
        params.delta = 1.0e-8;

        /* Start minimization */
        int ret = lbfgs::lbfgs_optimize(x,
                                        finalCost,
                                        costFunction,
                                        nullptr,
                                        monitorProgress,
                                        this,
                                        params);
        /* Report the result. */
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "L-BFGS Optimization Returned: " << ret << std::endl
                  << "Minimized Cost: " << finalCost << std::endl
                  << "Optimal Variables: " << std::endl
                  << x.transpose() << std::endl;
        
        
        
        return {x, finalCost};
    }

private:
    static double costFunction(void *instance,
                               const Eigen::VectorXd &x,
                               Eigen::VectorXd &grad)
    {
        double V = 0;

        // initialise gradient vector to zero
        for (int i = 0; i < 3 * N; i++)
            grad[i] = 0;

        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {
                // update energy
                double disx = x[3 * i] - x[3 * j];
                double disy = x[3 * i + 1] - x[3 * j + 1];
                double disz = x[3 * i + 2] - x[3 * j + 2];
                double r = sqrt(disx * disx + disy * disy + disz * disz);
                double ep = exp(P * (1 - r));
                double v = f(ep);
                V += v;

                // get gradients
                double deriv = -2 * P * ep * (ep - 1);
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

    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        g.cwiseAbs().maxCoeff();
        x.transpose();

        // std::cout << std::setprecision(4)
        //           << "================================" << std::endl
        //           << "Iteration: " << k << std::endl
        //           << "Function Value: " << fx << std::endl
        //           << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
        //           << "Variables: " << std::endl
        //           << x.transpose() << std::endl;
        return 0;
    }
};

int main(int argc, char **argv)
{
    /*
        Part A
                */
    double ep;
    std::array<int, 4> p = {3, 6, 10, 14};
    // outfiles for part a graphs
    std::ofstream p3("3.txt");
    std::ofstream p6("6.txt");
    std::ofstream p10("10.txt");
    std::ofstream p14("14.txt");

    for (auto i : p)
    {
        for (double r = 0.6; r <= 2.0; r += 0.001)
        {
            ep = exp(i * (1 - r));
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
        Part B
        Set global P values
                                */
    double minima;
    MinimizationExample part_b;
    cout << "n = " << N << ", p = " << P << endl;
    Eigen::VectorXd x = Eigen::VectorXd::Zero(N * 3), x_;
    set_position(x);
    
    cout << "Original particle position (x y z): ";
    int count = 0;
    for (auto i : x) 
    {
        if (count++ % 3 == 0) cout << endl;
        cout << i << "    ";
    }
    cout << endl;

    std::tie(x_, minima) = part_b.run(x);

    cout << "Minima = " << minima << endl << "New particle position (x y z): ";
    count = 0;
    for (auto i : x_) 
    {
        if (count++ % 3 == 0) cout << endl;
        cout << i << "    ";
    }
    cout << endl << endl;

/*
        Part C
        Set global N/P values
                                */
    std::vector<double> minima_count;
    for (int t = 0; t < 1000; t++)
    {
        MinimizationExample part_c;
        Eigen::VectorXd x2 = Eigen::VectorXd::Zero(N * 3), x2_;
        set_position(x2);

        std::tie(x2_, minima) = part_c.run(x2);
        cout << "Minima Found = " << minima << endl;

        // add minima to set of minima if not already found
        auto it = std::find_if(minima_count.begin(), minima_count.end(), 
            [minima](double b) { return fabs(minima - b) < 1e-6; });
        if (it == minima_count.end())
            minima_count.push_back(minima);
    }
    cout << "Count of distinc minima at n(" << N << "), p(" << 
        P << ") = " << minima_count.size() << endl;

    return 0;
   
}

void set_position(Eigen::VectorXd& x)
{
    // uniform random number gen
    // std::random_device rand_dev;
    // std::mt19937 generator(rand_dev());

    // for (int i = 0; i < N * 3; i++)
    // {
    //     std::uniform_real_distribution<double> gen(0, cbrt(N));
    //     double random = gen(generator);
    //     x[i] = random;
    // }
    double d = 2.2 * sqrt(3) * (pow(N, 0.333) - 1) * 
        (0.8 + double(rand()) / double(RAND_MAX));
    for (int i = 0; i < N * 3; i++)
        x[i] = (double(rand()) / double(RAND_MAX) - 0.5) * d;
}