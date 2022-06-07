#include "lbfgs.hpp" // https://github.com/ZJU-FAST-Lab/LBFGS-Lite
#include "Eigen/Eigen" // https://eigen.tuxfamily.org/index.php?title=Main_Page
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
int N, P;

void set_position(Eigen::VectorXd& x);

class MinimizationExample
{
public:
    std::vector<double> minima_count;
    double best = 0;
    Eigen::VectorXd global_config;
    int run(Eigen::VectorXd& x)
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

        // count distinct minima
        auto it = std::find_if(minima_count.begin(), minima_count.end(), 
            [finalCost](double b) { return fabs(finalCost - b) < 1e-6; });
        if (it == minima_count.end())
            minima_count.push_back(finalCost);

        // get global config and cost
        if (finalCost < best)
        {
            best = finalCost;
            global_config = x.transpose();
        }
        
        return ret;
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
        return 0;
    }
};

int main()
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
                    */
    N = 2;
    cout << "PART B" << endl;
    double minima;
    MinimizationExample part_b;

    for (auto i : p)
    {
        Eigen::VectorXd x = Eigen::VectorXd::Zero(N * 3);
        set_position(x);
        P = i;
        part_b.run(x);

        cout << "minima at n(" << N << "), p(" << 
            P << ") = " << part_b.best << endl;
    }

/*
        Part C - A
                    */

    cout << endl << endl << "Part C - A" << endl;
     int iter = 100000;
    std::ofstream c_a_results("c_a_results.txt");
    for (int n = 2; n <= 32; n++)
    {
        N = n;
        c_a_results << N;
        MinimizationExample part_ca;
        for (auto pval : p)
        {
            P = pval;

            // loops for trials changing atom positions per iter
            for (int t = 0; t < iter; t++)
            {
                if (t % 1000 == 0) 
                    cout << "Progress Part C - A, N = " << N << ",P = " << P << ": " << t << "/" << iter << endl;
                
                Eigen::VectorXd x2 = Eigen::VectorXd::Zero(N * 3);
                set_position(x2);

                part_ca.run(x2);
            }

            // save results to file
            c_a_results << " " << part_ca.minima_count.size();
            part_ca.minima_count.clear();
        }
        c_a_results << endl;
    }

    /*
        Part C - B
                    */

    cout << endl << endl << "Part C - B" << endl;
    std::ofstream c_b_results("c_b_results.txt");
    P = 6;
    for (int n = 2; n <= 32; n++)
    {
        MinimizationExample part_cb;
        N = n;

        // loops for trials changing atom positions per iter
        for (int t = 0; t < iter; t++)
        {
            if (t % 1000 == 0) 
                cout << "Progress Part C - B, N = " << N << ",P = 6: " <<  t << "/" << iter << endl;
            
            Eigen::VectorXd x3 = Eigen::VectorXd::Zero(N * 3);
            set_position(x3);

            part_cb.run(x3);
        }
        
        // save results to file
        c_b_results << N << " " << part_cb.best << " " << endl << part_cb.global_config << endl << endl;
        
    }

    return 0;
}

void set_position(Eigen::VectorXd& x)
{// set rand atomic positions based on n^1/3
    double d = 2.2 * sqrt(3) * (pow(N, 0.333) - 1) * 
        (0.8 + double(rand()) / double(RAND_MAX));

    for (int i = 0; i < N * 3; i++) 
        x[i] = (double(rand()) / double(RAND_MAX) - 0.5) * d;
}