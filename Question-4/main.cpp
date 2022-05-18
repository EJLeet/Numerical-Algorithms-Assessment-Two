#include <iostream>
#include <vector>
#include <cmath>

using std::endl;
using std::cout;

void multiple_linear_regression(std::vector< std::vector<double>>& lhs, 
                                std::vector<double>& rhs, std::vector<double>& d, 
                                std::vector<double>& s, std::vector<double>& f);
std::vector<double> gaussian_elimination(std::vector<std::vector<double>> lhs, 
                                         std::vector<double> rhs);
double F(double a, double b, double c, double d, double s)
{ return a * pow(d, b) * pow(s, c); }
double F_log(double a, double b, double c, double d, double s)
{ return log10(a) + b * d + c * s; }

int main()
{
    std::vector<std::vector<double>> lhs = {{0, 0, 0},
                                            {0, 0, 0},
                                            {0, 0, 0}};
    std::vector<double> d = {log10(1), log10(2), log10(3), log10(1), log10(2), log10(3), log10(1), log10(2), log10(3)},
                        s = {log10(0.001), log10(0.001), log10(0.001), log10(0.01), log10(0.01), log10(0.01), log10(0.05), log10(0.05), log10(0.05)},
                        f = {log10(1.4), log10(8.3), log10(24.2), log10(4.7), log10(28.9), log10(84.0), log10(11.1), log10(69.0), log10(200.0)},
                        rhs = {0, 0, 0}, vals;
                                            
    
    // get matrix values
    multiple_linear_regression(lhs, rhs, d, s, f);
    
    // find a/b/c using matrix values
    vals = gaussian_elimination(lhs, rhs);

    double a = pow(10, vals[0]),
           b = vals[1],
           c = vals[2];
    cout << "A  " << a << "   B " <<  b << "    C " << c << endl;

    // prove gaussian is correct
    for (int i = 0; i  < 9; i++)
        cout << "Confirming Experiment " << i + 1 << " = " 
             <<  pow(10, F_log(a, b, c, d[i], s[i])) << endl;
    
    // solve f at given values for d and s
    cout << "F = " << F(a, b, c, 2.5, 0.025) << " at D = 2.5 and S = 0.025" << endl;

    return 0;
}

void multiple_linear_regression(std::vector< std::vector<double>>& lhs, 
                                std::vector<double>& rhs, std::vector<double>& d, 
                                std::vector<double>& s, std::vector<double>& f)
{
    double x_1 = 0, x_2 = 0, x2_1 = 0, x2_2 = 0, 
           x12 = 0, y = 0, yx1 = 0, yx2 = 0;

    lhs[0][0] = d.size();
    for (int i = 0; i < d.size(); i++)
    {
        // get x1, x2 and y values
        x_1 = d[i];
        x2_1 = x_1 * x_1;
        x_2 = s[i];
        x2_2 = x_2 * x_2;
        x12 = x_1 * x_2;
        y = f[i];
        yx1 = y * x_1;
        yx2 = y * x_2;

        // fill multiple linear regression matrix
        lhs[0][1] += x_1;
        lhs[0][2] += x_2;
        lhs[1][0] += x_1;
        lhs[1][1] += x2_1;
        lhs[1][2] += x12;
        lhs[2][0] += x_2;
        lhs[2][1] += x12;
        lhs[2][2] += x2_2;

        // fill y matrix
        rhs[0] += y;
        rhs[1] += yx1;
        rhs[2] += yx2;
    }
}

std::vector<double> gaussian_elimination(std::vector<std::vector<double>> lhs, 
                                         std::vector<double> rhs)
{
    double sum = 0;

    // forward elimination
    for (int k = 0; k < rhs.size() - 1; k++)
        for (int i = k + 1; i < rhs.size(); i++)
        {
            float s = lhs[i][k] / lhs[k][k];
            for (int j = k; j < rhs.size(); j++)
                lhs[i][j] -= lhs[k][j] * s;
            rhs[i] = rhs[i] - rhs[k] * s;
        }

    // backward substitution
    rhs[2] = rhs[2] / lhs[2][2];
    for (int i = rhs.size() - 2; i >= 0; i--)
    {
        sum = rhs[i];
        for (int j = i+1; j < rhs.size(); j++)
            sum -=  lhs[i][j] * rhs[j];
        rhs[i] = sum / lhs[i][i];
    }
    return rhs;
}