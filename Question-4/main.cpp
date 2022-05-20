#include <iostream>
#include <vector>
#include <cmath>

using std::endl;
using std::cout;

void multiple_linear_regression(std::vector< std::vector<double>>& lhs, 
                                std::vector<double>& rhs, std::vector<double>& d, 
                                std::vector<double>& s, std::vector<double>& f);
std::vector<double> gaussian_elimination(std::vector<std::vector<double>>& lhs, 
                                         std::vector<double>& rhs);
// used to fill multiple matrix values                                       
void fill_matrix(std::vector<double>& row, const std::vector<double>& inc)
{ for (int i = 0; i < inc.size(); i++) row[i] += inc[i]; }
double F(double a, double b, double c, double d, double s)
{ return a * pow(d, b) * pow(s, c); }

int main()
{
    std::vector<std::vector<double>> lhs = {{0, 0, 0},
                                            {0, 0, 0},
                                            {0, 0, 0}};
    std::vector<double> d = {1, 2, 3, 1, 2, 3, 1, 2, 3},
                        s = {0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05},
                        f = {1.4, 8.3, 24.2, 4.7, 28.9, 84.0, 11.1, 69.0, 200.0},
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
             <<  F(a, b, c, d[i], s[i]) << endl;
    
    // solve f at given values for d and s
    cout << "F = " << F(a, b, c, 2.5, 0.025) << " at D = 2.5 and S = 0.025" << endl;

    return 0;
}

void multiple_linear_regression(std::vector< std::vector<double>>& lhs, 
                                std::vector<double>& rhs, std::vector<double>& d, 
                                std::vector<double>& s, std::vector<double>& f)
{
    double x1 = 0, x2 = 0, y = 0;

    lhs[0][0] = d.size();
    for (int i = 0; i < d.size(); i++)
    {
        // get x1, x2 and y values
        x1 = log10(d[i]);
        x2 = log10(s[i]);
        y = log10(f[i]);

        // fill multiple linear regression matrix
        fill_matrix(lhs[0], {0, x1, x2});
        fill_matrix(lhs[1], {x1, x1 * x1, x1 * x2});
        fill_matrix(lhs[2], {x2, x1 * x2, x2 * x2});

        // fill y matrix
        fill_matrix(rhs, {y, y * x1, y * x2});
    }
}

std::vector<double> gaussian_elimination(std::vector<std::vector<double>>& lhs, 
                                         std::vector<double>& rhs)
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