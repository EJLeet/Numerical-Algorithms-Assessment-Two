#include <iostream>
#include <cmath>
#include <vector>
#define x6 10 * 10 * 10 * 10 * 10 * 10

using std::cout;
using std::endl;

double moment(double wb, double theta);
double secant(double wb);
double interpolate(std::vector<double> x, std::vector<double> y, double xv);
void linear_regression(std::vector<double> x, std::vector<double> y);

int main()
{
    // part b
    double wb = 1.663 * x6, fuel = 1.300 * x6;
    double fuel_part = fuel / 127;

    cout << "Root of the moment equation: Theta (radians) = " << secant(wb) << endl;

    // part c
    std::vector<double> x, y;

    for (int i = 0; i < 128; i++)
    {
        x.push_back(i);
        y.push_back(secant(wb - fuel_part * i));
    }

    cout << endl;

    // do lagrange with first and last point - linear
    std::vector<double> x_ = {0, 127}, y_ = {y[0], y[127]};
    double temp = interpolate(x_, y_, 63);
    cout << "Lagrange estimate on point 63 using points {0, 127}:\n"
         << "Predicted Theta (rads): " << temp << "\nTrue Theta (rads): " << y[100] 
         << "\nTrue error: " << abs((y[63] - temp) / y[63]) << "%" << endl << endl;
    
    // linear regression
    cout << "Linear Regression 2 Points: " << endl;
    linear_regression(x_, y_);
    cout << endl;

    // do legrange with first last and middle point - quadratic
    std::vector<double> x__ = {0, 63, 127}, y__ = {y[0], y[63], y[127]};
    temp = interpolate(x__, y__, 100);
    cout << "Lagrange estimate on point 100 using points {0, 63, 127}:\n"
         << "Predicted Theta (rads): " << temp << "\nTrue Theta (rads): " << y[100] 
         << "\nTrue error: " << abs((y[100] - temp) / y[100]) << "%" << endl << endl;
    
    // linear regression
    cout << "Linear Regression 3 Points: " << endl;
    linear_regression(x__, y__);
    cout << endl;


    // lagrange at 4 points for cubic
    std::vector<double> x___ = {0, 42, 84, 127}, y___ = {y[0], y[42], y[84], y[127]};
    temp = interpolate(x___, y___, 100);
    cout << "Lagrange estimate on point 100 using points {0, 42, 84, 127}:\n"
         << "Predicted Theta (rads): " << temp << "\nTrue Theta (rads): " << y[100] 
         << "\nTrue error: " << abs((y[100] - temp) / y[100]) << "%" << endl << endl;
    
    // linear regression
    cout << "Linear Regression 4 Points: " << endl;
    linear_regression(x___, y___);
    cout << endl;

    // lagrange at 5 points for pow(4)
    std::vector<double> x____ = {0, 32, 64, 96, 127}, y____ = {y[0], y[32], y[64], y[96], y[127]};
    temp = interpolate(x____, y____, 100);
    cout << "Lagrange estimate on point 100 using points {0, 32, 64, 96, 127}:\n"
         << "Predicted Theta (rads): " << temp << "\nTrue Theta (rads): " << y[100] 
         << "\nTrue error: " << abs((y[100] - temp) / y[100]) << "%" << endl << endl;
    
    // linear regression
    cout << "Linear Regression 5 Points: " << endl;
    linear_regression(x____, y____);
    cout << endl;

    return 0;
}

double moment(double wb, double theta)
{
    double ws = 0.23 * x6, tb = 5.30 * x6, ts = 1.125 * x6;
    return (tb * 4 + ws * 24 + 38 * ts * std::sin(theta)) - 
           (wb * 4 + 24 * ts * std::cos(theta));
}

double secant(double wb)
{
    double a = 0, b = M_PI / 2;
    double c;
    while (fabs(b - a) >= 0.0005)
    {
        c = b - (moment(wb, b) * (b - a)) / (moment(wb, b) - moment(wb, a));
        a = b;
        b = c;
    }
    return c;
}

double interpolate(std::vector<double> x, std::vector<double> y, double xv)
{
    double result = 0;
    for (int i = 0; i < x.size(); i++)
    {
        double p = 1;
        for (int j = 0; j < x.size(); j++)
            if (j != i)
                p *= (xv - x[j]) / (x[i] - x[j]);
        result += p * y[i];
    }
    return result;
}

void linear_regression(std::vector<double> x, std::vector<double> y)
{
    double a0, a1, sumxy = 0, sumx = 0, sumy = 0, sumx2 = 0;
    
    for (int i = 0; i < x.size(); i++)
    {
        sumxy += x[i] * y[i];
        sumx += x[i];
        sumy += y[i];
        sumx2 += x[i] * x[i];
    }
    a1 = (x.size() * sumxy - sumx * sumy) / (x.size() * sumx2 - sumx * sumx);
    a0 = (sumy / x.size()) - a1 * (sumx / x.size());
    
    cout << "y = " << a1 << "x + " << a0 << endl;

    double st = 0, ybar = sumy / x.size();
    for (int i = 0; i < x.size(); i++)
        st += (y[i] - ybar) * (y[i] - ybar);

    double sr = 0;
    for (int i = 1; i < x.size(); i++)
        sr += (y[i] - a0 - a1 * x[i]) * (y[i] - a0 - a1 * x[i]);
        
    cout << "st = " << st << ", sr = " << sr << endl;
    double r2 = (st - sr) / st;
    cout << "r2 = " << r2 << ", Correlation coeff = " << sqrt(r2) << endl;
}
