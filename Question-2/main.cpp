#include <cmath>
#include <iostream>

using std::endl;
using std::cout;

// set fuinction to 0 and use secant
double f(double theta)
{
    return std::tan(theta) * 35 - (9.81 / (2 * 20 * 20 *
        (std::cos(theta) * std::cos(theta)))) * 35 * 35 + 2 - 1;
}

double bisection(double a, double b);

int main()
{
    cout << "Theta = " << bisection(0, M_PI / 2) << " degrees" << endl;

    return 0;
}

double bisection(double a, double b)
{
    double c;

    while (fabs(b - a) >= 0.00001)
    {
        c = (a + b) / 2;

        if (f(c) == 0.0)
            break;
        else if ((f(c) * f(a)) < 0)
            b = c;
        else
            a = c;
    }
    return c * 180 / M_PI;
}