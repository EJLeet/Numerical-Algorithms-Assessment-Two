#include <cmath>
#include <iostream>

using std::endl;
using std::cout;

// set function to 0 and use bisection
double f(double theta)
{
    return std::tan(theta) * 35 - (9.81 / (2 * 20 * 20 *
        (std::cos(theta) * std::cos(theta)))) * 35 * 35 + 2 - 1;
}

double bisection(double a, double b);

int main()
{   
    double theta1 = bisection(0, M_PI / 2);
    double theta2 = bisection(theta1 * M_PI / 180, M_PI / 2);
    cout << "Theta 1 = " << theta1 << " degrees" << endl
         << "Theta 2 = " << theta2 << " degrees" << endl;

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