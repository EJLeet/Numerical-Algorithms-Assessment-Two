#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

// function will return the length of the rod at some theta and alpha
double l(double w1, double w2, double theta, double alpha)
{ return ((w1 / std::sin(theta)) + (w2 / (std::sin(M_PI - alpha - theta)))); }

// perform golden section search
double golden_section(double w1, double w2, double alpha);

int main()
{
    double w1 = 2, w2 = 2;
    for (double alpha = M_PI / 2; alpha < M_PI; alpha += 0.01)
    {
        
        double temp = golden_section(w1, w2, alpha);
        double length = l(w1, w2, temp, alpha);

        cout << "Rod Length = " << length
             << " at alpha = " << alpha * 180 / M_PI << endl;
    }
    return 0;
}

double golden_section(double w1, double w2, double alpha)
{
    double ratio = 0.61803, xl = 0, xu = M_PI - alpha;
    double d = ratio * (xu - 0);
    double x1 = xl + d, x2 = xu - d;
    int count = 0;

    while (fabs((x1 - x2)) > 0.00001)
    {
        if (l(w1, w2, x1, alpha) > l(w1, w2, x2, alpha))
        {
            x2 = x1;
            x1 = xl + d;
        }
        else
        {
            x1 = x2;
            x2 = xu - d;
        }
    }
    return x1;
}
