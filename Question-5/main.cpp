#include <iostream>
#include <cmath>
#include <fstream>

using std::cout;
using std::endl;

// function will return the length of the rod at some theta and alpha
double l(double w1, double w2, double theta, double alpha)
{ return -((w1 / std::sin(theta)) + (w2 / (std::sin(M_PI - alpha - theta)))); }

// perform golden section search
double golden_section(double w1, double w2, double alpha);

int main()
{
    double w1 = 2, w2 = 2;
    std::ofstream file("output.txt");
    for (double alpha = 0.0; alpha < M_PI; alpha += 0.01)
    {
        
        double temp = golden_section(w1, w2, alpha);
        double length = l(w1, w2, temp, alpha);

        cout << "Rod Length = " << length * -1 << "m"
             << " at alpha = " << alpha * 180 / M_PI << " degrees" << endl;
        file << length << " " << alpha << endl;
    }
    file.close();
    return 0;
}

double golden_section(double w1, double w2, double alpha)
{/*
    This function performs the golden section search for the
    minimum of theta at some alpha.
                                                                */
    // initialise variables
    const double ratio = 0.61803;
    double xl = 0, xu = M_PI - alpha;
    double x1 = xl + ratio * (xu - xl), x2 = x2 = xu - ratio * (xu - xl);
    double fx1 = l(w1, w2, x1, alpha), fx2 = l(w1, w2, x2, alpha), 
           fxl = l(w1, w2, xl, alpha), fxu = l(w1, w2, xu, alpha);

    while (fabs(xu - xl) > 0.00001)
    {
        if (fx1 > fx2)
        {
            xl = x2;
            fxl = fx2;
            x2 = x1;
            fx2 = fx1;     
            x1 = xl + ratio * (xu - xl);
            fx1 = l(w1, w2, x1, alpha);
        }
        else
        {
            xu = x1;
            fxu = fx1;
            x1 = x2;
            fx1 = fx2;    
            x2 = xu - ratio * (xu - xl);
            fx2 = l(w1, w2, x2, alpha);
        }   
    }  
    return xl;  
    }
