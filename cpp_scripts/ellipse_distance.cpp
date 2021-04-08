#include <iostream>
#include <utility>
#include <stdlib.h>
#include <math.h>


double ellipse_distance_sqr(const double a, const double b, const double xp, 
                        const double yp, const double theta)
{
    const double distance_sqr = pow(a * cos(theta) - xp, 2.0) +
                                 pow(b * sin(theta) - yp, 2.0);
    return distance_sqr;
}


double ellipse_distance_sqr_div(const double a, const double b, const double xp, 
                                const double yp, const double theta)
{
    const double distance_sqr_div = (b*b - a*a) * sin(2.0*theta)
                                    + 2.0 * a * xp * sin(theta) 
                                    - 2.0 * b * yp * cos(theta);
    // std::cout << "a: " << a << ", b: " << b << ", xp: " << xp << ", yp: " << yp << ", theta: " << theta << std::endl;  // debug
    // std::cout << "distance_sqr_div: " << distance_sqr_div << std::endl;
    return distance_sqr_div;
}


double ellipse_distance_sqr_div_div(const double a, const double b, const double xp, 
                                    const double yp, const double theta)
{
    const double distance_sqr_div_div = 2.0 * (b*b - a*a) * cos(2*theta)
                                    + 2.0 * a * xp * cos(theta)
                                    + 2.0 * b * yp * sin(theta);
    return distance_sqr_div_div;
}


// Newton method
std::pair<double, int> 
ellipse_distance_sqr_shortest(const double a, const double b,
                              const double xp, const double yp, const double tolerance)
{
    double theta = M_PI / 4.0;
    double distance_sqr_div = ellipse_distance_sqr_div(a, b, xp, 
                                                       yp, theta);

    // use a newtonian method, looking for minimum derivative of distance squared
    // in a loop.
    int i = 0;
    while(abs(distance_sqr_div / (a * a)) > tolerance)
    {
        double distance_sqr_div_div = ellipse_distance_sqr_div_div(a, b, xp, 
                                                                   yp, theta);
        distance_sqr_div = ellipse_distance_sqr_div(a, b, xp, 
                                                    yp, theta);
        double theta_diff = - distance_sqr_div / distance_sqr_div_div;
        theta += theta_diff; 
        distance_sqr_div = ellipse_distance_sqr_div(a, b, xp, 
                                                    yp, theta);
        i++;
    }
    
    // return both the result and the number of iteration 
    return std::make_pair(theta, i);
}

// dual method
std::pair<double, int> 
ellipse_distance_sqr_shortest1(const double a, const double b,
                              const double xp, const double yp, const double tolerance)
{
    double theta0 = 0.0;
    double theta1 = M_PI / 2.0;
    double distance_sqr_div0 = ellipse_distance_sqr_div(a, b, xp, 
                                                       yp, theta0);
    double distance_sqr_div1 = ellipse_distance_sqr_div(a, b, xp, 
                                                       yp, theta1);
    
    int i = 0;
    while(abs(distance_sqr_div0 / (a * a)) > tolerance)
    {
       double theta2 = (theta0 + theta1) / 2.0; 
       double distance_sqr_div2 = ellipse_distance_sqr_div(a, b, xp, 
                                                           yp, theta2);
        if(distance_sqr_div2 > 0.0)
        {
            theta1 = theta2;
            distance_sqr_div1 = distance_sqr_div2;
        }
        else if(distance_sqr_div2 < 0.0)
        {
            theta0 = theta2;
            distance_sqr_div0 = distance_sqr_div2;
        }
        else
        {
            // 0.0, break and return
            theta0 = theta2;
            break;
        }
        i++;
    }

    // return both the result and the number of iteration 
    return std::make_pair(theta0, i);
}

void sample_outputs(const double a, const double b,
                    const double xp, const double yp)
{
    // sample outputs
    double theta = 0.0; 
    double distance = sqrt(ellipse_distance_sqr(a, b, xp, yp, theta));
    std::cout << "theta: " << theta << ", distance: " << distance << std::endl;
    
    theta = M_PI / 2.0; 
    distance = sqrt(ellipse_distance_sqr(a, b, xp, yp, theta));
    std::cout << "theta: " << theta << ", distance: " << distance << std::endl;
    
    theta = -3.22; 
    distance = sqrt(ellipse_distance_sqr(a, b, xp, yp, theta));
    std::cout << "theta: " << theta << ", distance: " << distance << std::endl;
    
    theta = -M_PI / 2.0; 
    distance = sqrt(ellipse_distance_sqr(a, b, xp, yp, theta));
    std::cout << "theta: " << theta << ", distance: " << distance << std::endl;
    
    theta = -3.0; 
    distance = sqrt(ellipse_distance_sqr(a, b, xp, yp, theta));
    std::cout << "theta: " << theta << ", distance: " << distance << std::endl;
    
    theta = -3.3; 
    distance = sqrt(ellipse_distance_sqr(a, b, xp, yp, theta));
    std::cout << "theta: " << theta << ", distance: " << distance << std::endl;

}


int main(int argc, char *argv[]){
    // we need the semi-axis of the ellipse, and the coordinates of the query point
    if(argc != 5){
        std::cerr << "Four input values are needed" << std::endl;
        return 1;
    }
    const double a = atof(argv[1]);
    const double b = atof(argv[2]);
    const double xp = atof(argv[3]);
    const double yp = atof(argv[4]);
    
    //const double distance = ellipse_distance_sqr(a, b, xp, yp, M_PI/5.0);
    double tolerance = 1e-5;
    std::pair<double, int> results = ellipse_distance_sqr_shortest1(a, b,
                                                                   xp, yp, tolerance);
    
    double distance_shortest = sqrt(ellipse_distance_sqr(a, b, xp, 
                                                        yp, results.first));
    
    std::cout << results.first << std::endl << distance_shortest 
              << std::endl << results.second << std::endl;

    //const double distance = ellipse_distance_sqr(a, b, xp, yp, M_PI/5.0);
    // sample_outputs(a, b, xp, yp);
    return 0;
}