#include "function.h"

Function::Function() { }

Function::~Function() { }

double Function::evaluate(std::vector<double>& x) {
    // return std::pow(x[0], 3.0) + std::pow(x[0], 2.0) * x[1] + 
    //        x[0] * std::pow(x[1], 2.0) + std::pow(x[1], 3.0) +
    //        x[0] * x[1] * x[2] + std::pow(x[2], 3.0) + 
    //        x[1] * std::pow(x[2], 2.0);

    return 0.5 * std::pow((x[0] - 2), 2.0) + std::pow((x[1] - 1), 2);
}

std::vector<double> Function::evaluateFirstDerivative(std::vector<double>& x) {
    double d_dx1 = x[0] - 2;
    double d_dx2 = 2 * x[1] - 2;

    std::vector<double> v;
    v.push_back(d_dx1);
    v.push_back(d_dx2);

    return v;
}