#include <iostream>

#include "function.h"

// sum of two vectors
template<class Type>
std::vector<Type> operator+(std::vector<Type> x, const std::vector<Type> &y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] += y[i];
    }
    return x;
}

// product between a vector and a scalar
template<class Type>
std::vector<Type> operator*(std::vector<Type> x, const Type &y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] *= y;
    }
    return x;
}

// product between a vector and a vector
template<class Type>
Type operator*(std::vector<Type> &x, const std::vector<Type> &y) {
    double inner_product = 0.0;
    for (int i = 0, sz = x.size(); i < sz; i++) {
        inner_product += x[i] * y[i];
    }
    return inner_product;
}

double armijo(std::vector<double> x_, std::vector<double> d, 
              double gamma, double eta, Function& f) {
    if (gamma <= 0 || gamma >= 1) {
        std::cout << "Gamma must be in (0,1)" << std::endl;
    }
    if (eta <= 0 || eta >= 1) {
        std::cout << "Eta must be in (0,1)" << std::endl;
    }

    double t = 1.0;
    std::vector<double> v = x_ + (d * t);
    std::vector<double> gradf = f.evaluateFirstDerivative(x_);
    double gradf_d = gradf * d;
    int it_count = 0;
    while (f.evaluate(v) > f.evaluate(x_) + eta * t * gradf_d) {
        t *= gamma;
        v = x_ + (d * t);
        it_count++;
    }

    std::cout << "ARMIJO iterations = " << it_count << std::endl;

    return t;
}

int main(int argc, const char* argv[]) {
    Function f;
    std::vector<double> x;
    x.push_back(1);
    x.push_back(1);
    x.push_back(2);
    std::vector<double> d;
    d.push_back(3);
    d.push_back(1);
    // d.push_back(1);
    // std::cout << armijo(x, d, 0.8, 0.25, f) << std::endl;
    // std::cout << f.evaluate(x) << std::endl;
    std::vector<double> gradf = f.evaluateFirstDerivative(x);
    for (int i = 0; i < gradf.size(); i++) {
        std::cout << gradf[i] << std::endl;
    }
    return 0;
}


