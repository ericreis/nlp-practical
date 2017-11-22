#include <iostream>

#include "function.h"

// sum of two vectors
template<class Type>
std::vector<Type> operator+(std::vector<Type> x, const std::vector<Type>& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] += y[i];
    }
    return x;
}

// sum of two vectors (+=)
template<class Type>
void operator+=(std::vector<Type>& x, const std::vector<Type>& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] += y[i];
    }
}

// product between a vector and a scalar
template<class Type>
std::vector<Type> operator*(std::vector<Type> x, const Type& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] *= y;
    }
    return x;
}

// product between a vector and a vector
template<class Type>
Type operator*(std::vector<Type>& x, const std::vector<Type>& y) {
    double inner_product = 0.0;
    for (int i = 0, sz = x.size(); i < sz; i++) {
        inner_product += x[i] * y[i];
    }
    return inner_product;
}

// inequality between vector and scalar
template<class Type>
bool operator!=(const std::vector<double>& x, const Type& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        if (x[i] != y) return 1;
    }
    return 0;
}

// negative of vector
template<class Type>
std::vector<Type> operator-(const std::vector<Type>& y) {
    std::vector<Type> x(y.size());
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] = -y[i];
    }
    return x;
}

// negative of vector
template<class Type>
std::ostream& operator<<(std::ostream& os, const std::vector<Type>& v) {
    os << "[";
    for (std::size_t i = 0; i < v.size(); i++) {
        os << v[i] << (i == v.size() - 1 ? "" : ",");
    }
    return os << "]";
}

double armijo(std::vector<double> x_, std::vector<double>& d, double gamma,
              double eta, double p, Function &f) {
    if (gamma <= 0 || gamma >= 1) {
        std::cout << "Gamma must be in (0,1)" << std::endl;
    }
    if (eta <= 0 || eta >= 1) {
        std::cout << "Eta must be in (0,1)" << std::endl;
    }

    double t = 1.0;
    std::vector<double> v = x_ + (d * t);
    std::vector<double> gradf = f.evaluateFirstDerivative(x_, p);
    double gradf_d = gradf * d;
    int it_count = 0;
    while (f.evaluate(v, p) > f.evaluate(x_, p) + eta * t * gradf_d) {
        t *= gamma;
        v = x_ + (d * t);
        it_count++;
    }

    // std::cout << "ARMIJO iterations = " << it_count << std::endl;

    return t;
}

bool gradient_has_converged(std::vector<double>& gradf, double tol=1e-5) {
    for (int i = 0, sz = gradf.size(); i < sz; i++) {
        // std::cout << (gradf[i] > tol || gradf[i] < -tol) << std::endl;
        if (gradf[i] > tol || gradf[i] < -tol) return 0;
    }
    return 1;
}

bool vector_has_converged(std::vector<double>& v1, std::vector<double>& v2,
                          double tol=1e-5) {
    int v1_sz = v1.size();
    int v2_sz = v2.size();
    if (v1_sz != v2_sz) {
        return 0;
    }
    for (int i = 0; i < v1_sz; i++) {
        if (v1[i] - v2[i] > tol || v1[i] - v2[i] < -tol) return 0;
    }
    return 1;
}

bool solution_has_converged(double sol1, double sol2, double tol=1e-5) {
    return !(sol1 - sol2 > tol || sol1 - sol2 < -tol);
}

std::vector<double> gradient(std::vector<double> x, double p, Function& f) {
    int it_count = 0;
    std::vector<double> gradf = f.evaluateFirstDerivative(x, p);
    while (!gradient_has_converged(gradf)) {
        std::vector<double> d = -f.evaluateFirstDerivative(x, p);
        double t = armijo(x, d, 0.8, 0.25, p, f);
        x += d * t;
        gradf = f.evaluateFirstDerivative(x, p);
        it_count++;
    }
    std::cout << "GRADIENT iterarions = " << it_count << std::endl;
    return x;
}

std::vector<double> external_penalty(std::vector<double> x, double& p, Function& f,
                                     std::vector<double> (*method)(std::vector<double> x, 
                                                                   double p, 
                                                                   Function& f),
                                     int beta=3) {
    int it_count = 0;
    std::vector<double> xk = (*method)(x, p, f);
    while (!solution_has_converged(f.evaluate(x, p), f.evaluate(xk, p), 1e-5)) {
        p *= beta;
        x = xk;
        xk = (*method)(x, p, f);
        it_count++;
    }

    std::cout << "EXTERNAL PENALTY iterations = " << it_count << std::endl;

    return x;
}

int main(int argc, const char* argv[]) {
    Function f;
    std::vector<double> x;
    x.push_back(1);
    x.push_back(1);
    x.push_back(1);
    double p = 1;
    std::vector<double> x_opt = external_penalty(x, p, f, gradient);
    std::cout << "x_opt = [";
    for (int i = 0; i< x_opt.size(); i++) {
        std::cout << x_opt[i] << " ";
    }
    std::cout << "]" << std::endl;
    std::cout << "f(x_opt) = " << f.evaluate(x_opt, p) << std::endl;
    return 0;
}


