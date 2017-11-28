#include <iostream>
#include <sstream> 
#include <string>
#include <limits>

#include "function.h"
#include "bfgs.h"

typedef std::numeric_limits<double> dbl;

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

// difference of two vectors
template<class Type>
std::vector<Type> operator-(std::vector<Type> x, const std::vector<Type>& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] -= y[i];
    }
    return x;
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
Type operator*(std::vector<Type> x, const std::vector<Type> y) {
    double inner_product = 0.0;
    for (int i = 0, sz = x.size(); i < sz; i++) {
        inner_product += x[i] * y[i];
    }
    return inner_product;
}

// product between a vector and a matrix
template<class Type>
std::vector<Type> operator*(const std::vector<Type>& x, const std::vector< std::vector<Type> >& m) {
    std::vector<double> result;
    for (int j = 0; j < m[0].size(); j++) {
        double sum = 0.0;
        for (int i = 0; i < x.size(); i++) {
            sum += x[i] * m[i][j];
        }
        result.push_back(sum);
    }
    return result;
}

// product between a scalar and a matrix
template<class Type>
std::vector< std::vector<Type> > operator*(const Type& k, std::vector< std::vector<Type> > m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            m[i][j] *= k;
        }
    }
    return m;
}

// inequality between vector and scalar
template<class Type>
bool operator!=(const std::vector<double>& x, const Type& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        if (x[i] != y) return 1;
    }
    return 0;
}

// inequality between vector and scalar
template<class Type>
bool operator==(const std::vector<double>& x, const Type& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        if (x[i] != y) return 0;
    }
    return 1;
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

// print vector
template<class Type>
std::ostream& operator<<(std::ostream& os, const std::vector<Type>& v) {
    os << "[";
    for (std::size_t i = 0; i < v.size(); i++) {
        os << v[i] << (i == v.size() - 1 ? "" : ",");
    }
    return os << "]";
}

// print matrix
template<class Type>
std::ostream& operator<<(std::ostream& os, 
                         const std::vector< std::vector<Type> >& m) {
    for (std::size_t i = 0; i < m.size(); i++) {
        os << m[i] << (i == m.size() - 1 ? "" : "\n");
    }
    return os;
}

std::vector< std::vector<double> > identity(int size) {
    std::vector< std::vector<double> > identity;
    for (int i = 0; i < size; i++) {
        std::vector<double> r;
        for (int j = 0; j < size; j++) {
            i == j ? r.push_back(1.0) : r.push_back(0.0);
        }
        identity.push_back(r);
    }
    return identity;
}

double armijo(std::vector<double> x_, std::vector<double>& d, double gamma,
              double eta, double p, Function &f) {
    if (eta < 0 || eta > 1) {
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

std::vector<double> quase_newton(std::vector<double> x, double p, Function& f) {
    int it_count = 0;
    std::vector< std::vector<double> > hk = identity(x.size());
    std::vector<double> gradf = f.evaluateFirstDerivative(x, p);
    std::vector<double> x_;
    std::vector<double> gradf_;
    BFGS bfgs;
    while (!gradient_has_converged(gradf)) {
        // std::cout << "gradf = " << gradf << std::endl;
        // std::cout << "hk = " << hk << std::endl;
        std::vector<double> d = -(gradf * hk);
        double t = armijo(x, d, 0.8, 0.25, p, f);
        x_ = x;
        // std::cout << "d = " << d << std::endl;
        // std::cout << "t = " << t << std::endl;
        x += d * t;
        // std::cout << "d * t = " << d * t << std::endl;
        // std::cout << "x_ = " << x_ << std::endl;
        // std::cout << "x = " << x << std::endl;
        // std::cout << std::endl;
        gradf_ = gradf;
        gradf = f.evaluateFirstDerivative(x, p);
        if (x - x_ == 0 || gradf - gradf_ == 0) {
            std::cout << "BFGS pk or qk == 0" << std::endl;
            std::cout << "x = " << x <<std::endl;
            std::cout << "delta_x = " << x - x_ <<std::endl;
            std::cout << "delta_gradf = " << gradf - gradf_ <<std::endl;
            break;
        }
        hk = bfgs.method(hk, x, x_, gradf, gradf_);
        it_count++;
    }
    std::cout << "QUASE NEWTON iterarions = " << it_count << std::endl;
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

static void show_usage(std::string name) {
    std::cerr << "Usage: " << name << " <option(s)>\n\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-a,--algo\t\tAlgorithm to use (gradient, quase-newton)\n"
              << "\t-x0,--initialx\t\tInitial point x (R3) for the iterative algorithms\n"
              << "\t-p0,--initialp\t\tInitial weight for the external penalty algorithm\n"
              << std::endl;
}

int main(int argc, const char* argv[]) {
    // std::cout.precision(50);

    if (argc < 2) {
        show_usage(argv[0]);
        return 1;
    }
    std::vector<double> x0;
    double p0;
    std::string algo;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string delimiter = ",";
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else if ((arg == "-a") || (arg == "--algo")) {
            algo = argv[++i];
        } else if ((arg == "-x0") || (arg == "--initialx")) {
            std::istringstream iss(argv[++i]);
            while(std::getline(iss, arg, ',')) {
                x0.push_back(std::stod(arg));
            }
        } else if ((arg == "-p0") || (arg == "--initialp")) {
            p0 = std::stod(argv[++i]);
        } else {
            std::cerr << "Unkown option " << arg << std::endl;
            show_usage(argv[0]);
        }
    }
    Function f;
    std::vector<double> x_opt;
    if (algo == "gradient") {
        x_opt = external_penalty(x0, p0, f, gradient);
    } else if (algo == "quase-newton") {
        x_opt = external_penalty(x0, p0, f, quase_newton);
    }
    std::cout << "x_opt = " << x_opt << std::endl;
    std::cout << "f(x_opt) = " << f.evaluate(x_opt, p0) << std::endl;

    return 0;
}


