#include<iostream>
#include<cmath> 

double f(double x1, double x2, double x3) {
    return std::pow(x1, 3.0) + std::pow(x1, 2.0) * x2 + x1 * std::pow(x2, 2.0) + std::pow(x2, 3.0) +
           x1 * x2 * x3 + std::pow(x3, 3.0) + x2 * std::pow(x3, 2.0);
}

int main(int argc, const char* argv[]) {
    std::cout << f(1,1,1) << std::endl;
    return 0;
}
