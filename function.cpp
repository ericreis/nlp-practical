#include "function.h"

Function::Function(double x1, double x2, double x3) {
    this->x1 = x1;
    this->x2 = x2;
    this->x3 = x3;
}

Function::~Function() {
    
}

double Function::evaluate() {
    return std::pow(this->x1, 3.0) + std::pow(this->x1, 2.0) * this->x2 + 
           this->x1 * std::pow(this->x2, 2.0) + std::pow(this->x2, 3.0) +
           this->x1 * this->x2 * this->x3 + std::pow(this->x3, 3.0) + 
           this->x2 * std::pow(x3, 2.0);
}