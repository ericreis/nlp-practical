#include <cmath>
#include <vector>

class Function {
public:
    Function();
    ~Function();
    double evaluate(std::vector<double>& x, double p);
    std::vector<double> evaluateFirstDerivative(std::vector<double>& x, 
                                                double p);
};
