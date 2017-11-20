#include <cmath>

class Function {
private:
    double x1, x2, x3;
public:
    Function(double x1, double x2, double x3);
    ~Function();
    double evaluate();
};
