#include <iostream>
#include <vector>

class BFGS
{
public:
    std::vector< std::vector<double> > method(std::vector< std::vector<double> > hk, 
                                              std::vector<double> x, 
                                              std::vector<double> x_,
                                              std::vector<double> gradf, 
                                              std::vector<double> gradf_);
    double vectorTXvector(std::vector<double> vector1, 
                          std::vector<double> vector2);
    std::vector< std::vector<double> > vectorXvectorT(std::vector<double> vector1,
                                                      std::vector<double> vector2);
    std::vector<double> vectorTXmatrix(std::vector<double> vector1, 
                                       std::vector< std::vector<double> > matrix);
    std::vector< std::vector<double> > matrixXmatrix (std::vector< std::vector<double> > matrix1,
                                                      std::vector< std::vector<double> > matrix2);
    std::vector<double> matrixXvector (std::vector< std::vector<double> > matrix, 
                                       std::vector<double> vector1);
};