#include "bfgs.h"

// difference of two vectors
template<class Type>
std::vector<Type> operator-(std::vector<Type> x, const std::vector<Type>& y) {
    for (int i = 0, sz = x.size(); i < sz; i++) {
        x[i] -= y[i];
    }
    return x;
}

std::vector< std::vector<double> > BFGS::method(std::vector< std::vector<double> > hk, 
                                                std::vector<double> x, 
                                                std::vector<double> x_,
                                                std::vector<double> gradf, 
                                                std::vector<double> gradf_) {
    std::vector<double> pk = x - x_;
    std::vector<double> qk = gradf - gradf_;
    std::vector< std::vector<double> > newHk;

    double denominator = this->vectorTXvector(pk, qk);

    double firstComponent = this->vectorTXvector(this->vectorTXmatrix(qk, hk), qk) / denominator;
    
    std::vector< std::vector<double> > secondComponent = this->vectorXvectorT(pk, pk);
    for (int i = 0; i < secondComponent.size(); i++) {
        for (int j = 0; j < secondComponent[0].size(); j++) {
            secondComponent[i][j] = (1 + firstComponent) * (secondComponent[i][j] / denominator);
        }
    }

    std::vector< std::vector<double> > thirdComponent1 = this->matrixXmatrix(this->vectorXvectorT(pk, qk), hk);
    std::vector< std::vector<double> > thirdComponent2 = this->vectorXvectorT(this->matrixXvector(hk, qk), pk);
    std::vector< std::vector<double> > thirdComponent;
    for (int i = 0; i < thirdComponent1.size(); i++) {
        std::vector<double> row;
        for (int j = 0; j < thirdComponent1[0].size(); j++) {
            row.push_back((thirdComponent1[i][j] + thirdComponent2[i][j]) / denominator);
        }
        thirdComponent.push_back(row);
    }

    for (int i = 0; i < hk.size(); i++) {
        std::vector<double> row;
        for (int j = 0; j < hk[0].size(); j++) {
            row.push_back(hk[i][j] + secondComponent[i][j] - thirdComponent[i][j]);
        }
        newHk.push_back(row);
    }

    return newHk;
}

double BFGS::vectorTXvector(std::vector<double> vector1, 
                            std::vector<double> vector2) {
    double result = 0.0;
    for (int i = 0; i < vector1.size(); i++) {
        result += vector1[i] * vector2[i];
    }
    return result;
}

std::vector< std::vector<double> > BFGS::vectorXvectorT(std::vector<double> vector1,
                                                        std::vector<double> vector2) {
    std::vector< std::vector<double> > result;
    for (int i = 0; i < vector1.size(); i++)
    {
        std::vector<double> row;
        for (int j = 0; j < vector2.size(); j++)
        {
            row.push_back(vector1[i] * vector2[j]);
        }
        result.push_back(row);
    }
    return result;
}

std::vector<double> BFGS::vectorTXmatrix(std::vector<double> vector1, 
                                         std::vector< std::vector<double> > matrix) {
    std::vector<double> result;
    for (int j = 0; j < vector1.size(); j++) {
        double value = 0.0;
        for (int i = 0; i < matrix[0].size(); i++) {
            value += matrix[i][j] * vector1[i];
        }
        result.push_back(value);
    }
    return result;
}

std::vector< std::vector<double> > BFGS::matrixXmatrix (std::vector< std::vector<double> > matrix1,
                                                        std::vector< std::vector<double> > matrix2) {
    std::vector< std::vector<double> > result;
    for (int i = 0; i < matrix1.size(); i++) {
        std::vector<double> row;
        for (int j = 0; j < matrix2[0].size(); j++) {
            double value = 0.0;
            for (int k = 0; k < matrix1[0].size(); k++) {
                value += matrix1[i][k] * matrix2[k][j];
            }
            row.push_back(value);
        }
        result.push_back(row);
    }
    return result;
}

std::vector<double> BFGS::matrixXvector (std::vector< std::vector<double> > matrix, 
                                         std::vector<double> vector1) {
    std::vector<double> result;
    for (int i = 0; i < matrix.size(); i++) {
        double value = 0.0;
        for (int j = 0; j < matrix[0].size(); j++) {
            value += matrix[i][j] * vector1[j];
        }
        result.push_back(value);
    }
    return result;
}