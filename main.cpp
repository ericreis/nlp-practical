#include <iostream>

#include "function.h"

int main(int argc, const char* argv[]) {
	Function f(1,1,1);
    std::cout << f.evaluate() << std::endl;
    return 0;
}


