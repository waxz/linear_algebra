#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <valarray>

using namespace std;

int main() {

    valarray<double> vec(10);
    for (size_t i = 0; i < vec.size(); i++)
        vec[i] = 0.3 * i;
    valarray<double> sinc = sin(vec);

    for (size_t i = 0; i < vec.size(); ++i)
        std::cout << vec[i] << "sin:" << sinc[i] << "sin x" << sin(vec[i]) << endl;

    return 0;
}