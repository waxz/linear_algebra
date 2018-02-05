#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <valarray>
#include "tm.h"
#include <numeric>


using namespace std;

void cbk(double a) {
    cout << "cbk" << a + 12 << endl;
}

int main() {
    //array convert to vector
    //vector convert to arrar
    //array convert to valarray
    double arr1[] = {1, 2, 3, 4};
    vector<double> vec(arr1, arr1 + 4);
    double *arr2 = &vec[0];
    valarray<double> valarr(arr2, 4);

    for (int i = 0; i < valarr.size(); i++) {
        cout << valarr[i];
    }

    M m;

    m.Combiner(Add, cbk, 12, 13);
    m.Combiner(Mult, cbk, 12, 13);

    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    double mean = sum / vec.size();

    double sq_sum = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / vec.size() - mean * mean);

    return 0;
}