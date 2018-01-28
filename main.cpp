#include "Matrix.h"

vector<double> array2vector(double *p,int size )
{
    return vector<double>(p,p+size);
}



int main() {

    double v[2] = {1,2,};
    vector<double> vv = array2vector(v,2);
    Matrix<double> m(vv, 2, 1);
    m.reshape(2,1);
    Matrix<double> mt(m.T());
    double dot_sum = mt.normal().cross(mt);
    double sum= (m*m).sum();


//    m.T();


    Matrix<double> mm(-m);


    Matrix<double> res = mm + mm;


    Matrix<bool> o = Matrix<bool>::ones(2, 3);
    int num = 4;
    bool b = (bool) num;
    Matrix<double> mo = Matrix<double>::ones_like(m);
    cout << endl << 233 << "bool" << b << endl;

    return 0;
}