//
// Created by waxz on 2018/1/28.
//

#ifndef LINEARMATH_MATRIX_H
#define LINEARMATH_MATRIX_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <valarray>
using namespace std;

template<class Typ>
class Matrix {
public:
    //construct

    Matrix(const Matrix &r_v);

    Matrix(vector<Typ>, int, int);

    explicit Matrix(vector<vector<Typ> >);

    ~Matrix();

    //member value

    pair<int, int> shape;

    int size();

    vector<vector<Typ> > vec;

    //function
    Matrix reshape(int, int);

    Matrix T();

    //static function
    static Matrix ones(int, int);

    static Matrix zeros(int, int);

    static Matrix ones_like(Matrix &r_v);

    static Matrix zeros_like(Matrix &r_v);

    static Matrix arange(Typ, Typ, Typ);

    static Matrix linspace(Typ, Typ, int);

    static Matrix logspace(Typ, Typ, int);

    static Matrix concat(Matrix, Matrix);

    static Matrix stack(Matrix, Matrix);

    static Matrix where(Matrix);


    //operate
    Matrix operator-();

    Matrix operator+();

    Matrix operator-(const Matrix &r_v);

    Matrix operator+(const Matrix &r_v);

    Matrix operator-(Typ r_v);

    Matrix operator+(Typ r_v);

    Matrix operator*(Typ r_v);

    Matrix operator/(Typ r_v);

    Matrix operator*(const Matrix &r_v);

    Matrix operator/(const Matrix &r_v);

    Matrix operator>(double &r_v);

    Matrix operator<(double &r_v);


    vector<Typ> operator[](int idx);


    //matrx compute
    Typ cross(const Matrix &r_v);

    Typ dot(const Matrix &r_v);

    //
    Typ sum();

    Typ max();

    Typ min();

    Matrix normal();

    vector<Typ> sum(int axis);

    vector<Typ> max(int axis);

    vector<Typ> min(int axis);

private:
    void descontruc();


};

//col vector,col num,row num
template<class Typ>
Matrix<Typ>::Matrix(vector<Typ> array, int col, int row) : vec(col, vector<Typ>(row, (Typ) 1)), shape(col, row) {
    assert(col != 0 && row != 0);
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            vec[i][j] = array[i * row + j];
        }
    }
}

template<class Typ>
Matrix<Typ>::Matrix(vector<vector<Typ> > init_vec) {
    vec = init_vec;
    shape.first = init_vec.size();
    shape.second = init_vec[0].size();
}

template<class Typ>
Matrix<Typ>::Matrix(const Matrix &m) {
    this->shape = m.shape;
    this->vec = m.vec;
}

template<class Typ>
void Matrix<Typ>::descontruc() {
    cout << "destruction!!!";
    vector<vector<Typ> >().swap(vec);
}

template<class Typ>
Matrix<Typ>::~Matrix() {
    descontruc();

}

template<class Typ>
int Matrix<Typ>::size() {
    return shape.first * shape.second;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::reshape(int col, int row) {
    assert(row * col == this->size());
    Matrix res(*this);

    vector<vector<double> > new_vec(col, vector<double>(row, 1));

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            new_vec[i][j] = vec[(i * row + j) / shape.second][(i * row + j) % shape.second];
//            new_vec[i][j] = vec[(i*col + j)/shape.second][(i*col +j)%shape.second];
        }
    }
    res.vec = new_vec;
    res.shape.first = col;
    res.shape.second = row;
    return res;

}

template<class Typ>
Matrix<Typ> Matrix<Typ>::T() {
    Matrix res(*this);

    vector<vector<double> > new_vec(shape.second, vector<double>(shape.first, 1));

    double element;
    for (int i = 0; i < shape.second; i++) {
        for (int j = 0; j < shape.first; j++) {
            new_vec[i][j] = vec[j][i];
        }
    }
    res.vec = new_vec;
    res.shape.first = shape.second;
    res.shape.second = shape.first;
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator-() {


    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = -1 * res.vec[i][j];
        }
    }
    return res;

}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator+() {


    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = 1 * res.vec[i][j];
        }
    }
    return res;

}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator+(Typ r_v) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] + r_v;
        }
    }
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator-(Typ r_v) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] - r_v;
        }
    }
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator+(const Matrix &r_v) {
    assert(r_v.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] + r_v.vec[i][j];
        }
    }
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator-(const Matrix &r_v) {
    assert(r_v.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] - r_v.vec[i][j];

        }
    }
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator*(const Matrix &r_v) {
    assert(r_v.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] * r_v.vec[i][j];
        }
    }
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator/(const Matrix &r_v) {
    assert(r_v.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            assert(r_v.vec[i][j] != 0);
            res.vec[i][j] = res.vec[i][j] / r_v.vec[i][j];
        }
    }
    return res;
}


template<class Typ>
Matrix<Typ> Matrix<Typ>::operator*(Typ r_v) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] * r_v;
        }
    }
    return res;
}

template<class Typ>
Matrix<Typ> Matrix<Typ>::operator/(Typ r_v) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            assert(r_v != 0);

            res.vec[i][j] = res.vec[i][j] / r_v;

        }
    }
    return res;
}

//norma
template<class Typ>
Matrix<Typ> Matrix<Typ>::normal() {
    Typ vec_sum = sqrt(((*this) * (*this)).sum());
    Matrix res = *this / vec_sum;
    return res;

}

//sum
template<class Typ>
Typ Matrix<Typ>::sum() {
    Typ res = 0;
    for (int i = 0; i < shape.first; i++) {
        for (int j = 0; j < shape.second; j++) {
            res += vec[i][j];
        }
    }
    return res;
}

template<class Typ>
Typ Matrix<Typ>::min() {
    Typ res = 0;
    for (int i = 0; i < shape.first; i++) {
        for (int j = 0; j < shape.second; j++) {
            if (res > vec[i][j])
                res = vec[i][j];
        }
    }
    return res;
}

template<class Typ>
Typ Matrix<Typ>::max() {
    Typ res = 0;
    for (int i = 0; i < shape.first; i++) {
        for (int j = 0; j < shape.second; j++) {
            if (res < vec[i][j])
                res = vec[i][j];
        }
    }
    return res;
}

//[]
template<class Typ>
vector<Typ> Matrix<Typ>::operator[](int idx) {
    return vec[idx];
}

//cross sum = x1*y2 - y1*x2
template<class Typ>
Typ Matrix<Typ>::cross(const Matrix &r_v) {
    assert(r_v.shape.first == 1 && shape.first == 1 && r_v.shape.second == 2 && shape.second == 2);
    Typ cross_sum = vec[0][0] * r_v.vec[0][1] - vec[0][1] * r_v.vec[0][0];
    return cross_sum;


}

//dot sum =x1*x2 + y1*y2 + z1*z2
template<class Typ>
Typ Matrix<Typ>::dot(const Matrix &r_v) {
    assert(r_v.shape.first == 1 && shape.first == 1);
    Typ dot_sum = (*this * r_v).sum();
    return dot_sum;

}


//static function
template<class Typ>
Matrix<Typ> Matrix<Typ>::ones(int col, int row) {
    vector<vector<Typ> > new_vec(col, vector<Typ>(row, (Typ) 1));
    Matrix res(new_vec);
    return res;


}

template<class Typ>
Matrix<Typ> Matrix<Typ>::zeros(int col, int row) {
    vector<vector<Typ> > new_vec(col, vector<Typ>(row, (Typ) 0));
    Matrix res(new_vec);
    return res;


}

template<class Typ>
Matrix<Typ> Matrix<Typ>::ones_like(Matrix &r_v) {
    cout << "ones_like";
    vector<vector<Typ> > new_vec(r_v.shape.first, vector<Typ>(r_v.shape.second, (Typ) 1));
    Matrix res(new_vec);
    return res;

}

template<class Typ>
Matrix<Typ> Matrix<Typ>::zeros_like(Matrix &r_v) {
    cout << "ones_like";
    vector<vector<Typ> > new_vec(r_v.shape.first, vector<Typ>(r_v.shape.second, (Typ) 0));
    Matrix res(new_vec);
    return res;

}

#endif //LINEARMATH_MATRIX_H
