#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
using namespace std;
class Matrix
{
public:
    //construct

    Matrix();
    Matrix(const Matrix & m);
    Matrix(vector<double>,int , int );
    ~Matrix();

    //member value

    pair<int,int> shape;
    int size();

    //function
    Matrix reshape(int,int);
    Matrix T();

    //operate
    Matrix operator -();
    Matrix operator +();
    Matrix operator -(const Matrix & right_value);
    Matrix operator +(const Matrix & right_value);
    Matrix operator -(double right_value);
    Matrix operator +(double right_value);
    Matrix operator *(double right_value);
    Matrix operator /(double right_value);
    Matrix operator *(const Matrix &right_value);
    Matrix operator /(const Matrix &right_value);

    vector<double> operator[](int idx);

    //matrx compute
    double cross(const Matrix &right_value);
    double dot(const Matrix &right_value);

    //
    double sum();
    double max();
    double min();
    Matrix normal();
    double sum(int axis);
    double max(int axis);
    double min(int axis);
private:
    void descontruc();
    vector<vector<double> > vec;


};

//col vector,col num,row num
Matrix::Matrix(vector<double> array, int col, int row) : vec(col,vector<double>(row,1)) ,shape(col,row){
    assert(col != 0 && row != 0);
    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            vec[i][j] = array[i*row + j];
        }
    }
}

Matrix::Matrix(const Matrix & m) {
    this->shape = m.shape;
    this->vec = m.vec;
}


void Matrix::descontruc() {
    cout<<"destruction!!!";
    vector<vector<double> > ().swap(vec);
}

Matrix::~Matrix() {
    descontruc();

}


int Matrix::size() {
    return shape.first*shape.second;
}
Matrix Matrix::reshape(int col, int row) {
    assert(row*col == this->size());
    Matrix res(*this);

    vector<vector<double> > new_vec(col,vector<double>(row,1));

    for (int i = 0; i < col; i++) {
        for (int j = 0; j < row; j++) {
            new_vec[i][j] = vec[(i*row + j)/shape.second][(i*row +j)%shape.second];
//            new_vec[i][j] = vec[(i*col + j)/shape.second][(i*col +j)%shape.second];
        }
    }
    res.vec = new_vec;
    res.shape.first = col;
    res.shape.second = row;
    return res;

}

Matrix Matrix::T() {
    Matrix res(*this);

    vector<vector<double> > new_vec(shape.second,vector<double>(shape.first,1));

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

//-
Matrix  Matrix::operator-() {


    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = -1*res.vec[i][j];
        }
    }
    return res;

}

Matrix  Matrix::operator+() {


    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = 1*res.vec[i][j];
        }
    }
    return res;

}

Matrix Matrix::operator+(double right_value) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] + right_value;
        }
    }
    return res;
}

Matrix Matrix::operator-(double right_value) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] - right_value;
        }
    }
    return res;
}

Matrix Matrix::operator+(const Matrix & right_value) {
    assert(right_value.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] + right_value.vec[i][j] ;
        }
    }
    return res;
}

Matrix Matrix::operator-(const Matrix &right_value) {
    assert(right_value.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] - right_value.vec[i][j] ;

        }
    }
    return res;
}

Matrix Matrix::operator *(const Matrix &right_value)
{
    assert(right_value.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j] * right_value.vec[i][j] ;
        }
    }
    return res;
}

Matrix Matrix::operator /(const Matrix &right_value)
{
    assert(right_value.shape == shape);
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            assert(right_value.vec[i][j] != 0);
            res.vec[i][j] = res.vec[i][j]/ right_value.vec[i][j] ;
        }
    }
    return res;
}



Matrix Matrix::operator *(double right_value)
{
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            res.vec[i][j] = res.vec[i][j]* right_value ;
        }
    }
    return res;
}

Matrix Matrix::operator/(double right_value) {
    Matrix res(*this);
    for (int i = 0; i < res.shape.first; i++) {
        for (int j = 0; j < res.shape.second; j++) {
            assert(right_value != 0);

            res.vec[i][j] = res.vec[i][j] /right_value ;

        }
    }
    return res;
}

//norma
Matrix Matrix::normal() {
    double vec_sum = sqrt(((*this)*(*this)).sum());
    Matrix res = *this/vec_sum;
    return res;

}

//sum
double Matrix::sum() {
    double res=0;
    for (int i = 0; i < shape.first; i++) {
        for (int j = 0; j < shape.second; j++) {
            res += vec[i][j];
        }
    }
    return res;
}

double Matrix::min() {
    double res=0;
    for (int i = 0; i < shape.first; i++) {
        for (int j = 0; j < shape.second; j++) {
            if(res >vec[i][j])
                res=vec[i][j];
        }
    }
    return res;
}
double Matrix::max() {
    double res=0;
    for (int i = 0; i < shape.first; i++) {
        for (int j = 0; j < shape.second; j++) {
            if(res <vec[i][j])
                res=vec[i][j];
        }
    }
    return res;
}

//[]
vector<double> Matrix::operator[](int idx){
    return vec[idx];
}
//cross sum = x1*y2 - y1*x2
double Matrix::cross(const Matrix &right_value)
{
    assert(right_value.shape.first == 1&&shape.first == 1 &&right_value.shape.second == 2 && shape.second == 2);
    double cross_sum = vec[0][0]*right_value.vec[0][1] -  vec[0][1]*right_value.vec[0][0];
    return cross_sum;


}
//dot sum =x1*x2 + y1*y2 + z1*z2
double Matrix::dot(const Matrix &right_value)
{
    assert(right_value.shape.first == 1 && shape.first == 1);
    double dot_sum = (*this*right_value).sum();
    return dot_sum;

}



vector<double> array2vector(double *p,int size )
{
    return vector<double>(p,p+size);
}



int main() {

    double v[2] = {1,2,};
    vector<double> vv = array2vector(v,2);
    Matrix m(vv,2,1);
    m.reshape(2,1);
    Matrix mt(m.T());
    double dot_sum = mt.normal().cross(mt);
    double sum= (m*m).sum();


//    m.T();


    Matrix mm(-m);


    Matrix res = mm + mm;







    return 0;
}