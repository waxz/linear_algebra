//
// Created by waxz on 2018/2/1.
//

#ifndef LINEARMATH_TM_H
#define LINEARMATH_TM_H

class M {
public:
    template<class T, class S>
    double Combiner(T func, S func1, double a, double b);

};


template<class T, class S>
double M::Combiner(T func, S func1, double a, double b) {
    func1(a);
    return func(a, b);
}


double Add(double a, double b) {
    return a + b;
}

double Mult(double a, double b) {
    return a * b;
}

#endif //LINEARMATH_TM_H
