#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <valarray>

using namespace std;

template<class element_type>
class matrix {
public:
    matrix(size_t width, size_t height) : m_stride(width), m_height(height), m_storage(width * height) {}

    element_type &operator()(size_t row, size_t column) {
        // column major
        return m_storage[slice(column, m_height, m_stride)][row];

        // row major
        return m_storage[slice(row, m_stride, m_height)][column];
    }

private:
    std::valarray<element_type> m_storage;
    size_t m_stride;
    size_t m_height;
};

int main() {
    matrix<double> m(20, 20);
    m(2, 3) = 5;

    valarray<double> vec(10);
    for (size_t i = 0; i < vec.size(); i++)
        vec[i] = 0.3 * i;
    valarray<double> sinc = sin(vec);
    vec.resize(2);

    for (size_t i = 0; i < vec.size(); ++i)
        std::cout << vec[i] << "sin:" << sinc[i] << "sin x" << sin(vec[i]) << endl;

    return 0;
}