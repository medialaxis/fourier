#include <vector>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <complex>
#include <iostream>

static const double eps = 0.01;

typedef std::complex<double> Complex;
typedef std::vector<Complex> Signal;

static Signal operator-(Signal const& a, Signal const& b)
{
    assert(a.size() == b.size());
    Signal result(a.size());

    for (size_t i = 0; i != result.size(); ++i) {
        result[i] = a[i] - b[i];
    }

    return result;
}

static Complex dot(Complex const& a, Complex const& b)
{
    return a*std::conj(b);
}

static Complex dot(Signal const& a, Signal const& b)
{
    assert(a.size() == b.size());

    Complex result = 0;
    for (size_t i = 0; i != a.size(); ++i) {
        result += dot(a[i], b[i]);
    }
    return result;
}

static Signal dft(Signal const& signal)
{
    return signal;
}

static Signal idft(Signal const& spectrum)
{
    return spectrum;
}

static double error(Signal const& a, Signal const& b)
{
    auto e = a - b;
    return sqrt(std::real(dot(e, e)));
}

int main()
{
    Signal signal(1024);
    Signal spectrum = dft(signal);
    Signal signal2 = idft(spectrum);
    if (error(signal, signal2) < eps) {
        std::cout << "PASS\n";
    }
    else {
        std::cout << "FAIL\n";
    }

    return 0;
}
