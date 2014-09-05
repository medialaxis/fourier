#include <vector>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <complex>
#include <iostream>

typedef std::complex<double> Complex;
typedef std::vector<Complex> Signal;

static const double eps = 0.01;
static const Complex i(0, 1);

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
    Signal result(signal.size());

    for (size_t k = 0; k != result.size(); ++k) {
        Complex c(0, 0);

        for (size_t n = 0; n != signal.size(); ++n) {
            c += signal[n]*exp(-i*2.0*M_PI*(double) k*(double) n/(double) signal.size());
        }

        result[k] = c;
    }

    return result;
}

static Signal idft(Signal const& spectrum)
{
    Signal result(spectrum.size());

    for (size_t n = 0; n != result.size(); ++n) {
        Complex c(0, 0);

        for (size_t k = 0; k != spectrum.size(); ++k) {
            c += spectrum[k]*exp(i*2.0*M_PI*(double) k*(double) n/(double) spectrum.size());
        }

        result[n] = c/(double) spectrum.size();
    }

    return result;
}

static double error(Signal const& a, Signal const& b)
{
    auto e = a - b;
    return sqrt(std::real(dot(e, e)));
}

int main()
{
    Signal signal(1024, 1);
    Signal spectrum = dft(signal);
    Signal signal2 = idft(spectrum);
    double e = error(signal, signal2);
    std::cout << "e=" << e << "\n";
    if (e < eps) {
        std::cout << "PASS\n";
    }
    else {
        std::cout << "FAIL\n";
    }

    return 0;
}
