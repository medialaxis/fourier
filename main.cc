#include <vector>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>

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

static Signal fft(Signal const& signal)
{
    return signal;
}

static Signal ifft(Signal const& spectrum)
{
    return spectrum;
}

static double error(Signal const& a, Signal const& b)
{
    auto e = a - b;
    return sqrt(std::real(dot(e, e)));
}

static double prop_inverse_dft(Signal const& test_signal)
{
    return error(test_signal, idft(dft(test_signal)));
}

static double prop_inverse_fft(Signal const& test_signal)
{
    return error(test_signal, ifft(fft(test_signal)));
}

static double prop_dft_equal_fft(Signal const& test_signal)
{
    return error(dft(test_signal), fft(test_signal));
}

#define TEST(signal) test(#signal, signal)

static void test(const char* test_name, double residue)
{
    if (residue < eps) {
        std::cout << test_name << ": PASS: residue=" << residue << std::endl;
    }
    else {
        std::cout << test_name << ": FAIL: residue=" << residue << std::endl;
    }
}

static Signal random_signal(size_t size)
{
    std::default_random_engine generator(0);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    Signal result(size);
    for (size_t i = 0; i != result.size(); ++i) {
        result[i] = distribution(generator);
    }
    return result;
}

int main()
{
    TEST(prop_inverse_dft(Signal(1024, 1)));
    TEST(prop_inverse_dft(random_signal(1024)));
    TEST(prop_inverse_fft(Signal(1024, 1)));
    TEST(prop_inverse_fft(random_signal(1024)));
    TEST(prop_dft_equal_fft(Signal(1024, 1)));
    TEST(prop_dft_equal_fft(random_signal(1024)));

    return 0;
}
