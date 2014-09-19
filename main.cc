#include <vector>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <boost/noncopyable.hpp>

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

static inline Complex w(int k, int N)
{
    return exp(-i*2.0*M_PI*(double) k/(double) N);
}

static inline Complex q(int n, int N)
{
    return exp(i*2.0*M_PI*(double) n/(double) N);
}

static Signal fft(Signal const& signal)
{
    size_t const N = signal.size();
    Signal result(N);

    for (size_t i = 0; i != N/2; ++i) {
        result[2*i] = signal[i];
        result[2*i+1] = signal[N/2+i];
    }

    size_t transform_count = N/2;
    size_t sample_count = 2;
    while (sample_count <= N) {
        assert(transform_count*sample_count == N);

        for (size_t transform = 0; transform != transform_count; ++transform) {
            for (size_t sample = 0; sample != sample_count/2; ++sample) {
                size_t offset = transform*sample_count;
                size_t sample1 = sample;
                size_t sample2 = sample + sample_count/2;

                Complex even = result[offset+sample1];
                Complex odd = result[offset+sample2];

                result[offset+sample1] = even + w(sample1, sample_count)*odd;
                result[offset+sample2] = even + w(sample2, sample_count)*odd;
            }
        }

        transform_count >>= 1;
        sample_count <<= 1;
    }

    return result;
}

static Signal ifft(Signal const& spectrum)
{
    size_t const N = spectrum.size();
    Signal result(N);

    for (size_t i = 0; i != N/2; ++i) {
        result[2*i] = spectrum[i];
        result[2*i+1] = spectrum[N/2+i];
    }

    size_t transform_count = N/2;
    size_t sample_count = 2;
    while (sample_count <= N) {
        assert(transform_count*sample_count == N);

        for (size_t transform = 0; transform != transform_count; ++transform) {
            for (size_t sample = 0; sample != sample_count/2; ++sample) {
                size_t offset = transform*sample_count;
                size_t sample1 = sample;
                size_t sample2 = sample + sample_count/2;

                Complex even = result[offset+sample1];
                Complex odd = result[offset+sample2];

                result[offset+sample1] = Complex(2,0)*(even + q(sample1, sample_count)*odd);
                result[offset+sample2] = Complex(2,0)*(even + q(sample2, sample_count)*odd);
            }
        }

        transform_count >>= 1;
        sample_count <<= 1;
    }

    return result;
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
    TEST(prop_dft_equal_fft(random_signal(4)));
    TEST(prop_dft_equal_fft(random_signal(1024)));
    TEST(prop_dft_equal_fft(Signal(2, 1)));
    TEST(prop_dft_equal_fft(Signal(4, 1)));
    TEST(prop_dft_equal_fft(Signal(8, 1)));

    return 0;
}
