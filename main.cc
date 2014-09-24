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

static size_t reverse_bits(size_t n, size_t max)
{
    size_t result = 0;

    for (size_t i = 1; i != max; i <<= 1) {
        result <<= 1;
        result |= (n & 1);
        n >>= 1;
    }

    return result;
}

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

static inline Complex W(int k, int N)
{
    return exp(-i*2.0*M_PI*(double) k/(double) N);
}

static inline Complex Q(int n, int N)
{
    return exp(i*2.0*M_PI*(double) n/(double) N);
}

static void fft_step(Complex* spectrum, size_t spectrumSize)
{
    for (size_t i = 0; i != spectrumSize/2; ++i) {
        size_t sample1 = i;
        size_t sample2 = i + spectrumSize/2;

        Complex even = spectrum[sample1];
        Complex odd = spectrum[sample2];

        spectrum[sample1] = even + W(sample1, spectrumSize)*odd;
        spectrum[sample2] = even + W(sample2, spectrumSize)*odd;
    }
}

static void print_signal(std::string name, Signal const& signal) __attribute__((unused));
static void print_signal(std::string name, Signal const& signal)
{
    std::cout << name << "=";
    for (size_t i = 0; i != signal.size(); ++i) {
        std::cout << signal[i] << ",";
    }
    std::cout << "\n";
}

static Signal fft(Signal const& signal)
{
    size_t const N = signal.size();
    Signal result(N);

    for (size_t i = 0; i != N; ++i) {
        result[i] = signal[reverse_bits(i, N)];
    }

    size_t transform_count = N/2;
    while (transform_count >= 1) {
        size_t sample_count = N/transform_count;

        for (size_t transform = 0; transform != transform_count; ++transform) {
            fft_step(&result[transform*sample_count], sample_count);
        }

        transform_count >>= 1;
    }

    return result;
}

static void ifft_step(Complex* spectrum, size_t spectrumSize)
{
    for (size_t i = 0; i != spectrumSize/2; ++i) {
        size_t sample1 = i;
        size_t sample2 = i + spectrumSize/2;

        Complex even = spectrum[sample1];
        Complex odd = spectrum[sample2];

        spectrum[sample1] = 0.5*(even + Q(sample1, spectrumSize)*odd);
        spectrum[sample2] = 0.5*(even + Q(sample2, spectrumSize)*odd);
    }
}

static Signal ifft(Signal const& spectrum)
{
    size_t const N = spectrum.size();
    Signal result(N);

    for (size_t i = 0; i != N; ++i) {
        result[i] = spectrum[reverse_bits(i, N)];
    }

    size_t sample_count = 2;
    size_t transform_count = N/sample_count;
    while (sample_count <= N) {
        assert(transform_count*sample_count == N);

        for (size_t transform = 0; transform != transform_count; ++transform) {
            ifft_step(&result[transform*sample_count], sample_count);
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

static double prop_idft_equal_ifft(Signal const& test_signal)
{
    return error(idft(test_signal), ifft(test_signal));
}

static double prop_fft_is_decomposed_dft(Signal const& test_signal)
{
    Signal even_samples;
    Signal odd_samples;

    // Partition even and odd samples.
    for (size_t i = 0; i != test_signal.size()/2; ++i) {
        even_samples.push_back(test_signal[2*i]);
        odd_samples.push_back(test_signal[2*i + 1]);
    }

    Signal even_spectrum = dft(even_samples);
    Signal odd_spectrum = dft(odd_samples);

    Signal intermediate_spectrum;
    intermediate_spectrum.insert(intermediate_spectrum.end(), even_spectrum.begin(), even_spectrum.end());
    intermediate_spectrum.insert(intermediate_spectrum.end(), odd_spectrum.begin(), odd_spectrum.end());

    fft_step(&intermediate_spectrum[0], intermediate_spectrum.size());

    return error(dft(test_signal), intermediate_spectrum);
}

static bool prop_reverse_bits(size_t n, size_t max, size_t correct)
{
    return reverse_bits(n, max) == correct;
}

#define TEST_RESIDUE(signal) test_residue(#signal, signal)

static void test_residue(const char* test_name, double residue)
{
    if (residue < eps) {
        std::cout << test_name << ": PASS: residue=" << residue << std::endl;
    }
    else {
        std::cout << test_name << ": FAIL: residue=" << residue << std::endl;
    }
}

#define TEST(prop) test(#prop, prop)

static void test(const char* test_name, bool result)
{
    if (result) {
        std::cout << test_name << ": PASS" << std::endl;
    }
    else {
        std::cout << test_name << ": FAIL" << std::endl;
    }
}

static Signal random_signal(size_t size)
{
    std::default_random_engine generator(0);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    Signal result(size);
    for (size_t i = 0; i != result.size(); ++i) {
        result[i] = Complex(distribution(generator), distribution(generator));
    }
    return result;
}

int main()
{
    TEST_RESIDUE(prop_inverse_dft(Signal(1024, 1)));
    TEST_RESIDUE(prop_inverse_dft(random_signal(1024)));
    TEST_RESIDUE(prop_inverse_fft(Signal(2, 1)));
    TEST_RESIDUE(prop_inverse_fft(Signal(1024, 1)));
    TEST_RESIDUE(prop_inverse_fft(random_signal(1024)));
    TEST_RESIDUE(prop_dft_equal_fft(Signal(1024, 1)));
    TEST_RESIDUE(prop_dft_equal_fft(random_signal(4)));
    TEST_RESIDUE(prop_dft_equal_fft(random_signal(1024)));
    TEST_RESIDUE(prop_dft_equal_fft(Signal{7,6,5,4,3,2,i,0}));
    TEST_RESIDUE(prop_idft_equal_ifft(random_signal(1024)));
    TEST_RESIDUE(prop_idft_equal_ifft(Signal{1.1,i,2.1,3}));
    TEST_RESIDUE(prop_fft_is_decomposed_dft(Signal{7,6,5,4,3,2,i,0}));
    TEST_RESIDUE(prop_fft_is_decomposed_dft(random_signal(1024)));
    TEST(prop_reverse_bits(0xAA, 0x100, 0x55));
    TEST(prop_reverse_bits(0xA5, 0x100, 0xA5));

    return 0;
}
