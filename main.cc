#include <vector>
#include <cassert>
#include <cstddef>
#include <cmath>
#include <complex>
#include <map>
#include <iostream>
#include <fstream>
#include <random>
#include <boost/noncopyable.hpp>
#include <CL/opencl.h>
#include <iomanip>

typedef float Float;
typedef std::complex<Float> Complex;
typedef std::vector<Complex> Signal;

static const Float eps = 0.01;
static const Complex i(0, 1);

static bool operator==(cl_float2 a, cl_float2 b)
{
    return a.s[0] == b.s[0] && a.s[1] == b.s[1];
}

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
            c += signal[n]*exp(-i*(Float) 2.0*(Float) M_PI*(Float) k*(Float) n/(Float) signal.size());
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
            c += spectrum[k]*exp(i*(Float) 2.0*(Float) M_PI*(Float) k*(Float) n/(Float) spectrum.size());
        }

        result[n] = c/(Float) spectrum.size();
    }

    return result;
}

static inline Complex W(int k, int N)
{
    return exp(-i*(Float) 2.0*(Float) M_PI*(Float) k/(Float) N);
}

static inline Complex Q(int n, int N)
{
    return exp(i*(Float) 2.0*(Float) M_PI*(Float) n/(Float) N);
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

        spectrum[sample1] = (Float) 0.5*(even + Q(sample1, spectrumSize)*odd);
        spectrum[sample2] = (Float) 0.5*(even + Q(sample2, spectrumSize)*odd);
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

// RMS of error signal.
static Float error(Signal const& a, Signal const& b)
{
    auto e = a - b;
    return sqrt(std::real(dot(e, e))/e.size());
}

static Float prop_inverse_dft(Signal const& test_signal)
{
    return error(test_signal, idft(dft(test_signal)));
}

static Float prop_inverse_fft(Signal const& test_signal)
{
    return error(test_signal, ifft(fft(test_signal)));
}

static Float prop_dft_equal_fft(Signal const& test_signal)
{
    return error(dft(test_signal), fft(test_signal));
}

static Float prop_idft_equal_ifft(Signal const& test_signal)
{
    return error(idft(test_signal), ifft(test_signal));
}

static Float prop_fft_is_decomposed_dft(Signal const& test_signal)
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

static void test_residue(const char* test_name, Float residue)
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
    std::uniform_real_distribution<Float> distribution(0.0, 1.0);

    Signal result(size);
    for (size_t i = 0; i != result.size(); ++i) {
        result[i] = Complex(distribution(generator), distribution(generator));
    }
    return result;
}

static void notify(char const* errinfo, void const* private_info, size_t cb, void* user_data) __attribute__((unused));
static void notify(char const* errinfo, void const* private_info, size_t cb, void* user_data)
{
    std::cout << "OpenCL error: " << errinfo << "\n";
}

static void fatal(std::string const& msg) __attribute__((noreturn));
static void fatal(std::string const& msg)
{
    std::cout << "ERROR: " << msg << "\n";
    abort();
}

static std::string error_code_to_string(cl_int ec) __attribute__((unused));
static std::string error_code_to_string(cl_int ec)
{
    switch (ec) {
        case CL_INVALID_COMMAND_QUEUE: return "CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_CONTEXT: return "CL_INVALID_CONTEXT";
        case CL_INVALID_DEVICE: return "CL_INVALID_DEVICE";
        case CL_INVALID_EVENT_WAIT_LIST: return "CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_GLOBAL_OFFSET: return "CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_GLOBAL_WORK_SIZE: return "CL_INVALID_GLOBAL_WORK_SIZE";
        case CL_INVALID_IMAGE_SIZE: return "CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_KERNEL_ARGS: return "CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_KERNEL: return "CL_INVALID_KERNEL";
        case CL_INVALID_PROGRAM_EXECUTABLE: return "CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_VALUE: return "CL_INVALID_VALUE";
        case CL_INVALID_WORK_DIMENSION: return "CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE: return "CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE: return "CL_INVALID_WORK_ITEM_SIZE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET: return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case CL_OUT_OF_HOST_MEMORY: return "CL_OUT_OF_HOST_MEMORY";
        case CL_OUT_OF_RESOURCES: return "CL_OUT_OF_RESOURCES";
        case CL_SUCCESS: return "CL_SUCCESS";
        default: return "<UNKNOWN>";
    }
}

static std::map<cl_platform_id, std::pair<std::string, std::string>> get_platforms()
{
    std::map<cl_platform_id, std::pair<std::string, std::string>> result;

    std::vector<cl_platform_id> platforms(16);
    cl_uint platform_count;
    if (clGetPlatformIDs(platforms.size(), &platforms[0], &platform_count) != CL_SUCCESS) {
        fatal("Could not get platform IDs.");
    }
    platforms.resize(platform_count);

    for (auto& platform_id : platforms) {
        std::vector<char> platformName(64);
        if (clGetPlatformInfo(
                    platform_id,
                    CL_PLATFORM_NAME,
                    platformName.size(),
                    &platformName[0],
                    NULL) != CL_SUCCESS) {
            fatal("Could not get platform name.");
        }

        std::vector<char> platformVersion(64);
        if (clGetPlatformInfo(
                    platform_id,
                    CL_PLATFORM_VERSION,
                    platformVersion.size(),
                    &platformVersion[0],
                    NULL) != CL_SUCCESS) {
            fatal("Could not get platform version.");
        }

        result[platform_id] = std::make_pair(&platformName[0], &platformVersion[0]);
    }

    return result;
}

static std::map<cl_device_id, std::string> get_devices(cl_platform_id platform)
{
    std::map<cl_device_id, std::string> result;

    std::vector<cl_device_id> devices(16);
    cl_uint device_count;
    if (clGetDeviceIDs(
                platform,
                CL_DEVICE_TYPE_ALL,
                devices.size(),
                &devices[0],
                &device_count) != CL_SUCCESS) {
        fatal("Could not get device IDs.");
    }
    devices.resize(device_count);

    for (auto& device_id : devices) {
        std::vector<char> name(64);
        if (clGetDeviceInfo(
                device_id,
                CL_DEVICE_NAME,
                name.size(),
                &name[0],
                NULL) != CL_SUCCESS) {
            fatal("Could not get device name.");
        }

        result[device_id] = &name[0];
    }

    return result;
}

static void print_platforms()
{
    for (auto& platform : get_platforms()) {
        std::cout
            << platform.first
            << ": name='" << platform.second.first
            << "', version='" << platform.second.second
            << "'\n";

        for (auto& device : get_devices(platform.first)) {
            std::cout
                << "  "
                << device.first
                << ": name='" << device.second
                << "'\n";
        }
    }
}

static cl_platform_id get_platform(std::string const& name)
{
    for (auto& platform : get_platforms()) {
        if (platform.second.first == name) return platform.first;
    }

    fatal("Could not find platform.");
}

static cl_device_id get_device(cl_platform_id platform, std::string const& name)
{
    for (auto& device : get_devices(platform)) {
        if (device.second == name) return device.first;
    }

    fatal("Could not find device.");
}

static std::vector<char> read_program(std::string const& name)
{
    std::fstream f(name);
    if (!f) {
        fatal("Could not read file: " + name);
    }

    f.seekg(0, std::ios_base::end);
    size_t size = f.tellg();
    f.seekg(0, std::ios_base::beg);

    std::vector<char> result(size);

    f.read(&result[0], result.size());
    result.push_back(0);
    return result;
}

static void run_opencl()
{
    cl_int ec;

    print_platforms();

    cl_platform_id platform = get_platform("NVIDIA CUDA");
    cl_device_id device = get_device(platform, "GeForce GTX 550 Ti");

    cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties) platform, 0 };
    cl_context context = clCreateContext(properties, 1, &device, notify, NULL, NULL);
    if (0 == context) fatal("Could not create contex.");

    cl_command_queue queue = clCreateCommandQueue(context, device, 0, NULL);
    if (0 == queue) fatal("Could not create command queue.");

    std::vector<std::vector<char>> sources;
    sources.push_back(read_program("fourier.cl"));

    std::vector<char const*> sources_raw;
    for (auto& source : sources) {
        sources_raw.push_back(&source[0]);
    }

    cl_program program = clCreateProgramWithSource(
            context,
            sources_raw.size(),
            &sources_raw[0],
            NULL,
            NULL);
    if (program == NULL) fatal("Could not create program.");

    if (clBuildProgram(
                program,
                1,
                &device,
                "",
                NULL,
                NULL) != CL_SUCCESS) {
        std::vector<char> build_log(1024);
        size_t build_log_size;
        if (clGetProgramBuildInfo(
                    program,
                    device,
                    CL_PROGRAM_BUILD_LOG,
                    build_log.size(),
                    &build_log[0],
                    &build_log_size)
                != CL_SUCCESS) {
            fatal("Could not get program build info.");
        }
        build_log.resize(build_log_size);

        std::cout << &build_log[0];

        fatal("Could not build program.");
    }

    cl_kernel init_kernel = clCreateKernel(program, "fft_init", NULL);
    if (init_kernel == NULL) fatal("Could not create init kernel.");

    cl_kernel step_kernel = clCreateKernel(program, "fft_step", NULL);
    if (step_kernel == NULL) fatal("Could not create step kernel.");

    int exponent_n = 10;
    cl_float2 dummy;
    dummy.s[0] = 43;
    dummy.s[1] = 37;
    std::vector<cl_float2> x_buffer(1 << exponent_n, dummy);
    size_t byte_count = x_buffer.size()*sizeof(cl_float2);

    cl_mem x_mem = clCreateBuffer(
            context,
            CL_MEM_READ_ONLY,
            byte_count,
            NULL,
            NULL);
    if (x_mem == NULL) fatal("Could not create X buffer.");

    cl_mem y1_mem = clCreateBuffer(
            context,
            CL_MEM_WRITE_ONLY,
            byte_count,
            NULL,
            NULL);
    if (y1_mem == NULL) fatal("Could not create Y1 buffer.");

    cl_mem y2_mem = clCreateBuffer(
            context,
            CL_MEM_WRITE_ONLY,
            byte_count,
            NULL,
            NULL);
    if (y2_mem == NULL) fatal("Could not create Y2 buffer.");

    if (clEnqueueWriteBuffer(
                queue,
                x_mem,
                CL_TRUE,
                0,
                byte_count,
                &x_buffer[0],
                0,
                NULL,
                NULL) != CL_SUCCESS) {
        fatal("Coule not write to X buffer.");
    }

    if (clSetKernelArg(
                init_kernel,
                0,
                sizeof(x_mem),
                &x_mem
                ) != CL_SUCCESS) {
        fatal("Could not set kernel argument 0.");
    }

    cl_uint exponent_n_arg = exponent_n;
    if (clSetKernelArg(
                init_kernel,
                1,
                sizeof(exponent_n_arg),
                &exponent_n_arg
                ) != CL_SUCCESS) {
        fatal("Could not set kernel argument 1.");
    }

    if (clSetKernelArg(
                init_kernel,
                2,
                sizeof(y1_mem),
                &y1_mem
                ) != CL_SUCCESS) {
        fatal("Could not set kernel argument 2.");
    }

    size_t global_work_size = x_buffer.size();
    ec = clEnqueueNDRangeKernel(
                queue,
                init_kernel,
                1,
                NULL,
                &global_work_size,
                NULL,
                0,
                NULL,
                NULL);
    if (ec != CL_SUCCESS) {
        std::cout << error_code_to_string(ec) << "\n";
        fatal("Could not enqueue kernel.");
    }

    cl_mem y_old_mem = y1_mem;
    cl_mem y_new_mem = y2_mem;
    cl_uint B = 2;
    while (B != x_buffer.size()) {
        if (clSetKernelArg(
                    step_kernel,
                    0,
                    sizeof(y_old_mem),
                    &y_old_mem
                    ) != CL_SUCCESS) {
            fatal("Could not set kernel argument 0.");
        }

        if (clSetKernelArg(
                    step_kernel,
                    1,
                    sizeof(B),
                    &B
                    ) != CL_SUCCESS) {
            fatal("Could not set kernel argument 1.");
        }

        if (clSetKernelArg(
                    step_kernel,
                    2,
                    sizeof(y_new_mem),
                    &y_new_mem
                    ) != CL_SUCCESS) {
            fatal("Could not set kernel argument 2.");
        }

        ec = clEnqueueNDRangeKernel(
                queue,
                step_kernel,
                1,
                NULL,
                &global_work_size,
                NULL,
                0,
                NULL,
                NULL);
        if (ec != CL_SUCCESS) {
            std::cout << error_code_to_string(ec) << "\n";
            fatal("Could not enqueue kernel.");
        }

        std::swap(y_old_mem, y_new_mem);
        B <<= 1;
    }

    std::vector<cl_float2> y1_buffer(x_buffer.size());
    ec = clEnqueueReadBuffer(
                queue,
                y_old_mem,
                CL_TRUE,
                0,
                byte_count,
                &y1_buffer[0],
                0,
                NULL,
                NULL);
    if (ec != CL_SUCCESS) {
        std::cout << error_code_to_string(ec) << "\n";
        fatal("Could not read Y1 buffer.");
    }

    if (clFinish(queue) != CL_SUCCESS) fatal("Could not finish.");

    if (x_buffer != y1_buffer) {
        fatal("X buffer is not equal Y1 buffer.");
    }

    if (clReleaseMemObject(y2_mem) != CL_SUCCESS) fatal("Could not release Y2 buffer.");
    if (clReleaseMemObject(y1_mem) != CL_SUCCESS) fatal("Could not release Y1 buffer.");
    if (clReleaseMemObject(x_mem) != CL_SUCCESS) fatal("Could not release X buffer.");
    if (clReleaseKernel(step_kernel) != CL_SUCCESS) fatal("Could not release step kernel.");
    if (clReleaseKernel(init_kernel) != CL_SUCCESS) fatal("Could not release init kernel.");
    if (clReleaseProgram(program) != CL_SUCCESS) fatal("Could not release program");
    if (clUnloadCompiler() != CL_SUCCESS) fatal("Could not unload compiler.");
    if (clReleaseCommandQueue(queue) != CL_SUCCESS) fatal("Could not release command queue.");
    if (clReleaseContext(context) != CL_SUCCESS) fatal("Could not release context.");
}

static void print_reverse_bits_table() __attribute((unused));
static void print_reverse_bits_table()
{
    std::stringstream ss;

    for (size_t row = 0; row != 16; ++row) {
        for (size_t column = 0; column != 16; ++column) {
            ss << std::hex << std::setw(2) << std::setfill('0');
            ss << reverse_bits(row*16+column, 256) << " ";
        }
        ss << "\n";
    }

    std::cout << ss.str();
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

//    print_reverse_bits_table();
    run_opencl();

    return 0;
}
