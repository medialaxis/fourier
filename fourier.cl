// vim:filetype=opencl

typedef float2 Complex;

static uint index(uint n, uint B, uint k)
{
    return n*B + k;
}

static Complex W(uint k, uint K)
{
    float angle = -2.0*M_PI_F*(float) k/(float) K;
    return (cos(angle), sin(angle));
}

static Complex mult(Complex a, Complex b)
{
    return (a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

kernel void fft_step(Complex global const* Y, uint B, Complex global* Y_)
{
    uint i = get_global_id(0);
    uint B_ = B/2;
    uint n_ = i/B_;
    uint k_ = i%B_;

    Y_[index(n_, B_, k_)] = Y[index(n_*2, B, k_%B)] + mult(W(k_, B_), Y[index(n_*2+1, B, k_%B)]);
}
