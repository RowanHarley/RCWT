#include "transforms.h"
#include <cmath>

// #define TIME_DEBUGGING

#ifdef TIME_DEBUGGING
#include <iostream>
#endif
#include <chrono>
#include "omp.h"

using namespace std::chrono;

std::vector<double> approximateRicker(int samples){
    double stepSize = 1.0 * (upper_bound - lower_bound)/samples;
    std::vector<double> int_psi(samples);
    double coef = 1.0 * (upper_bound - lower_bound)/(2 * samples);
    int_psi[0] = coef * ricker(static_cast<double>(lower_bound));
    for(int i = 1; i < samples - 1; i++){
        int_psi[i] = int_psi[i - 1] + 2 * coef * ricker(static_cast<double>(lower_bound) + i * stepSize);
    }
    int_psi[samples - 1] = int_psi[samples - 2] + coef * ricker(static_cast<double>(upper_bound));
    return int_psi;
}
std::vector<double> stepSize(int noSteps){
    double stepSize = (upper_bound - lower_bound)/static_cast<double>(noSteps);
    std::vector<double> steps(noSteps);
    for(int i = 0; i < noSteps; i++){
        steps[i] = lower_bound + i * stepSize;
    }
    return steps;
}
/*std::vector<int> CreateJ(double scale, double step, const std::vector<double>& x){
    int range = scale * (x.back() - x.front()) + 1; // Truncation here. To ensure correct, we add 1 in loop
    std::vector<int> j(range + 1);
    for(int i = 0; i < range + 1; i++){
        j[i] = i/(scale * step);
    }
    if(j.back() >= x.size()){
        j.erase(std::remove_if(j.begin(), j.end(),
                               [&x](int a) {
            return a >= x.size();
        }), j.end());
    }
    return j;
}*/
std::vector<double> CreatePsiScale(const double& scale, const double& step, const std::vector<double>& int_psi, const int& N){
    // Using 16 here as we have upper and lower bounds of +- 8
    int range = scale * 16 + 1; // Truncation here. To ensure correct, we add 1 in loop

    // int len = ((range - 1)/(scale * step) >= N ? floor(N*(scale * step)) - 1 : range);
    std::vector<double> psiSc(range + 1);
    // std::cout << "Len: "<< len << ", N: " << N << std::endl;
    for(int i = range; i > -1; i--){
        int j = i/(scale * step);
        if (j < N)
            psiSc[range - i] = int_psi[j];
    }
    return psiSc;
}
inline int next_fast_len(int& n){
    return powf(2.0f, ceilf(log2f((float)n)));
}

std::vector<double> CreateIntPsiScale(const std::vector<int>& j, const std::vector<double>& int_psi){
    std::vector<double> int_psi_scale(j.size());
    for(int i = 0; i < j.size(); i++){
        int_psi_scale[i] = int_psi[j[j.size() - i - 1]];
    }
    return int_psi_scale;
}
std::vector<double> convolve(std::vector<double>& x, const std::vector<double>& y){
    const int N = static_cast<int>(y.size());
    std::vector<double> convolution(x.size() + y.size() - 1);
    // std::reverse(x.begin(), x.end());
    int k = 0;
    for(double i : x){
        for(double j : y){
            convolution[k] += i * j;
            k++;
        }
        k -= N - 1;
    }
    return convolution;
}
std::vector<double> shortConvolve(const std::vector<double>& x, const std::vector<double>& y){
    // Exploits the fact that we only need from floor(d):len(conv) - ceil(d)
    const int N = static_cast<int>(y.size());
    const int M = static_cast<int>(x.size());
    const int coefLen = (M + N - 1) - 1; // Differencing will remove 1
    const double d = (coefLen - M)/2.;
    const int dFloor = floor(d);
    const int len =  coefLen - ceil(d) - dFloor;

    std::vector<double> convolution(len + 1); // Add 1 back as we need to remove it when differencing
    // std::reverse(x.begin(), x.end());
    // double MeanTime = 0.0;
    int k = -dFloor;
    for(int i = 0; i < M; i++){
        /*int kDiff = -k;
        if (kDiff >= 0) {
            for (int j = kDiff; j < N; j++) {
                convolution[j - kDiff] += x[i] * y[j];
            }
        }else{
            for (double j : y) {
                if(k >= len + 1)
                    break;
                convolution[k] += x[i] * j;
                k++;
            }
        }*/
        auto start = high_resolution_clock::now();
        for(double j : y){
            if(k >= len + 1)
                break;
            if(k >= 0)
                convolution[k] += x[i] * j;
            k++;
        }
        /*k = -dFloor + i + 1;
        if (k >= len + 1)
            break;
        auto finish = high_resolution_clock::now();
        auto diff = finish - start;
        MeanTime += (double)diff.count()/M;
        if(diff.count() > 5000)
            std::cout << "Loop took " << diff.count() << " units." << std::endl;*/
    }
    // std::cout << "Took average time of " << MeanTime << " units. M = " << M << ", N = " << N <<  std::endl;
    return convolution;
}
void shortConvolveArray(double *x, const int& M, const std::vector<double>& data, const std::vector<double>& y){
    const int N = static_cast<int>(y.size());
    const int coefLen = (M + N - 1) - 1; // Differencing will remove 1
    const double d = (coefLen - M)/2.;
    const int len =  coefLen - ceil(d) - floor(d);
    // std::reverse(x.begin(), x.end());
    int k = -floor(d);
    int kFilled = -1;
    for(int i = 0; i < M - 1; i++){
        for(double j : y){
            if(k >= len + 1)
                break;
            if(k >= 0 && k <= kFilled)
                x[k] += data[i] * j;
            else if(k >= 0){
                x[k] = data[i] * j;
                kFilled++;
            }
            k++;
        }
        k = i + 1 - floor(d);
        if (k >= len + 1)
            break;
    }
}

std::vector<double> operator*(const std::vector<double> &b, const double &a) {
    std::vector<double> c(b.size());
    for (int i = 0; i < b.size(); ++i)
        c[i] = b[i] / a;
    return c;
}
std::vector<double> diff(const std::vector<double>& x){
    std::vector<double> diffX(x.size() - 1);
    for(int i = 1; i < x.size(); i++){
        diffX[i - 1] = x[i] - x[i-1];
    }
    return diffX;
}
std::vector<double> getCoef(const std::vector<double>& x, const double& scale){
    std::vector<double> diffX(x.size() - 1);
    const double coef = -sqrt(scale);
    for(int i = 1; i < x.size(); i++){
        diffX[i - 1] = coef * (x[i] - x[i-1]);
    }
    return diffX;
}
std::vector<double> getCoefArr(double *x, const int& N, const double& scale){
    std::vector<double> diffX(N - 1);
    const double coef = -sqrt(scale);
    for(int i = 1; i < N; i++){

        diffX[i - 1] = coef * (x[i] - x[i-1]);
    }
    return diffX;
}

std::vector<std::vector<double>> cwt(std::vector<double>& x, std::vector<double>& scales, int samples, int threads = 1) {
    const int N = static_cast<int>(scales.size());
    std::vector<std::vector<double>> out(N);
    const std::vector<double> int_psi = approximateRicker(samples);
    const double step = 16./samples;

    omp_set_num_threads( threads );                     // OpenMP
    #pragma omp parallel for
    for(int i = 0; i < scales.size(); i++){
#ifdef TIME_DEBUGGING
        auto start = high_resolution_clock::now();
#endif
        const int range = scales[i] * 16 + 1;
        std::vector<double> int_psi_scale(range + 1);
        int_psi_scale = CreatePsiScale(scales[i], step, int_psi, samples);
#ifdef TIME_DEBUGGING
        auto finish = high_resolution_clock::now();
        auto totTime = finish - start;
        std::cout << "Took " << totTime.count() << "s to create psi scale" << std::endl;
        start = high_resolution_clock::now();
#endif
        // double *conv = (double*)malloc(sizeof(double)*(M + 1));
        std::vector<double> conv(N + 1);
        conv = shortConvolve(x , int_psi_scale);
#ifdef TIME_DEBUGGING
        finish = high_resolution_clock::now();
        totTime = finish - start;
        std::cout << "Took " << totTime.count() << "s to create convolve!" << std::endl;
        start = high_resolution_clock::now();
#endif
        // coef = diff(conv) * - (sqrt(scales[i])/(i*scalesStep + 1)); // So dividing by i*(Scale Step) + 1 solves the mismatch (I have no clue why)

#ifdef TIME_DEBUGGING
        finish = high_resolution_clock::now();
        totTime = finish - start;
        std::cout << "Took " << totTime.count() << "s to multiply vector by const." << std::endl;
#endif
        out[i] = getCoef(conv, scales[i]);

    }
    return out;
}
