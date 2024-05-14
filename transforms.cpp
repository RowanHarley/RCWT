#include "transforms.h"
#include <cmath>
#include "Scales.h"
#include "Ricker.h"

// #define TIME_DEBUGGING
// #define PROFILING
#include <iostream>
#include <cstdint>
#ifndef PROFILING
#include "omp.h"
#endif



std::vector<float> approximateRicker(int samples){
    float stepSize = 16.f/(float)samples;
    std::vector<float> int_psi(samples);
    float coef = 1.0f * 16.f/(2.0f * (float)samples);
    int_psi[0] = coef * ricker(-8.f);
    for(int i = 1; i < samples - 1; i++){
        int_psi[i] = int_psi[i - 1] + 2 * coef * ricker(-8.f + (float)i * stepSize);
    }
    int_psi[samples - 1] = int_psi[samples - 2] + coef * ricker(static_cast<float>(upper_bound));
    return int_psi;
}
float approximateRickerArr(float* out, int samples){
    const float stepSize = 16.0f/(float)samples;
    float coef = 16.0f/(float)(2 * samples);
    out[0] = coef * ricker(-8.0);
    for(int i = 1; i < samples - 1; i++){
        out[i] = out[i - 1] + 2 * coef * ricker(-8.0f + (float)i * stepSize);
    }
    out[samples - 1] = out[samples - 2] + coef * ricker(8.0);
    return stepSize;
}

std::vector<float> stepSize(int noSteps){
    float stepSize = (upper_bound - lower_bound)/static_cast<float>(noSteps);
    std::vector<float> steps(noSteps);
    for(int i = 0; i < noSteps; i++){
        steps[i] = lower_bound + (float)i * stepSize;
    }
    return steps;
}
std::vector<float> CreatePsiScale(const float& scale, const float& step, const std::vector<float>& int_psi, const int& N){
    // Using 16 here as we have upper and lower bounds of +- 8
    int range = scale * 16 + 1; // Truncation here. To ensure correct, we add 1 in loop

    // int len = ((range - 1)/(scale * step) >= N ? floor(N*(scale * step)) - 1 : range);
    std::vector<float> psiSc(range + 1);
    // std::cout << "Len: "<< len << ", N: " << N << std::endl;
    for(int i = range; i > -1; i--){
        int j = i/(scale * step);
        if (j < N)
            psiSc[range - i] = int_psi[j];
    }
    return psiSc;
}
void CreatePsi(const float& scale, const float& step, float *out, float* int_psi, const int& N){
    // Using 16 here as we have upper and lower bounds of +- 8
    int range = scale * 16; // Truncation here
    // First row will take 16/(16/Samples) = Samples, 15/16 * Samples, ...
    // Second row will take 32/(2 * 16/Samples) = Samples, 31/32 * Samples, ...
    // In general, we get samples from Samples to (Scale * 16 - 1)/(Scale * 16) * Samples,..., 0
    for(int i = range; i > -1; i--){
        int j = i/(scale * step);
        if (j < N) // This should be always true except when floor(scale * 16) = scale * 16
            out[range - i] = int_psi[j]; // Reverse the int_psi
    }
}
void JointPsiConvolve(const std::vector<float>& x, const std::vector<float>& scales, std::vector<float>& out, const float* int_psi, const float& step, const int& samples){
    const int M = (int)x.size();
    const int SL = (int)scales.size();
    const int OL = M * SL;
    int ConvL;
    short int row = -1; // Initialise at -1 since the for loop will increment on first run
    float d;
    int dFloor;
    int range = scales[0] * 16;
    int* N = &range;
    int coefLen;
    float* coefs = (float*)malloc(sizeof(float) * SL);
    float firstElement = 0.0;
    int k;
    int i;
    for(i = 0; i < SL; ++i){
        coefs[i] = -sqrtf(scales[i]);
    }
    for(i = 0; i <= OL; ++i){
        short int modM = i % M;
        if(modM == 0){
            if(k < M){
                for(k = int_max(k, 1);k <= M; k++){
                    float secondElement = out[row * M + k - 1];
                    out[row * M + k - 1] -= firstElement;
                    firstElement = secondElement;
                }
            }
            if (i == OL) {
                return;
            }
            firstElement = 0.0;
            row++;
            range = scales[row] * 16;
            if((int)(*N/(scales[row] * step)) >= samples)
                range--;

            coefLen = (M + *N - 1) - 1;
            d = (coefLen - M)/2.f;
            dFloor = floorf(d);
            ConvL = coefLen - (int)ceilf(d) - dFloor;
            k = -dFloor;

        }
        if(k > 0){
            float secondElement = out[row * M + k - 1];
            out[row * M + k - 1] -= firstElement;
            firstElement = secondElement;
            // out[row * M + k - 1] += coefs[row] * x[modM] * ;
        }
        for(int j = 0; j < *N; j++){
            int psi_val =  (1 + j)/(scales[row] * step);
            if(k >= ConvL + 1)
                break;
            else if(k > 0){
                out[row * M + k - 1] += coefs[row] * x[modM] * -int_psi[psi_val];
            }
            else if(k == 0)
                firstElement += coefs[row] * x[modM] * -int_psi[psi_val];
            k++;
        }
        k = -dFloor + modM + 1;

    }
}
std::vector<float> convolve(std::vector<float>& x, const std::vector<float>& y){
    const int N = static_cast<int>(y.size());
    std::vector<float> convolution(x.size() + y.size() - 1);
    // std::reverse(x.begin(), x.end());
    int k = 0;
    for(float i : x){
        for(float j : y){
            convolution[k] += i * j;
            k++;
        }
        k -= N - 1;
    }
    return convolution;
}
std::vector<float> shortConvolve(const std::vector<float>& x, const std::vector<float>& y){
    // Exploits the fact that we only need from floor(d):len(conv) - ceil(d)
    const int N = static_cast<int>(y.size());
    const int M = static_cast<int>(x.size());
    const int coefLen = (M + N - 1) - 1; // Differencing will remove 1
    const float d = (coefLen - M)/2.;
    const int dFloor = floor(d);
    const int len =  coefLen - ceil(d) - dFloor;
    // const float coef = -sqrt(scale);
    std::vector<float> convolution(len); // Add 1 back as we need to remove it when differencing
    // std::reverse(x.begin(), x.end());
    // float MeanTime = 0.0;
    int k = -dFloor;
    // float firstElement = 0.0;
    for(int i = 0; i < M; i++){
        /*if(k > 0){
            float secondElement = convolution[k - 1];
            convolution[k - 1] -= firstElement;
            firstElement = secondElement;
        }*/
        for(float j : y){
            if(k >= len + 1)
                break;
            if(k >= 0)
                convolution[k] += x[i] * j;
            /*else if(k == 0)
                firstElement += coef * x[i] * j;*/

            k++;
        }
        k = -dFloor + i + 1;
    }
    /*if(k < M){
        for(k = int_max(k, 1);k <= M; k++){
            float secondElement = convolution[k - 1];
            convolution[k - 1] -= firstElement;
            firstElement = secondElement;
        }
    }*/
    // std::cout << "Took average time of " << MeanTime << " units. M = " << M << ", N = " << N <<  std::endl;
    return convolution;
}

void shortConvolveCoefArr(const std::vector<float>& x, float* y, const int N, const float scale, std::vector<float>& out, const int row, const float step){
    // Exploits the fact that we only need from floor(d):len(conv) - ceil(d)
    // std::cout << "Size: " << out->size() << std::endl;
    const int M = (int)x.size();
    const int coefLen = (M + N - 1) - 1; // Differencing will remove 1
    const float d = (coefLen - M)/2.;
    const int dFloor = floorf(d);
    const int len =  coefLen - (int)ceilf(d) - dFloor;
    const float coef = -sqrtf(scale);
    // std::vector<float> convolution(len); // Add 1 back as we need to remove it when differencing
    // std::reverse(x.begin(), x.end());
    // float MeanTime = 0.0;
    int k = -dFloor;
    float firstElement = 0.0;
    float* currRow = (float*)malloc(sizeof(float) * M);
    int init = -1;
    for(int i = 0; i < M; i++){
        if(k > 0){
            float secondElement = currRow[k - 1]; //out[row * M + k - 1];

            currRow[k - 1] -= firstElement;
            //out[row * M + k - 1] -= firstElement;
            firstElement = secondElement;
        }
        if(k < 0)
            k = 0;
        const float xi = x[i];
        for(int j = int_max(0, dFloor - i - 1); j < N; j++){
            float tmp = (1 + j)/(scale * step);
            if(k > 0){
                if(k > init) {
                    currRow[k - 1] = coef * xi * -y[(int) tmp];
                    ++init;
                }else{
                    currRow[k - 1] += coef * xi * -y[(int) tmp];
                }
                // out[row * M + k - 1] += coef * x[i] * -y[(int)tmp];
            }
            else
                firstElement += coef * x[i] * -y[(int)tmp];
            k++;

        }
        k = -dFloor + i + 1;
    }
    if(k < M){
        for(k = int_max(k, 1);k <= M; k++){
            float secondElement = currRow[k - 1]; //out[row * M + k - 1];
            currRow[k - 1] -= firstElement;
            /*float secondElement = out[row * M + k - 1];
            out[row * M + k - 1] -= firstElement;*/
            firstElement = secondElement;
        }
    }
    out.insert(out.begin() + row * M, currRow, currRow + M);
    // std::cout << "Done 1 conv" << std::endl;
}
inline float GetRicker(const int j, const float sc, const float coef, Ricker* ricker) {
    // const float div = sc * stepSize;
    // int idx = (int)((1.f + (float)j)/div);
    return - coef * ricker->RickerInt(-8.f + 16.f * (float)(1 + j)/(sc * 16.f));
}
inline void FillRickerArray(const float sc, const float coef, const int N, Ricker* ricker, float* out){
    for(int i = 0; i < N; i++){
        out[i] = - coef * ricker->RickerInt(-8.f + 16.f * (float)(1 + i)/(sc * 16.f));
    }
}

/*inline void IdxRange(const float sc, const float stepSize, const int N, const float coef, Ricker* ricker, float* out){
    const float div = sc * stepSize;
    const int samps = ricker->samples;
    for(int i = 0; i < N; i++){
        auto tmp = (int)((1.f + (float)i)/div);
        if(tmp < samps)
            out[i] = - coef * ricker->mother[tmp];
        else
            out[i] = coef * ricker->mother[2 * samps - tmp];
    }
}*/

void ShortConvolve(const std::vector<float>& x, Scales* scale, Ricker* ricker, float* out, const int row){
    // Exploits the fact that we only need from floor(d):len(conv) - ceil(d)
    // std::cout << "Size: " << out->size() << std::endl;
    const float sc = scale->scale[row];
    const int N = sc * 16 + 1;
    const int M = (int)x.size();
    const int coefLen = (M + N - 1) - 1; // Differencing will remove 1
    const float d = (float)(coefLen - M)/2.f;
    const int dFloor = floorf(d);
    const int len =  coefLen - (int)ceilf(d) - dFloor;
    const float coef = -sqrtf(sc);
    auto* multF = (float*)malloc(sizeof(float) * N);
    // IdxRange(sc, SS, N, coef, ricker, multF);
    FillRickerArray(sc, coef, N, ricker, multF);
    // std::reverse(x.begin(), x.end());
    // float MeanTime = 0.0;

    int k = -dFloor;
    float firstElement = 0.0;
    // float* currRow = (float*)malloc(sizeof(float) * M);
    int init = -1;
    for(int i = 0; i < M; ++i){
        const float xi = x[i];
        if(k > 0){
            /*float secondElement = currRow[k - 1]; //out[row * M + k - 1];
            currRow[k - 1] -= firstElement;*/
            const float secondElement = out[row * M + k - 1];
            // std::cout << "Taking " << firstElement << " from " << out[row * M + k - 1] << std::endl;
            out[row * M + k - 1] -= firstElement;
            firstElement = secondElement;
        }
        if(k < 0)
            k = 0;

        for(int j = int_max( -1, dFloor - i - 1); j < N; ++j){
            // const int tmp = ScaleIdx(j, sc, SS);
            const int idx = row * M + k - 1;
            // const float mother = -ricker->mother[tmp];
            // const float o = coef * xi * mother;
                /*if(k > init) {
                    currRow[k - 1] = coef * x[i] * -ricker->mother[(int) tmp];
                    ++init;
                }else{
                    currRow[k - 1] += coef * x[i] * -ricker->mother[(int) tmp];
                }*/
            if (k > 0) {
                if (idx > init) {
                    out[idx] = xi * multF[j];//GetRicker(j, sc, coef, ricker);
                    init = idx;
                } else {
                    out[idx] += xi * multF[j];//GetRicker(j, sc, coef, ricker);//xi * multF[j];
                }
            }else
                firstElement += xi * multF[j];//GetRicker(j, sc, coef, ricker);

            if(k - 1 >= len)
                break;
            else
                ++k;
        }
        k = -dFloor + i + 1;
    }

    if(k < M){
        for(k = int_max(k, 1);k <= M; ++k){
            const int idx = row * M + k - 1;
                /*float secondElement = currRow[k - 1]; //out[row * M + k - 1];
                currRow[k - 1] -= firstElement;*/
            const float secondElement = out[idx];
            out[idx] -= firstElement;
            firstElement = secondElement;
        }
    }
    // out.insert(out.begin() + row * M, currRow, currRow + M);
    // delete currRow;
    // std::cout << "Done 1 conv" << std::endl;
}


std::vector<float> operator*(const std::vector<float> &b, const float &a) {
    std::vector<float> c(b.size());
    for (int i = 0; i < b.size(); ++i)
        c[i] = b[i] / a;
    return c;
}
std::vector<float> getCoef(const std::vector<float>& x, const float& scale){
    std::vector<float> diffX(x.size() - 1);
    const float coef = -sqrtf(scale);
    for(int i = 1; i < x.size(); i++){
        diffX[i - 1] = coef * (x[i] - x[i-1]);
    }
    return diffX;
}
std::vector<float> getCoefArr(float *x, const int& N, const float& scale){
    std::vector<float> diffX(N - 1);
    const float coef = -sqrtf(scale);
    for(int i = 1; i < N; i++){

        diffX[i - 1] = coef * (x[i] - x[i-1]);
    }
    return diffX;
}

void cwt(std::vector<float>& x, Scales* scales, Ricker* ricker, float* out, int threads, int ndim){
#ifndef PROFILING
    omp_set_num_threads( threads );                     // OpenMP
#pragma omp parallel default(none) shared(scales, ricker, out, x)
    {
#pragma omp for nowait schedule(dynamic)
#endif
    for (int i = 0; i < scales->n; ++i) {
        // std::cout << i << std::endl;
#ifdef TIME_DEBUGGING
        auto start = high_resolution_clock::now();
#endif
        // float *int_psi_sc = (float *) malloc(sizeof(float) * (range + 1));

        // CreatePsi(scales[i], step, int_psi_sc, int_psi, samples);
        ShortConvolve(x, scales, ricker, out, i);
        // shortConvolveCoefArr(x, ricker->mother, range + 1, scales->scale[i], out, i, ricker->stepSize);
        // delete int_psi_sc;
        continue;
        // std::vector<float> int_psi_scale(range + 1);


        // int_psi_scale = CreatePsiScale(scales[i], step, int_psi, samples);
#ifdef TIME_DEBUGGING
        auto finish = high_resolution_clock::now();
            auto totTime = finish - start;
            std::cout << "Took " << totTime.count() << "s to create psi scale" << std::endl;
            start = high_resolution_clock::now();
#endif
        // shortConvolveCoef(x, int_psi_scale, scales[i], out, i);
#ifdef TIME_DEBUGGING
        finish = high_resolution_clock::now();
            totTime = finish - start;
            std::cout << "Took " << totTime.count() << "s to create convolve!" << std::endl;
#endif
    }
#ifndef PROFILING
    }
#endif
}


void cwt(std::vector<float>& x, std::vector<float>& scales, std::vector<float>& out, int samples, int threads, int ndim) {
    float* int_psi = (float*)malloc(sizeof(float) * samples);
    const float step = approximateRickerArr(int_psi, samples);
    //const std::vector<float> int_psi = approximateRicker(samples);
    //JointPsiConvolve(x, scales, out, int_psi, step, samples);
    //delete int_psi;
#ifndef PROFILING
    omp_set_num_threads( threads );                     // OpenMP
#pragma omp parallel default(none) shared(scales, step, int_psi, samples, out, x)
    {
#pragma omp for nowait
#endif
        for (int i = 0; i < scales.size(); ++i) {
            // std::cout << i << std::endl;
#ifdef TIME_DEBUGGING
            auto start = high_resolution_clock::now();
#endif
            const int range = scales[i] * 16 + 1;
            // float *int_psi_sc = (float *) malloc(sizeof(float) * (range + 1));

            // CreatePsi(scales[i], step, int_psi_sc, int_psi, samples);
            shortConvolveCoefArr(x, int_psi, range + 1, scales[i], out, i, step);
            // delete int_psi_sc;
            // std::vector<float> int_psi_scale(range + 1);


            // int_psi_scale = CreatePsiScale(scales[i], step, int_psi, samples);
#ifdef TIME_DEBUGGING
            auto finish = high_resolution_clock::now();
            auto totTime = finish - start;
            std::cout << "Took " << totTime.count() << "s to create psi scale" << std::endl;
            start = high_resolution_clock::now();
#endif
            // shortConvolveCoef(x, int_psi_scale, scales[i], out, i);
#ifdef TIME_DEBUGGING
            finish = high_resolution_clock::now();
            totTime = finish - start;
            std::cout << "Took " << totTime.count() << "s to create convolve!" << std::endl;
#endif
        }
#ifndef PROFILING
    }
#endif
}
std::vector<std::vector<std::vector<float>>> cwt(std::vector<std::vector<float>>& x, std::vector<float>& scales, int samples, int threads = 1) {
    const int N = static_cast<int>(scales.size());
    std::vector<std::vector<std::vector<float>>> out(N);
    const std::vector<float> int_psi = approximateRicker(samples);
    const float step = 16./samples;
    const int batchSize = N/threads;
    const bool roundedDown = N - threads * (int)(N/threads) != 0;

    // omp_set_num_threads( threads );                     // OpenMP
    // #pragma omp parallel for
    for(int th = 0; th < threads; th++) {
        int start = batchSize * th;
        int end = batchSize * (th + 1);
        if (th == threads - 1 && roundedDown)
            end = N;
        for (int i = start; i < end; i++) {
#ifdef TIME_DEBUGGING
            auto start = high_resolution_clock::now();
#endif
            const int range = scales[i] * 16;
            std::vector<float> int_psi_scale(range + 1);
            int_psi_scale = CreatePsiScale(scales[i], step, int_psi, samples);
#ifdef TIME_DEBUGGING
            auto finish = high_resolution_clock::now();
            auto totTime = finish - start;
            std::cout << "Took " << totTime.count() << "s to create psi scale" << std::endl;
            start = high_resolution_clock::now();
#endif
            // float *conv = (float*)malloc(sizeof(float)*(M + 1));
            std::vector<std::vector<float>> conv(N + 1);
            for (int j = 0; j < x[0].size(); j++) {
                conv[j] = shortConvolve(x[j], int_psi_scale);
                out[i][j] = getCoef(conv[j], scales[i]);
            }

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
        }
    }
    return out;
}