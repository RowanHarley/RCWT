#include "transforms.h"
#include <cmath>
#include "Scales.h"
#include "Ricker.h"

// #define TIME_DEBUGGING
// #define PROFILING
#include <iostream>
#ifndef PROFILING
#include "omp.h"
#endif

inline void FillRickerArray(const float sc, const int N, Ricker* ricker, float* out){
    const float coef = sqrtf(sc);
    for(int i = 0; i < N; i++){
        out[i] = coef * ricker->RickerInt(-8.f + 16.f * (float)(1 + i)/(sc * 16.f));
    }
}

void ShortConvolve(float* x, const int M, Scales* scale, Ricker* ricker, float* out, const int row){
    // Exploits the fact that we only need from floor(d):len(conv) - ceil(d)
    // std::cout << "Size: " << out->size() << std::endl;
    const float sc = scale->scale[row];
    const int N = sc * 16 + 1;
    const int coefLen = (M + N - 1) - 1; // Differencing will remove 1
    const float d = (float)(coefLen - M)/2.f;
    const int dFloor = floorf(d);
    const int len =  coefLen - (int)ceilf(d) - dFloor;
    auto* multF = (float*)malloc(sizeof(float) * N);
    FillRickerArray(sc, N, ricker, multF);
    int k = -dFloor;
    float firstElement = 0.0;
    // int init = -1;
    for(int i = 0; i < M; ++i, ++k){
        const float xi = x[i];
        int tempK = k;
        if (xi == 0.0f) {
            // k = -dFloor + i + 1;
            continue;
        }
        int j = int_max( -1, -k - 1); // Works, not sure why it's -1 here
        if(tempK > 0){
            const float secondElement = out[tempK - 1];
            out[tempK - 1] -= firstElement;
            firstElement = secondElement;
        }else {
            firstElement += xi * multF[j];//GetRicker(j, sc, coef, ricker);
            tempK = 1;
            ++j;
        }

        // for(int j = int_max( -1, -k + 1); j < N; ++j)
        while(j < N && tempK < len + 1){
            const int idx = tempK - 1;
            const float temp = multF[j];
            out[idx] += xi * temp;
            ++j;
            ++tempK;
        }

        /*for(; j < N && tempK < len + 1; ++j, ++tempK){
            const int idx = tempK - 1;
            const float temp = multF[j];
            out[idx] += xi * temp;
            //if (tempK > 0) {
            *//*if (idx > init) {
                 //GetRicker(j, sc, coef, ricker);
                init = idx;
            } else {
                out[idx] += xi * temp;//GetRicker(j, sc, coef, ricker);//xi * multF[j];
            }*//*
            //}else
            //    firstElement += xi * temp;//GetRicker(j, sc, coef, ricker);


        }*/
        // k = -dFloor + i + 1;
    }

    if(k < M){
        for(k = int_max(k, 1);k <= M; ++k){
            const int idx = k - 1;
            const float secondElement = out[idx];
            out[idx] -= firstElement;
            firstElement = secondElement;
        }
    }
}
void cwt(std::vector<float>& x, Scales* scales, Ricker* ricker, float* out, int threads, int ndim){
    const int X = (int)x.size();
#ifndef PROFILING
    omp_set_num_threads( threads );                     // OpenMP
#pragma omp parallel default(none) shared(scales, ricker, out, x, X)
    {
#pragma omp for nowait schedule(dynamic)
#endif
    for (int i = 0; i < scales->n; ++i) {
        ShortConvolve(&x[0], X, scales, ricker, &out[i * X], i);
    }
#ifndef PROFILING
    }
#endif
}
