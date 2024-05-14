#ifndef CWT_TRANSFORMS_H
#define CWT_TRANSFORMS_H

#include <vector>
#include <cmath>
#include "Ricker.h"
#include "Scales.h"

const float rickercoef = 0.8673250705840775;
const int lower_bound = -8;
const int upper_bound = 8;


inline int int_max(int a, int b) { return (a >= b ? a : b); }
inline int int_min(int a, int b) { return (a <= b ? a : b); }

inline float ricker(float t) {
    return rickercoef * (1 - t * t) * exp(-0.5 * t * t);
}
std::vector<float> CreatePsiScale(const float& scale, const float& step, const std::vector<float>& int_psi, const int& N);
std::vector<float> convolve(std::vector<float>& x, const std::vector<float>& y);
void cwt(std::vector<float>& x, Scales* scales, Ricker* ricker, float* out, int threads, int ndim);
void cwt(std::vector<float>& x, std::vector<float>& scales, std::vector<float>& out, int samples, int threads = 1, int ndim = 1);
std::vector<std::vector<std::vector<float>>> cwt(std::vector<std::vector<float>>& x, std::vector<float>& scales, int samples, int threads);
std::vector<float> approximateRicker(int precision);
#endif //CWT_TRANSFORMS_H
