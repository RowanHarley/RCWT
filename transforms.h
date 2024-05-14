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

void cwt(std::vector<float>& x, Scales* scales, Ricker* ricker, float* out, int threads, int ndim);

#endif //CWT_TRANSFORMS_H
