#ifndef CWT_TRANSFORMS_H
#define CWT_TRANSFORMS_H

#include <vector>
#include <cmath>
const double rickercoef = 0.8673250705840775;
const int lower_bound = -8;
const int upper_bound = 8;


inline int int_max(int a, int b) { return (a >= b ? a : b); }
inline int int_min(int a, int b) { return (a <= b ? a : b); }

inline double ricker(double t) {
    return rickercoef * (1 - t * t) * exp(-0.5 * t * t);
}
std::vector<double> CreatePsiScale(const double& scale, const double& step, const std::vector<double>& int_psi, const int& N);
std::vector<double> convolve(std::vector<double>& x, const std::vector<double>& y);
std::vector<std::vector<double>> cwt(std::vector<double>& x, std::vector<double>& scales, int samples, int threads);
std::vector<double> approximateRicker(int precision);
std::vector<std::vector<double>> cwt(std::vector<double>& x, std::vector<int>& scales);
#endif //CWT_TRANSFORMS_H
