//
// Created by rowan on 11/05/2024.
//

#ifndef CWT_RICKER_H
#define CWT_RICKER_H

#include <cmath>
#define log2e 1.442695040888963407

class Ricker {
public:
    Ricker() = default;
    [[nodiscard]] inline float RickerInt(float t) const{
        return coefexp * t * expf(t * t);
    }
private:
    const float rickercoef = 0.8673250705840775;
    const float coefexp = 0.5260592472466668754;
    [[nodiscard]] inline float ricker(float t) const {
        return rickercoef * (1.f - t * t) * exp2f(-0.5f * t * t);
    }
};


#endif //CWT_RICKER_H
