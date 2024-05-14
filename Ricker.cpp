//
// Created by rowan on 11/05/2024.
//

#include <cstdlib>
#include "Ricker.h"

/*
Ricker::Ricker(int s) : stepSize(8.0f/(float)s), samples(s) {
    float coef = stepSize/2.f;
    mother = (float*)malloc(sizeof(float) * samples);
    mother[0] = coef * ricker(-8.0);
    for(int i = 1; i < samples; i++){ // 1 sided so we can use the halfway point and mirror
        mother[i] = mother[i - 1] + 2 * coef * ricker(-8.0f + (float)i * stepSize);
    }
    mother[samples - 1] = mother[samples - 2] + coef * ricker(0.f);
}*/
