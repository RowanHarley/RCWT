//
// Created by rowan on 11/05/2024.
//

#include <cstdlib>
#include "Scales.h"

Scales::Scales(float f0, float f1, int nS): n(nS){
    scale = (float*)malloc(sizeof(float) * n);
    for(int i = 0; i < n; i++){
        scale[i] = f0 + i * (float)(f1 - f0)/n;
    }
}