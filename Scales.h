//
// Created by rowan on 11/05/2024.
//

#ifndef CWT_SCALES_H
#define CWT_SCALES_H

class Scales {

public:
    const int n;
    float* scale = nullptr;

    Scales(float f0, float f1, int n);
};


#endif //CWT_SCALES_H
