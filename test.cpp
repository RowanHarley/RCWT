//
// Created by rowan on 04/05/2024.
//
// #include "cwt.h"
#include "transforms.h"
#include <cmath>
#include <iostream>
#include <iomanip>

std::vector<float> testY(){
    std::vector<float> y(15);
    // std::cout << "[ ";
    for(int i = 0; i < y.size(); i++){
        y[i] = sin(2 * M_PI * i/32);
        // std::cout << y[i] << " ";
    }
    // std::cout << "]" << std::endl;
    return y;
}
std::vector<float> testX(){
    std::vector<float> y(15);
    // std::cout << "[ ";
    for(int i = 0; i < y.size(); i++){
        y[i] = 2 * sin(2 * M_PI * i/64 - 1);
        // std::cout << y[i] << " ";
    }
    // std::cout << "]" << std::endl;
    return y;
}
int main(){

    std::vector<float> y = testY();
    std::vector<float> scales(30);
    for(int i = 2; i <= 31; i++){
        scales[i-2] = i/2.;
    }
    std::vector<float> y2 = testX();
    Scales scale(1, 16, 30);
    Ricker ricker(1000);

    std::vector<float> res(30 * 15);
    cwt(y, &scale, &ricker, &res[0], 5, 1);
    // auto res = CWT.transform(1000);
    std::cout << std::fixed;
    std::cout << std::setprecision(3);
    std::cout << "[ ";
    for(int i = 0; i < 30 * 15; i++){
        std::cout << res[i] << " ";
        if(i % 15 == 14) {
            std::cout << "]" << std::endl;
            if(i != 30 * 15 - 1)
                std::cout << "[ ";
        }
    }
    return 0;
}