//
// Created by rowan on 04/05/2024.
//
// #include "cwt.h"
#include "transforms.h"
#include <cmath>
#include <iostream>
#include <iomanip>

std::vector<double> testY(){
    std::vector<double> y(15);
    // std::cout << "[ ";
    for(int i = 0; i < y.size(); i++){
        y[i] = sin(2 * M_PI * i/32);
        // std::cout << y[i] << " ";
    }
    // std::cout << "]" << std::endl;
    return y;
}

int main(){

    std::vector<double> y = testY();
    std::vector<double> scales(30);
    for(int i = 2; i <= 31; i++){
        scales[i-2] = i/2.;
    }
    auto res = cwt(y, scales, 1000, 5);
    // auto res = CWT.transform(1000);
    std::cout << std::fixed;
    std::cout << std::setprecision(3);
    for(auto & re : res){
        std::cout << "[ ";
        for(double j : re){
            std::cout << j << " ";
        }
        std::cout << "]" << std::endl;
    }
}