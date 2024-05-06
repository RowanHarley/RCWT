//
// Created by rowan on 05/05/2024.
//
#include <iostream>
#include "../transforms.h"
#include <chrono>
#include <thread>
using namespace std;

//Calculate and print mean and variance of times array
static void show_stats(chrono::duration<double> *times, int runs)
{
    double total = 0.0;
    double mean = 0.0;
    double std = 0.0;
    for(int i=0; i<runs; i++) {
        total += times[i].count();
    }
    mean = total/runs;
    total = 0.0;

    for(int i=0; i<runs; i++) {
        total += (times[i].count() - mean)*(times[i].count() - mean);
    }
    std = sqrt(total/(runs-1));

    cout << " | elapsed avg time: " << mean << "s (sd: " << std << "s) on " << runs << " runs\n";

    cout << "[";
    for(int i=0; i<runs; i++) {
        cout << times[i].count() << ",";
    }
    cout << "]\n";
}

int main(int argc, char * argv[]){
    int size = 10000;
    int nthreads = 8;
    const int fs = 64;
    const int f0 = 1;
    const int f1 = 32;
    const int noct = 6;
    const int nvoi = 50;
    int fn = noct*nvoi;
    const int sigoutsize = size*fn*2;
    float c0 = 2*M_PI;
    float hz = 1;
    int runs = 5;
    if(argc > 0){
        size = stoi(argv[1]);
        fn = stoi(argv[2]);
    }

    double *sig1d = (double*)malloc(sizeof(double)*size);
    double *sig2d = (double*)malloc(sizeof(double)*size);
    double *sig3d = (double*)malloc(sizeof(double)*size);
    double *sigoutd = (double*)malloc(sizeof(double)*sigoutsize);

    for(int i=0; i<size; i++) {
        //Sig1: dynamic sine wave with varying frequency from 1Hz-7Hz, sampling rate of 64Hz.
        sig1d[i] = cos((2.0*3.1415*(hz+((double)(7*i)/size)))*((double)i/(float)fs));

        //Sig2: random numbers between 0-10.
        sig2d[i] = ((double)(rand() % 1000))/100.0;

        //Sig3: repeating non-smooth function.
        sig3d[i] = (i%10==0);
    }

    auto start = chrono::high_resolution_clock::now();
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed;

    //Initialize timing array
    chrono::duration<double> *times = (chrono::duration<double>*)malloc(sizeof(chrono::duration<double>)*runs);

    cout << "=========== BENCHMARKING RCWT ============" << endl;
    //RCWT sig1
    cout << "----- First test -----" << endl;
    cout << "Testing with sig1: dynamic sine wave (1Hz-7Hz)" << endl;
    cout << "Sample sig1: [";
    for(int n=0; n< 10; n++) {
        cout << sig1d[n] << ",";
    }
    cout << "...]" << endl;

    std::vector<double> sig1(sig1d, sig1d + size);
    std::vector<double> sig2(sig2d, sig2d + size);
    std::vector<double> sig3(sig3d, sig3d + size);
    std::vector<double> scale(fn);

    for(int i = 0; i < fn; i++){
        scale[i] = f0 + i * (double)(f1 - f0)/fn;
    }
    std::vector<std::vector<double>> sigout(fs);

    for(int k=0; k<runs; k++) {
        cout << ".";

        start = chrono::high_resolution_clock::now();

        sigout = cwt(sig1, scale,1000, nthreads);
        finish = chrono::high_resolution_clock::now();
        times[k] = finish - start;
        this_thread::sleep_for(chrono::microseconds(10000000));

    }
    cout << endl;
    cout << "RCWT on sig1 with length N: " << size;
    show_stats(times,runs);


    //FCWT sig2
    cout << "----- Second test -----" << endl;
    cout << "Testing with sig2: random floats between 0-10" << endl;
    cout << "Sample sig2: [";
    for(int n=0; n< 10; n++) {
        cout << sig2[n] << ",";
    }
    cout << "...]" << endl;

    for(int k=0; k<runs; k++) {
        cout << ".";
        start = chrono::high_resolution_clock::now();

        sigout = cwt(sig2, scale,1000, nthreads);

        finish = chrono::high_resolution_clock::now();
        times[k] = finish - start;

        this_thread::sleep_for(chrono::microseconds(10000000));
    }
    cout << endl;
    cout << "RCWT on sig2 with length N: " << size;
    show_stats(times,runs);


    //FCWT sig3
    cout << "----- Third test -----" << endl;
    cout << "Testing with sig3: Repeating non-smooth function (x%10==0)" << endl;
    cout << "Sample sig3: [";
    for(int n=0; n< 15; n++) {
        cout << sig3[n] << ",";
    }
    cout << "...]" << endl;

    for(int k=0; k<runs; k++) {
        cout << ".";
        start = chrono::high_resolution_clock::now();

        sigout = cwt(sig3, scale,1000, nthreads);

        finish = chrono::high_resolution_clock::now();
        times[k] = finish - start;

        this_thread::sleep_for(chrono::microseconds(10000000));
    }
    cout << endl;
    cout << "RCWT on sig3 with length N: " << size;
    show_stats(times,runs);

    cout << "=========== BENCHMARKING END ============" << endl;
    delete sig1d;
    delete sig2d;
    delete sig3d;
    delete sigoutd;
    delete times;

    return 0;
}