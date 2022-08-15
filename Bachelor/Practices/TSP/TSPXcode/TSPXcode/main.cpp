//
//  main.cpp
//  Travelling Salesman
//
//  Created by Арсений Токарев on 07.05.2021.
//
#ifndef RoutesGenerator
#define RoutesGenerator

#include <iostream>
#include <fstream>
#include "Helpers/Algorithms.cpp"

int main(int argc, const char * argv[]) {
    NNAlgorithm().run();
    TwoOptAlgorithm().run();
    ThreeOptAlgorithm().run();
    return 0;
}

#endif
