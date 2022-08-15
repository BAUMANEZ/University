//
//  DataModule.cpp
//  ConvexHull
//
//  Created by Арсений Токарев on 17.12.2020.
//

#ifndef DataStructures
#define DataStructures

#define FILE_COORDINATES_PATH "/Users/arsenytokarev/Desktop/ConvexHull_BMSTU/ConvexHull/Helpers/Text files/coordinates.txt"

#define FILE_SAVE_PATH_BRUTE_FORCE "/Users/arsenytokarev/Desktop/ConvexHull_BMSTU/ConvexHull/Helpers/Text files/ConvexHullPoints(MethodI).txt"

#define FILE_SAVE_PATH_KIRKPATRICK_SEIDEL "/Users/arsenytokarev/Desktop/ConvexHull_BMSTU/ConvexHull/Helpers/Text files/ConvexHullPoints(MethodII).txt"

#include <iostream>
#include <vector>
#include "PointPair.cpp"
#include "HandleFiles.cpp"
using namespace std;

enum Result { isEmpty, isFilled };
enum PointPosition { oneSide, otherSide, within };

#endif
