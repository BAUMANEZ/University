#ifndef DataStructures
#define DataStructures

#define FILE_COORDINATES_PATH "\\Mac\Home\Desktop\ConvexHull_BMSTU\Visual studio\ConvexHull\ConvexHull\Text files\coordinates.txt"

#define FILE_SAVE_PATH_BRUTE_FORCE "\\Mac\Home\Desktop\ConvexHull_BMSTU\Visual studio\ConvexHull\ConvexHull\Text files\ConvexHullPoints(MethodI).txt"

#define FILE_SAVE_PATH_KIRKPATRICK_SEIDEL "\\Mac\Home\Desktop\ConvexHull_BMSTU\Visual studio\ConvexHull\ConvexHull\Text files\ConvexHullPoints(MethodII).txt"

#include <iostream>
#include <vector>
#include <algorithm>
#include "PointPair.cpp"
#include "HandleFiles.cpp"
using namespace std;

enum Result { isEmpty, isFilled };
enum PointPosition { oneSide, otherSide, within };

#endif
