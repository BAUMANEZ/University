//
//  PointPair.cpp
//  ConvexHull
//
//  Created by Арсений Токарев on 17.12.2020.
//
#ifndef PointPairFile
#define PointPairFile


#include "Point.cpp"

struct PointPair {
public:
	Point first;
	Point second;
	
	PointPair(const Point firstPoint, const Point secondPoint) {
		first = firstPoint;
		second = secondPoint;
	}
};

#endif
