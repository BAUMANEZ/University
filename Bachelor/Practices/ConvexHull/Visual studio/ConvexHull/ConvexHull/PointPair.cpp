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
