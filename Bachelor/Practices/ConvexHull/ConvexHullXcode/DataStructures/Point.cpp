//
//  Point.cpp
//  ConvexHull
//
//  Created by Арсений Токарев on 17.12.2020.
//

#ifndef PointFile
#define PointFile
#define EPS 10e-10

#include <iostream>

struct Point {
	double x;
	double y;
	
	bool operator== (const Point& coordinate) const {
		const bool xAreEqual = abs(coordinate.x - this->x) <= EPS;
		const bool yAreEqual = abs(coordinate.y - this->y) <= EPS;
		return (xAreEqual && yAreEqual);
	}
	
	Point(double coordX = 0.0, double coordY = 0.0) {
		x = coordX;
		y = coordY;
	}
};

#endif
