//
//  NonEffectiveAlgorithm.cpp
//  ConvexHull
//
//  Created by Арсений Токарев on 17.12.2020.
//

#include "../DataStructures/DataModule.cpp"

struct NonEffectiveAlgorithm {
private:
	vector<Point> points, convexHullPoints;
	
	const auto giveIteratorOfPointWithMaxX() {
		auto iteratorOfElementWithMaxX = points.begin();
		for (auto iterator = points.begin() + 1;
			 iterator < points.end();
			 ++iterator) {
			if (iterator->x > iteratorOfElementWithMaxX->x)
				iteratorOfElementWithMaxX = iterator;
		}
		return iteratorOfElementWithMaxX;
	}

	PointPosition determinePointPosition(const Point &lineStartPoint,
										 const Point &lineEndPoint,
										 const Point &pointToCheck) {
		const double result = (pointToCheck.x - lineStartPoint.x) * (lineEndPoint.y - lineStartPoint.y) - (pointToCheck.y - lineStartPoint.y) * (lineEndPoint.x - lineStartPoint.x);

		PointPosition positionFromTheLine;
		if (result < 0) {
			positionFromTheLine = oneSide;
		}
		else if (result > 0) {
			positionFromTheLine = otherSide;
		}
		else {
			positionFromTheLine = within;
		}

		return positionFromTheLine;
	}

	void fill() {
		if (points.empty() || convexHullPoints.empty()) {
			cout << "There are no points to form a convex hull!\n";
			exit(-2);
		}
		
		Point formerPoint = convexHullPoints.back();
		for (auto iterator = points.begin() + 1;
			 iterator < points.end();
			 ++iterator) {
			const Point latterPoint = *iterator; // A candidate to be the next point of the Convex Hull

			bool areOnTheSameSide = true;
			auto nestedIterator = points.begin();
			PointPosition comparator = within;
			
			/*
					First iteration: find the first non "within the line" point.
					Second iteration: compare remaining points with this comparator.
			*/
			for (; nestedIterator < points.end(); ++nestedIterator) {
				comparator = determinePointPosition(formerPoint,
													latterPoint,
													*nestedIterator);
				if (comparator != within) {
					++nestedIterator;
					break;
				}
			}
			
			for (; nestedIterator < points.end(); ++nestedIterator) {
				PointPosition currentPointPosition =
					determinePointPosition(formerPoint,
										   latterPoint,
										   *nestedIterator);

				if ((currentPointPosition != comparator) && (currentPointPosition != within)) {
					areOnTheSameSide = false;
					break;
				}
			}
			
			if (areOnTheSameSide) {
				if (find_if(convexHullPoints.begin(),
							convexHullPoints.end(),
							[latterPoint](const Point& pointOfConvexHull) {
								return latterPoint == pointOfConvexHull;
							}) != convexHullPoints.end()) {
					continue;
				}
				convexHullPoints.emplace_back(latterPoint);
				formerPoint = latterPoint;
				iterator = points.begin();
			}
		}
	}
	void run() {
		/*  FIND CONVEX HULL  */
		auto iteratorOfElementWithMaxX = giveIteratorOfPointWithMaxX();
		iter_swap(points.begin(), iteratorOfElementWithMaxX);
		convexHullPoints.emplace_back(points.at(0));
		fill();
		HandleFiles::writeToFile(FILE_SAVE_PATH_BRUTE_FORCE,
								 convexHullPoints);
		
	}
public:
	void operator()() {
		this->run();
	}
	
	NonEffectiveAlgorithm() {
		HandleFiles::fillVectorWithPointsFromTxt(FILE_COORDINATES_PATH,
												 points);
		convexHullPoints = {};
	}
};
