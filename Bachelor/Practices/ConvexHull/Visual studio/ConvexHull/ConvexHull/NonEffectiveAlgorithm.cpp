#include "DataModule.cpp"

struct NonEffectiveAlgorithm {
private:
	vector<Point> points, convexHullPoints;

	size_t giveElementWithGreatestXIn(const vector<Point>& points) {
		size_t indexOfElementWithMaxXCoord = 0;
		for (size_t i = 1; i < points.size() - 1; ++i) {
			if (points.at(i).x > points.at(indexOfElementWithMaxXCoord).x) {
				indexOfElementWithMaxXCoord = i;
			}
		}
		return indexOfElementWithMaxXCoord;
	}

	PointPosition determinePointPosition(const Point& lineStartPoint,
		const Point& lineEndPoint,
		Point pointToCheck) {
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
		for (size_t i = 1; i < points.size(); ++i) {
			const Point latterPoint = points.at(i); // A candidate to be the next point of the Convex Hull

			bool areOnTheSameSide = true;
			size_t indexInAllPoints = 0;
			PointPosition comparator = within;
			/*
					First iteration: find the first non "within the line" point.
					Second iteration: compare remaining points with this comparator.
			*/
			for (; indexInAllPoints < points.size(); ++indexInAllPoints) {
				comparator = determinePointPosition(formerPoint,
					latterPoint,
					points.at(indexInAllPoints));
				if (comparator != within) {
					++indexInAllPoints;
					break;
				}
			}

			for (; indexInAllPoints < points.size(); ++indexInAllPoints) {
				PointPosition currentPointPosition =
					determinePointPosition(formerPoint,
						latterPoint,
						points.at(indexInAllPoints));

				// If a point lays within a line we can skip it, because it conforms to our condition that all points should be on the side side from the line

				// If the point is on the other side from the line we break a cycle and go to the next candidate
				if ((currentPointPosition != comparator) && (currentPointPosition != within)) {
					areOnTheSameSide = false;
					break;
				}
			}

			// If all points are on the side side from the line, current candidate is a point of the Convex Hull. That's why we add it to our array of the CH Points
			if (areOnTheSameSide) {
				//if the element already forms a convexHull we skip the iteration
				if (find_if(convexHullPoints.begin(),
					convexHullPoints.end(),
					[latterPoint](const Point& pointOfConvexHull) {
						return latterPoint == pointOfConvexHull;
					}) != convexHullPoints.end()) {
					continue;
				}
				convexHullPoints.push_back(latterPoint);
				formerPoint = latterPoint;
				i = 0;
			}
		}
	}
	void run() {
		/*  FIND CONVEX HULL  */
		size_t indexOfElementWithMaxXCoord = giveElementWithGreatestXIn(points);
		swap(points.at(0), points.at(indexOfElementWithMaxXCoord));
		convexHullPoints.push_back(points.at(0));
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
