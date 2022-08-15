//
//  MyAlgorithm.cpp
//  Kirkpatrick-Seidel
//
//  Created by Арсений Токарев on 15.12.2020.
//
#include "../DataStructures/DataModule.cpp"
#include <tuple>
#include <optional>

using namespace std;

class KirkpatrickSeidelAlgorithm {
private:
	vector<Point> points;
	vector<Point> convexHullPoints;
	
	
	void sortXCoordinatesIn(vector<Point> &pointsToBeSorted) {
		sort(pointsToBeSorted.begin(),
			 pointsToBeSorted.end(),
			 [](const Point a, const Point b) {
					return a.x < b.x;
			 }
		);
	}
	
	double findMedianIn(vector<double> values,
					    const size_t numberOfElements) {
	  
		if (numberOfElements % 2 == 0) {
	  
			nth_element(values.begin(),
						values.begin() + numberOfElements / 2,
						values.end());
	  
			nth_element(values.begin(),
						values.begin() + (numberOfElements - 1) / 2,
						values.end());
	  
			return (double)(values[(numberOfElements - 1) / 2]
							+ values[numberOfElements / 2])
				   / 2.0;
		}
		else {
	
			nth_element(values.begin(),
						values.begin() + numberOfElements / 2,
						values.end());

			return (double)values[numberOfElements / 2];
		}
	}
	vector<PointPair> makePointPairs(const vector<Point> &setOfPoints,
									 vector<Point> &candidates) {
		
		vector<PointPair> pointPairs {};
		size_t startIndex = 0;
		if (setOfPoints.size() % 2 != 0) {
			candidates.emplace_back(*(setOfPoints.begin()));
			startIndex = 1;
		}
		for (auto iterator = setOfPoints.begin() + startIndex;
			 iterator < setOfPoints.end();
			 iterator += 2) {
			if (iterator->x < (iterator + 1)->x) {
				pointPairs.emplace_back(PointPair(*iterator,
												  *(iterator + 1)));
			} else {
				pointPairs.emplace_back(PointPair(*(iterator + 1),
												  *(iterator)));
			}
		}
		
		return pointPairs;
	}
	
	void getSlopeCharacterPairs(vector<PointPair> &smallPairs,
								vector<PointPair> &equalPairs,
								vector<PointPair> &largePairs,
								const vector<double> &slopes,
								const vector<PointPair> &pointPairs,
								const double medianSlope) {
		smallPairs.reserve(pointPairs.size() / 3);
		equalPairs.reserve(pointPairs.size() / 3);
		largePairs.reserve(pointPairs.size() / 3);
		for (int i = 0; i < pointPairs.size(); ++i) {
			const double slope = slopes[i];
			const PointPair pointPair = pointPairs[i];
			if (slope < medianSlope) {
				smallPairs.emplace_back(pointPair);
			}
			else if (slope > medianSlope) {
				largePairs.emplace_back(pointPair);
			}
			else if (abs(slope - medianSlope) <= EPS){
				equalPairs.emplace_back(pointPair);
			}
		}
	}
	
	vector<double> calculateSlopes(vector<PointPair> &pointPairs,
								   vector<Point> &candidates) {
		vector<double> slopes {};
		for (auto iterator = pointPairs.begin();
			 iterator < pointPairs.end();
			 ++iterator) {
			const auto pair = *iterator;
			if (pair.first.x == pair.second.x) {
				pointPairs.erase(iterator--);
				candidates.emplace_back(pair.first.y > pair.second.y ? pair.first : pair.second);
			} else {
				const double slope = (pair.first.y - pair.second.y) / (pair.first.x - pair.second.x);
				slopes.emplace_back(slope);
			}
		}
		return slopes;
	}
	
	optional<PointPair> continueForUpperBridge(const double medianX,
											   const double medianSlope,
											   const vector<PointPair> &largePairs,
											   const vector<PointPair> &equalPairs,
											   const vector<PointPair> &smallPairs,
											   const vector<Point> &remainingPoints,
											   vector<Point> &candidates) {
		
		double height = (remainingPoints[0].y - medianSlope * remainingPoints[0].x);
		vector<Point> pointsOnTheMaximizedLine;
		for(auto &currentPoint: remainingPoints) {
			const double currentHeight = (currentPoint.y - medianSlope * currentPoint.x);
			if (currentHeight > height) {
				pointsOnTheMaximizedLine.clear();
				pointsOnTheMaximizedLine.emplace_back(currentPoint);
				height = currentHeight;
			}
			else if (abs(currentHeight - height) <= EPS) {
				pointsOnTheMaximizedLine.emplace_back(currentPoint);
			}
		 }
		sortXCoordinatesIn(pointsOnTheMaximizedLine);
		
		Point maximizingMin = *(pointsOnTheMaximizedLine.begin()),
			  maximizingMax = *(pointsOnTheMaximizedLine.end() - 1);
		
		if (maximizingMin.x <= medianX && maximizingMax.x > medianX) {
			return PointPair(maximizingMin,
							 maximizingMax);
		}
		
		if (maximizingMax.x <= medianX) {
			for (auto &pointPair : largePairs)
				candidates.emplace_back(pointPair.second);
			for (auto &pointPair : equalPairs)
				candidates.emplace_back(pointPair.second);
			
			for (auto &pointPair : smallPairs) {
				candidates.emplace_back(pointPair.first);
				candidates.emplace_back(pointPair.second);
			}
		}
		
		if (maximizingMin.x > medianX) {
			for (auto &pointPair : smallPairs)
				candidates.emplace_back(pointPair.first);
			for (auto &pointPair : equalPairs)
				candidates.emplace_back(pointPair.first);
			
			for (auto &pointPair : largePairs) {
				candidates.emplace_back(pointPair.first);
				candidates.emplace_back(pointPair.second);
			}
		}
		return {};
	}
	
	optional<PointPair> continueForLowerBridge(const double medianX,
											   const double medianSlope,
											   const vector<PointPair> &largePairs,
											   const vector<PointPair> &equalPairs,
											   const vector<PointPair> &smallPairs,
											   const vector<Point> &remainingPoints,
											   vector<Point> &candidates) {
		double height = (remainingPoints[0].y - medianSlope * remainingPoints[0].x);
		vector<Point> pointsOnTheMinimizedLine;
		for(auto &currentPoint: remainingPoints) {
			const double currentHeight = (currentPoint.y - medianSlope * currentPoint.x);
			if (currentHeight < height) {
				pointsOnTheMinimizedLine.clear();
				pointsOnTheMinimizedLine.emplace_back(currentPoint);
				height = currentHeight;
			}
			else if (abs(currentHeight - height) <= EPS) {
				pointsOnTheMinimizedLine.emplace_back(currentPoint);
			}
		 }
		
		sortXCoordinatesIn(pointsOnTheMinimizedLine);
		
		Point minimizedMin = *(pointsOnTheMinimizedLine.begin()),
			  minimizedMax = *(pointsOnTheMinimizedLine.end() - 1);
		
		if (minimizedMin.x <= medianX && minimizedMax.x > medianX) {
			return PointPair(minimizedMin,
							 minimizedMax);
		}
		
		if (minimizedMax.x <= medianX) {
			for (auto &pointPair : smallPairs)
				candidates.emplace_back(pointPair.second);
			for (auto &pointPair : equalPairs)
				candidates.emplace_back(pointPair.second);
			
			for (auto &pointPair : largePairs) {
				candidates.emplace_back(pointPair.first);
				candidates.emplace_back(pointPair.second);
			}
		}
		
		if (minimizedMin.x > medianX) {
			for (auto &pointPair : largePairs)
				candidates.emplace_back(pointPair.first);
			for (auto &pointPair : equalPairs)
				candidates.emplace_back(pointPair.first);
			
			for (auto &pointPair : smallPairs) {
				candidates.emplace_back(pointPair.first);
				candidates.emplace_back(pointPair.second);
			}
		}
		return {};
		
	}
	
	PointPair getBridge(vector<Point> &remainingPoints,
						double medianX,
						bool isUpperBridge) {
		vector<Point> candidates {};
		const size_t numberOfRemainingPoints = remainingPoints.size();
		if (numberOfRemainingPoints == 2) {
			return PointPair(*remainingPoints.begin(),
							 *(remainingPoints.begin() + 1));
		}
		
		vector<PointPair> pointPairs = makePointPairs(remainingPoints,
													  candidates);
		
		auto slopes = calculateSlopes(pointPairs,
									  candidates);
		
		const double medianSlope = findMedianIn(slopes, slopes.size());
		
		vector<PointPair> smallPairs {}, equalPairs {}, largePairs {};
		getSlopeCharacterPairs(smallPairs,
							   equalPairs,
							   largePairs,
							   slopes,
							   pointPairs,
							   medianSlope);
		const auto result = isUpperBridge ? continueForUpperBridge(medianX,
																	medianSlope,
																	largePairs,
																	equalPairs,
																	smallPairs,
																    remainingPoints,
																	candidates)
		 :
											  continueForLowerBridge(medianX,
																	 medianSlope,
																	 largePairs,
																	 equalPairs,
																	 smallPairs,
																	 remainingPoints,
																	 candidates);
		if (const auto bridge = result) {
			return bridge.value();
		}
		return getBridge(candidates, medianX, isUpperBridge);
	}
	
	vector<double> getXCoordinatesIn(const vector<Point> &points) {
		vector<double> result;
		result.reserve(points.size());
		for (auto &point : points)
			result.emplace_back(point.x);
		return result;
	}
	
	void connect(bool isUpperHull,
				 const Point &firstPoint,
				 const Point &secondPoint,
				 vector<Point> &setOfPoints) {
		
		const double medianX = findMedianIn(getXCoordinatesIn(setOfPoints), setOfPoints.size());
		auto bridge = getBridge(setOfPoints, medianX, isUpperHull);
		if (bridge.first.x > bridge.second.x) {
			const Point tempPoint = bridge.first;
			bridge.first = bridge.second;
			bridge.second = tempPoint;
		}
		
		vector<Point> leftSetOfPoints ;
		leftSetOfPoints.reserve(setOfPoints.size() / 2);
		leftSetOfPoints.emplace_back(bridge.first);
		for (auto &point : setOfPoints) {
			if (point.x < bridge.first.x)
				leftSetOfPoints.emplace_back(point);
		}

		vector<Point> rightSetOfPoints;
		rightSetOfPoints.reserve(setOfPoints.size() / 2);
		rightSetOfPoints.emplace_back(bridge.second);
		for (auto &point : setOfPoints) {
			if (point.x > bridge.second.x)
				rightSetOfPoints.emplace_back(point);
		}
		
		if (bridge.first == firstPoint)
			convexHullPoints.emplace_back(bridge.first);
		else
			connect(isUpperHull, firstPoint, bridge.first, leftSetOfPoints);

		if (bridge.second == secondPoint)
			convexHullPoints.emplace_back(bridge.second);
		else
			connect(isUpperHull, bridge.second, secondPoint, rightSetOfPoints);
	}
	
	void getPartOfHull(bool isUpperHull,
					   const Point minPoint,
					   const Point maxPoint) {
		if (minPoint == maxPoint) {
			convexHullPoints.emplace_back(minPoint);
			return;
		}
		vector<Point> hullPoints ;
		hullPoints.reserve(points.size());
		hullPoints.emplace_back(minPoint);
		for (auto &point : points ) {
			if (point.x > minPoint.x && point.x < maxPoint.x)
				hullPoints.emplace_back(point);
		}
		hullPoints.emplace_back(maxPoint);
		
		connect(isUpperHull, minPoint, maxPoint, hullPoints);
	}
	
	tuple<
		Point,
		Point,
		Point,
		Point
	> getMaximumAndMiniumPointsForBothHulls() {
		Point maxUpperPoint,
			  minUpperPoint,
			  maxLowerPoint,
			  minLowerPoint;
		
		maxUpperPoint = maxLowerPoint = Point(INT_MIN, INT_MIN);
		minUpperPoint = minLowerPoint = Point(INT_MAX, INT_MAX);
		
		for (auto &point : points) {
			if (abs(point.x - minUpperPoint.x) <= EPS) {
				if (point.y > minUpperPoint.y)
					minUpperPoint = point;
				else
					minLowerPoint = point;
			}
			else if (point.x <= minUpperPoint.x)
				minUpperPoint = minLowerPoint = point;
			
			if (abs(point.x - maxUpperPoint.x) <= EPS) {
				if (point.y > maxUpperPoint.y)
					maxUpperPoint = point;
				else
					maxLowerPoint = point;
			}
			else if (point.x >= maxUpperPoint.x)
				maxUpperPoint = maxLowerPoint = point;
		}
		return make_tuple(maxUpperPoint, minUpperPoint, maxLowerPoint, minLowerPoint);
	}
	
	void saveToTxt() {
		HandleFiles::writeToFile(FILE_SAVE_PATH_KIRKPATRICK_SEIDEL,
								 convexHullPoints);
	}
	
	void connectHullPoints() {
		
	}
	
	void run() {
		const auto [ maxUpperPoint,
					 minUpperPoint,
					 maxLowerPoint,
					 minLowerPoint ] = getMaximumAndMiniumPointsForBothHulls();
		getPartOfHull(true, minUpperPoint, maxUpperPoint);
		if (minUpperPoint == minLowerPoint)
			convexHullPoints.erase(convexHullPoints.begin());
		reverse(convexHullPoints.begin(), convexHullPoints.end());
		getPartOfHull(false, minLowerPoint, maxLowerPoint);
		if (!(maxLowerPoint == maxUpperPoint))
			convexHullPoints.emplace_back(maxUpperPoint);
		else
			convexHullPoints.erase(convexHullPoints.end() - 1);
		saveToTxt();
	}
	
public:
	void operator()() {
		this->run();
	}
	
	KirkpatrickSeidelAlgorithm(const string &filePath = FILE_COORDINATES_PATH) {
		HandleFiles::fillVectorWithPointsFromTxt(filePath,
												 points);
		convexHullPoints = {};
	}

};
