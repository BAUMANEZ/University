#include "DataModule.cpp"
#include <tuple>
#include <optional>

using namespace std;

class KirkpatrickSeidelAlgorithm {
private:
	vector<Point> points;
	vector<Point> convexHullPoints;


	void sortXCoordinatesIn(vector<Point>& pointsToBeSorted) {
		sort(pointsToBeSorted.begin(),
			pointsToBeSorted.end(),
			[](const Point a, const Point b) {
				return a.x < b.x;
			}
		);
	}

	template<
		typename T,
		typename = typename enable_if<is_arithmetic<T>::value, T>::type
	>
		double getMedianIn(vector<T> values) {
		sort(values.begin(), values.end());
		const size_t numberOfValues = values.size();
		switch (numberOfValues % 2) {
		case 0:
			return (values[numberOfValues / 2] + values[numberOfValues / 2 - 1]) / 2;
		default:
			return values.at(numberOfValues / 2);
		}
	}

	vector<PointPair> makePointPairs(const vector<Point>& setOfPoints,
		vector<Point>& candidates) {

		vector<PointPair> pointPairs{};
		size_t startIndex = 0;
		if (setOfPoints.size() % 2 != 0) {
			candidates.push_back(*(setOfPoints.begin()));
			startIndex = 1;
		}
		for (auto iterator = setOfPoints.begin() + startIndex;
			iterator < setOfPoints.end();
			iterator += 2) {
			pointPairs.push_back(PointPair(*iterator,
				*(iterator + 1)));
		}

		return pointPairs;
	}

	tuple<
		vector<PointPair>,
		vector<PointPair>,
		vector<PointPair>
	> getSlopeCharacterPairs(const vector<double>& slopes,
		const vector<PointPair>& pointPairs,
		const double medianSlope) {
		if (slopes.size() != pointPairs.size()) {
			cout << "Slopes.size != pointPairs.size" << endl;
			exit(-4);
		}
		vector<PointPair> smallPairs{}, equalPairs{}, largePairs{};
		smallPairs.reserve(pointPairs.size() / 3);
		equalPairs.reserve(pointPairs.size() / 3);
		largePairs.reserve(pointPairs.size() / 3);
		for (int i = 0; i < pointPairs.size(); ++i) {
			const double slope = slopes.at(i);
			const PointPair pointPair = pointPairs.at(i);
			if (slope < medianSlope) {
				smallPairs.push_back(pointPair);
			}
			else if (slope > medianSlope) {
				largePairs.push_back(pointPair);
			}
			else if (abs(slope - medianSlope) <= EPS) {
				equalPairs.push_back(pointPair);
			}
		}
		return make_tuple(smallPairs, equalPairs, largePairs);
	}

	vector<double> calculateSlopes(vector<PointPair>& pointPairs,
		vector<Point>& candidates) {
		vector<double> slopes{};
		for (auto iterator = pointPairs.begin();
			iterator < pointPairs.end();
			++iterator) {
			const auto pair = *iterator;
			if (pair.first.x == pair.second.x) {
				pointPairs.erase(iterator--);
				candidates.push_back(pair.first.y > pair.second.y ? pair.first : pair.second);
			}
			else {
				const double slope = (pair.first.y - pair.second.y) / (pair.first.x - pair.second.x);
				slopes.push_back(slope);
			}
		}
		return slopes;
	}

	optional<PointPair> continueForUpperBridge(const double medianX,
		const double medianSlope,
		const vector<PointPair>& largePairs,
		const vector<PointPair>& equalPairs,
		const vector<PointPair>& smallPairs,
		const vector<Point>& remainingPoints,
		vector<Point>& candidates) {

		double height = (remainingPoints.at(0).y - medianSlope * remainingPoints.at(0).x);
		vector<Point> pointsOnTheMaximizedLine;
		for (auto& currentPoint : remainingPoints) {
			const double currentHeight = (currentPoint.y - medianSlope * currentPoint.x);
			if (currentHeight > height) {
				pointsOnTheMaximizedLine.clear();
				pointsOnTheMaximizedLine.push_back(currentPoint);
				height = currentHeight;
			}
			else if (abs(currentHeight - height) <= EPS) {
				pointsOnTheMaximizedLine.push_back(currentPoint);
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
			vector<PointPair> largeAndEqualPairs;
			largeAndEqualPairs.reserve(largePairs.size() + equalPairs.size());
			largeAndEqualPairs.insert(largeAndEqualPairs.end(),
				largePairs.begin(),
				largePairs.end());
			largeAndEqualPairs.insert(largeAndEqualPairs.end(),
				equalPairs.begin(),
				equalPairs.end());

			for (auto& pointPair : largeAndEqualPairs)
				candidates.push_back(pointPair.second);

			for (auto& pointPair : smallPairs) {
				candidates.push_back(pointPair.first);
				candidates.push_back(pointPair.second);
			}
		}

		if (maximizingMin.x > medianX) {
			vector<PointPair> smallAndEqualPairs;
			smallAndEqualPairs.reserve(smallPairs.size() + equalPairs.size());
			smallAndEqualPairs.insert(smallAndEqualPairs.end(),
				smallPairs.begin(),
				smallPairs.end());
			smallAndEqualPairs.insert(smallAndEqualPairs.end(),
				equalPairs.begin(),
				equalPairs.end());

			for (auto& pointPair : smallAndEqualPairs)
				candidates.push_back(pointPair.first);

			for (auto& pointPair : largePairs) {
				candidates.push_back(pointPair.first);
				candidates.push_back(pointPair.second);
			}
		}
		return {};
	}

	optional<PointPair> continueForLowerBridge(const double medianX,
		const double medianSlope,
		const vector<PointPair>& largePairs,
		const vector<PointPair>& equalPairs,
		const vector<PointPair>& smallPairs,
		const vector<Point>& remainingPoints,
		vector<Point>& candidates) {
		double height = (remainingPoints.at(0).y - medianSlope * remainingPoints.at(0).x);
		vector<Point> pointsOnTheMinimizedLine;
		for (auto& currentPoint : remainingPoints) {
			const double currentHeight = (currentPoint.y - medianSlope * currentPoint.x);
			if (currentHeight < height) {
				pointsOnTheMinimizedLine.clear();
				pointsOnTheMinimizedLine.push_back(currentPoint);
				height = currentHeight;
			}
			else if (abs(currentHeight - height) <= EPS) {
				pointsOnTheMinimizedLine.push_back(currentPoint);
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
			vector<PointPair> smallAndEqualPairs;
			smallAndEqualPairs.reserve(smallPairs.size() + equalPairs.size());
			smallAndEqualPairs.insert(smallAndEqualPairs.end(),
				smallPairs.begin(),
				smallPairs.end());
			smallAndEqualPairs.insert(smallAndEqualPairs.end(),
				equalPairs.begin(),
				equalPairs.end());

			for (auto& pointPair : smallAndEqualPairs)
				candidates.push_back(pointPair.second);

			for (auto& pointPair : largePairs) {
				candidates.push_back(pointPair.first);
				candidates.push_back(pointPair.second);
			}
		}

		if (minimizedMin.x > medianX) {
			vector<PointPair> largeAndEqualPairs;
			largeAndEqualPairs.reserve(largePairs.size() + equalPairs.size());
			largeAndEqualPairs.insert(largeAndEqualPairs.end(),
				largePairs.begin(),
				largePairs.end());
			largeAndEqualPairs.insert(largeAndEqualPairs.end(),
				equalPairs.begin(),
				equalPairs.end());

			for (auto& pointPair : largeAndEqualPairs)
				candidates.push_back(pointPair.first);

			for (auto& pointPair : smallPairs) {
				candidates.push_back(pointPair.first);
				candidates.push_back(pointPair.second);
			}
		}
		return {};

	}

	PointPair getBridge(vector<Point>& remainingPoints,
		double medianX,
		bool isUpperBridge) {
		vector<Point> candidates{};
		const size_t numberOfRemainingPoints = remainingPoints.size();
		if (numberOfRemainingPoints == 2) {
			return PointPair(remainingPoints.at(0),
				remainingPoints.at(1));
		}

		vector<PointPair> pointPairs = makePointPairs(remainingPoints,
			candidates);

		const auto slopes = calculateSlopes(pointPairs,
			candidates);

		const double medianSlope = getMedianIn(slopes);

		// C++17 and higher!!!
		const auto [smallPairs, equalPairs, largePairs] = getSlopeCharacterPairs(slopes,
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
		sortXCoordinatesIn(candidates);
		return getBridge(candidates, medianX, isUpperBridge);
	}

	vector<double> getXCoordinatesIn(const vector<Point>& points) {
		vector<double> result;
		result.reserve(points.size());
		for (auto& point : points) result.push_back(point.x);
		return result;
	}

	void connect(bool isUpperHull,
		const Point& firstPoint,
		const Point& secondPoint,
		vector<Point>& setOfPoints) {
		sortXCoordinatesIn(setOfPoints);
		const double medianX = getMedianIn(getXCoordinatesIn(setOfPoints));
		auto bridge = getBridge(setOfPoints, medianX, isUpperHull);
		if (bridge.first.x > bridge.second.x) {
			const Point tempPoint = bridge.first;
			bridge.first = bridge.second;
			bridge.second = tempPoint;
		}

		vector<Point> leftSetOfPoints;// { bridge.first };
		leftSetOfPoints.reserve(setOfPoints.size() / 2);
		leftSetOfPoints.push_back(bridge.first);
		for (auto& point : setOfPoints) {
			if (point.x < bridge.first.x)
				leftSetOfPoints.push_back(point);
		}

		vector<Point> rightSetOfPoints; //{ bridge.second };
		rightSetOfPoints.reserve(setOfPoints.size() / 2);
		rightSetOfPoints.push_back(bridge.second);
		for (auto& point : setOfPoints) {
			if (point.x > bridge.second.x)
				rightSetOfPoints.push_back(point);
		}

		if (bridge.first == firstPoint)
			convexHullPoints.push_back(bridge.first);
		else
			connect(isUpperHull, firstPoint, bridge.first, leftSetOfPoints);

		if (bridge.second == secondPoint)
			convexHullPoints.push_back(bridge.second);
		else
			connect(isUpperHull, bridge.second, secondPoint, rightSetOfPoints);
	}

	void getPartOfHull(bool isUpperHull,
		const Point minPoint,
		const Point maxPoint) {
		if (minPoint == maxPoint) {
			convexHullPoints.push_back(minPoint);
			return;
		}
		vector<Point> hullPoints;//{ minPoint };
		hullPoints.reserve(points.size());
		hullPoints.push_back(minPoint);
		for (auto& point : points) {
			if (point.x > minPoint.x && point.x < maxPoint.x)
				hullPoints.push_back(point);
		}
		hullPoints.push_back(maxPoint);

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

		for (auto& point : points) {
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
		const auto [maxUpperPoint,
			minUpperPoint,
			maxLowerPoint,
			minLowerPoint] = getMaximumAndMiniumPointsForBothHulls();
		getPartOfHull(true, minUpperPoint, maxUpperPoint);
		if (minUpperPoint == minLowerPoint)
			convexHullPoints.erase(convexHullPoints.begin());
		reverse(convexHullPoints.begin(), convexHullPoints.end());
		getPartOfHull(false, minLowerPoint, maxLowerPoint);
		if (!(maxLowerPoint == maxUpperPoint))
			convexHullPoints.push_back(maxUpperPoint);
		saveToTxt();
	}

public:
	void operator()() {
		this->run();
	}

	KirkpatrickSeidelAlgorithm(const string& filePath = FILE_COORDINATES_PATH) {
		HandleFiles::fillVectorWithPointsFromTxt(filePath,
			points);
		convexHullPoints = {};
	}

};
