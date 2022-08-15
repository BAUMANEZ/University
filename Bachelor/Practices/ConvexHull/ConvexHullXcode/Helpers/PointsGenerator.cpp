//
//  PointsGenerator.cpp
//  ConvexHull
//
//  Created by Арсений Токарев on 17.12.2020.
//
#ifndef PointsGeneratorFile
#define PointsGeneratorFile

#include "../DataStructures/DataModule.cpp"

using namespace std;

struct PointsGenerator {
	
	template <typename T>
	static void checkIfInputIsCorrect(T &value,
									  const string &message) {
		while (cout << message << "\n>> " && !(cin >> value)) {
			cin.clear();
			cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			cout << "Invalid input; press ENTER to try again.\n";
			cin.get();
		}
	}

	static void generatePoints() {
		ofstream writeToTxt;
		writeToTxt.open(FILE_COORDINATES_PATH);

		if (writeToTxt.is_open()) {
			bool isRandomFilling;
			checkIfInputIsCorrect(isRandomFilling, "Random(1) or Manual(0) filling?");

			vector<Point> coordinates;

			int numberOfPoints = 0;
			checkIfInputIsCorrect(numberOfPoints, "How many points do you need?");
			writeToTxt << numberOfPoints << endl;
			
			double xStart, xEnd, yStart, yEnd;
			checkIfInputIsCorrect(xStart, "Min X: ");
			checkIfInputIsCorrect(xEnd, "Max X");
			checkIfInputIsCorrect(yStart, "Min Y: ");
			checkIfInputIsCorrect(yEnd, "Max Y");
			
			if (isRandomFilling) {
				for (int i = 0; i < numberOfPoints; ++i) {
					double coordX = xStart +(double)rand() / RAND_MAX * (xEnd - xStart);
					double coordY = yStart + (double)rand() / RAND_MAX * (yEnd - yStart);
					Point temporaryPoint = { coordX, coordY };
					coordinates.push_back(temporaryPoint);
					writeToTxt << temporaryPoint.x << " " << temporaryPoint.y << endl;
				}

				return;
			}

			for (int i = 0; i < numberOfPoints; ++i) {
				cout << "POINT" << i + 1 << endl;
				double coordX = 0.0, coordY = 0.0;
				checkIfInputIsCorrect(coordX, "Coordinate X");
				checkIfInputIsCorrect(coordY, "Coordinate Y");
				const Point temporaryPoint = { coordX, coordY };
				coordinates.push_back(temporaryPoint);
				writeToTxt << temporaryPoint.x << " " << temporaryPoint.y << endl;
				cout << "\n------------------\n\n";
			}
			return;
		}
	}
};

#endif

