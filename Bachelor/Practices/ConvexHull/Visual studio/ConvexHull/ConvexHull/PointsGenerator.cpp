#ifndef PointsGeneratorFile
#define PointsGeneratorFile

#include "DataModule.cpp"

using namespace std;

struct PointsGenerator {

	template <typename T>
	static void checkIfInputIsCorrect(T& value,
		const string& message) {
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

			if (isRandomFilling) {
				for (int i = 0; i < numberOfPoints; ++i) {
					double coordX = ((double)(-1 + (rand() % 2) * 2)) * (double)rand() / RAND_MAX * 1000000;
					double coordY = ((double)(-1 + (rand() % 2) * 2)) * (double)rand() / RAND_MAX * 1000000;
					Point temporaryPoint = { coordX, coordY };
					while (find(coordinates.begin(), coordinates.end(), temporaryPoint) != coordinates.end()) {
						temporaryPoint.x = ((double)(-1 + (rand() % 2) * 2)) * (double)rand() / RAND_MAX * 1000000;
						temporaryPoint.y = ((double)(-1 + (rand() % 2) * 2)) * (double)rand() / RAND_MAX * 1000000;
					}
					coordinates.push_back(temporaryPoint);
					writeToTxt << temporaryPoint.x << " " << temporaryPoint.y << endl;
				}

				return;
			}

			for (int i = 0; i < numberOfPoints; ++i) {
				//Now I have no idea whether we have to make all these points unique. I am not going to handle
				//duplicates check in case a user chooses manual input of points.

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

