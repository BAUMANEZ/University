//
//  HandleFiles.cpp
//  ConvexHull
//
//  Created by Арсений Токарев on 17.12.2020.
//

#include <iostream>
#include <string>
#include <fstream>
#include "vector"
#include "Point.cpp"

using namespace std;

struct HandleFiles {
	static void fillVectorWithPointsFromTxt(const string &filePath,
											vector<Point> &points) {
		size_t numberOfPoints = 0;
		ifstream file(filePath);
		if (file.is_open()) {
			file >> numberOfPoints;
			points.reserve(numberOfPoints);

			double xCoord, yCoord;
			while (file >> xCoord >> yCoord) {
				points.push_back(Point(xCoord, yCoord));
			}
			file.close();
		}
		
		if (points.empty()) {
			cout << "Could not find the file or read data from it\n";
			exit(-1);
		}
	}
	
	static void writeToFile(const string filePath,
							const vector<Point> &convexHullPoints) {
	  ofstream saveFile;
	  saveFile.open(filePath);
	  if (saveFile.is_open()) {
		  saveFile << convexHullPoints.size() << endl;
		  for (int i = 0; i < convexHullPoints.size(); ++i) {
			  saveFile << convexHullPoints.at(i).x << " " << convexHullPoints.at(i).y << endl;
		  }
	  } else {
		  cout << "Could not write points to file\n";
		  exit(-1);
	  }
  }
};
