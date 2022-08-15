//
//  main.cpp
//  Kirkpatrick-Seidel
//
//  Created by Арсений Токарев on 14.12.2020.
//
#include "../Kirkpatrick-Seidel/Kirkpatrick-Seidel Algorithm.cpp"
#include "../NonEffectiveAlgorithm/NonEffectiveAlgorithm.cpp"
#include "../Helpers/PointsGenerator.cpp"
#include <chrono>

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
	template<typename F, typename ...Args>
	static typename TimeT::rep execution(F func, Args&&... args)
	{
		auto start = std::chrono::high_resolution_clock::now();
 
		// Now call the function with all the parameters you need.
		func(std::forward<Args>(args)...);
 
		auto duration = std::chrono::duration_cast< TimeT>
							(std::chrono::high_resolution_clock::now() - start);
 
		return duration.count();
	}
};


int main(int argc, const char *argv[]) {
	PointsGenerator::generatePoints();

	cout << "Kirkpatrick-Seidel algorithm: " << (double) measure<>::execution(KirkpatrickSeidelAlgorithm()) / 1000 << endl;

	cout << "Non effective algorirthm: " << (double)measure<>::execution(NonEffectiveAlgorithm()) / 1000 << endl;
	return 0;
}
