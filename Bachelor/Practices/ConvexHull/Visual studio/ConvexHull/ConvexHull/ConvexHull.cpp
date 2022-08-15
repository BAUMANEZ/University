#include "KirkpatrickSeidelAlgorithm.cpp"
#include "NonEffectiveAlgorithm.cpp"
#include "PointsGenerator.cpp"
#include <chrono>

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
	template<typename F, typename ...Args>
	static typename TimeT::rep execution(F func, Args&&... args)
	{
		auto start = std::chrono::system_clock::now();

		// Now call the function with all the parameters you need.
		func(std::forward<Args>(args)...);

		auto duration = std::chrono::duration_cast<TimeT>
			(std::chrono::system_clock::now() - start);

		return duration.count();
	}
};


int main(int argc, const char* argv[]) {
		//PointsGenerator::generatePoints();
	cout << "Non effective algorirthm: " << (double)measure<>::execution(NonEffectiveAlgorithm()) / 100 << endl;
	cout << "Kirkpatrick-Seidel algorithm: " << (double)measure<>::execution(KirkpatrickSeidelAlgorithm()) / 100 << endl;
	return 0;
}
