//
//  2Opt.cpp
//  Travelling Salesman
//
//  Created by Арсений Токарев on 08.05.2021.
//

#define INPUT_FILE_PATH "/Users/arsenytokarev/Desktop/BMSTU_Practices/TSP/TSPXcode/TSPXcode/Helpers/costMatrix.txt"

#include "Tour.cpp"
#include <float.h>

class TSPAlgorithm {
protected:
    virtual Tour giveBestTour(const Tour& initTour) { return {}; }
    
public:
    virtual void run() {
        auto tour = giveBestTour(Tour(INPUT_FILE_PATH));
        for (size_t i = 0; i < tour.count(); ++i) {
            std::cout << tour[i] + 1 << " -> ";
        }
        std::cout << "1";
        std::cout << std::endl;
        std::cout << tour.tripLength() << std::endl << std::endl;
    }
};

class NNAlgorithm : public TSPAlgorithm {
    Tour giveBestTour(const Tour& initTour) {
        auto bestTour = initTour;
        for (size_t i = 0; i < bestTour.count() - 1; ++i) {
            double minRouteCost = DBL_MAX , nextRouteIndex = 0;
            for (size_t j = i + 1; j < bestTour.count(); ++j) {
                if (bestTour.distance(i,j) < minRouteCost )
                {
                    minRouteCost = bestTour.distance(i, j);
                    nextRouteIndex = j;
                }
            }
            bestTour.swap(i + 1, nextRouteIndex);
        }
        return bestTour;
    }
};

class TwoOptAlgorithm: public TSPAlgorithm {
    Tour giveBestTour(const Tour& initTour) override {
        bool isImproved;
        auto bestTour = initTour;
        do {
            isImproved = false;
            for (size_t i = 1; i < bestTour.count() - 1; ++i) {
                for (size_t j = i + 1; j < bestTour.count(); ++j) {
                    const auto ijChange = bestTour.reversed(i, j);
                    if (ijChange.tripLength() < bestTour.tripLength()) {
                        isImproved = true;
                        bestTour.reverse(i, j);
                        break;
                    }
                }
            }
        } while (isImproved);
        return bestTour;
    }
};
class ThreeOptAlgorithm: public TSPAlgorithm {
    Tour giveBestTour(const Tour& initTour) override {
        auto bestTour = initTour;
        bool isImproved = true;
        do {
            isImproved = false;
            for (size_t i = 1; i < bestTour.count() - 2; ++i) {
                for (size_t j = i + 1; j < bestTour.count() - 1 ; ++j) {
                    for (size_t k = j + 1; k < bestTour.count(); ++k) {
                        if (isBetterWhenReversed(bestTour, i, j, k)) {
                            isImproved = true;
                            break;
                        }
                    }
                }
            }
        } while (isImproved);
        return bestTour;
    }
    
    bool isBetterWhenReversed(Tour& tour, size_t i, size_t j, size_t k) {
        const Tour __ij__Reverse = tour.reversed(i, j);
        const Tour __ik__Reverse = tour.reversed(i, k);
        const Tour __jk__Reverse = tour.reversed(j, k);
        const Tour __ij_jk__Reverse = __ij__Reverse.reversed(j, k);
        const Tour __jk_ij__Reverse = __jk__Reverse.reversed(i, j);
        const Tour __ik_ij__Reverse = __ik__Reverse.reversed(i, j);
        const Tour __ik_jk__Reverse = __ik__Reverse.reversed(j, k);
        
        const std::initializer_list<Tour> compList = {tour, __ij__Reverse, __ik__Reverse, __jk__Reverse, __ij_jk__Reverse, __jk_ij__Reverse, __ik_ij__Reverse, __ik_jk__Reverse};
        
        auto lambda = [](const Tour& one, const Tour& two){
            return one.tripLength() < two.tripLength();
        };
        
        const auto locallyShortestTour = std::min(compList, lambda);
        if (locallyShortestTour.tripLength() < tour.tripLength()) {
            tour = locallyShortestTour;
            return true;
        }
        
        return  false;
    }
};



/* For future improvements */

//                    const double oldDelta = bestTour.cost(i - 1, j + 1);
//                    const double newDelta = bestTour.distance(i - 1, j) + bestTour.cost(j, i) + bestTour.distance(i, j + 1);
//                    if (newDelta < oldDelta) {
//                        isImproved = true;
//                        bestTour.reverse(i, j);
//                        break;
//                    }
