
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include<string>

typedef std::vector<double> NumericVector;
typedef std::vector<NumericVector> NumericMatrix;

struct Tour {
public:
    double operator[](size_t i) const {
        return nodes[i];
    }
    
    void operator=(Tour& other) {
        this->nodes = other.nodes;
    }
    
    void operator=(const Tour& other) {
        this->nodes = other.nodes;
    }
    
    void add(double route) {
        this->nodes.push_back(route);
    }
    
    void replace(const size_t index, const double route) {
        this->nodes[index] = route;
    }
    
    void swap(size_t i, size_t j) {
        std::swap(this->nodes[i], this->nodes[j]);
    }
    
    double distance(const size_t from, const size_t to) const {
        return travelCostMatrix[nodes[from]][nodes[to]];
    }
    
    size_t count() const {
        return nodes.size();
    }
    
    double cost(const size_t start, const size_t end) const {
        double totalCost = 0;
        if (start <= end) {
            for (auto i = start; i < end; ++i) {
                const size_t from = (*this)[i];
                const size_t to = (*this)[ (i + 1) % this->count() ];
                totalCost += travelCostMatrix[from][to];
            }
        } else {
            for (auto i = end; i > start; --i) {
                const size_t from = (*this)[i];
                const size_t to = (*this)[ (i - 1) % this->count() ];
                totalCost += travelCostMatrix[from][to];
            }
        }
        return totalCost;
    }
    
    Tour reversed(size_t start, size_t end) const {
        Tour resultTour {*this};
        while (start < end) {
            std::swap(resultTour.nodes[start], resultTour.nodes[end]);
            ++start; --end;
        }
        return resultTour;
    }
    
    void reverse(size_t start, size_t end) {
        while (start < end) {
            std::swap(this->nodes[start], this->nodes[end]);
            ++start; --end;
        }
    }
    
    double tripLength() const {
        return cost(0, nodes.size());
    }
    

    Tour(const std::string& filePath = "") {
        /*
         
         Initializing a file to assign routes and the number of nodes
         
         Note: originally a 'filePath' variable is string typed and it's default value can be edited in main.cpp file
         
        */
        
        size_t numberOfNodes = 0;
        std::ifstream input(filePath);
        if (input.is_open()) {
            input >> numberOfNodes;
            travelCostMatrix.reserve(numberOfNodes);
            for (size_t i = 0; i < numberOfNodes; ++i) {
                NumericVector temp;
                temp.reserve(numberOfNodes);
                for (size_t j = 0; j < numberOfNodes; ++j) {
                    double value;
                    input >> value;
                    temp.push_back(value);
                }
                travelCostMatrix.push_back(temp);
            }
        }
        input.close();
        
        /*
            An initial route is usually represented in an arithmetic order (from the first node to the last one)
         
        */
        
        nodes.reserve(numberOfNodes);
        for (size_t node = 0; node < numberOfNodes; ++node) {
            nodes.push_back(node); // 0 is the first node
        }
    }
    
private:
    NumericVector nodes;
    NumericMatrix travelCostMatrix;
};
