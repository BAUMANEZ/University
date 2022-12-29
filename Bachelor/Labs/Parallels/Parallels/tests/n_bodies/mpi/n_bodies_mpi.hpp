//
//  n_bodies_mpi.hpp
//  Parallels
//
//  Created by Арсений Токарев on 27.12.2022.
//

#ifndef n_bodies_mpi_hpp
#define n_bodies_mpi_hpp

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include "mpi.h"
#include "../../../helpers/helpers.hpp"

struct n_bodies_mpi {
    
public:
    struct Coordinate {
        double x1 = 0;
        double x2 = 0;
        double x3 = 0;
    };
    
    struct Body {
        Coordinate r;
        Coordinate v;
        double m = 0;
    };
    
    static void run();
    
private:
    const int n = 4;
    const string input_file_name = "4body.txt";
    const string generated_data_file_name = "100body.txt";
    const string output_file_name = "result.txt";
    
    const double G = 6.67 * 1e-11;
    const double T = 2.0; // конечное время
    
    void load_file();    
    double norm(const Coordinate& lhs, const Coordinate& rhs);
};

#endif /* n_bodies_mpi_hpp */
