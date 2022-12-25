//
//  test_helmholtz_mpi.hpp
//  Parallels
//
//  Created by Арсений Токарев on 25.12.2022.
//

#ifndef test_helmholtz_mpi_hpp
#define test_helmholtz_mpi_hpp

#include <iostream>
#include <vector>
#include <math.h>
#include "mpi.h"
#include "../../../helpers/helpers.hpp"

struct test_helmholtz_mpi {
    
public:
    enum test_case {
        jacobi_send_recv,
        jacobi_sendrecv,
        jacobi_isend_irecv,
        red_black_send_recv,
        red_black_sendrecv,
        red_black_isend_irecv,
    };
    
    static void run(int* argc, char*** argv, int i, std::vector<test_case> tests);
    
private:
    static void prepare(
        std::vector<double>& matrix,
        std::vector<int>& rows_per_process
        const test_helmholtz_mpi& configuration
    );
    
    void run_jacobi_send_recv();
    void run_jacobi_sendrecv();
    void run_jacobi_isend_irecv();
    
    void run_red_black_send_recv();
    void run_red_black_sendrecv();
    void run_red_black_isend_irecv();
    
    int n;
    double h;
    double k;
    
    double f(double x, double y);
    double analytical(double x, double y);
    double norm(const std::vector<double>& lhs, const std::vector<double>& rhs);
    
    test_helmholtz_mpi(int i);
};

#endif /* test_helmholtz_mpi_hpp */
