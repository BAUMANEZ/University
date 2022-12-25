//
//  test_helmholtz_mpi.cpp
//  Parallels
//
//  Created by Арсений Токарев on 25.12.2022.
//

#include "test_helmholtz_mpi.hpp"

#define eps = 1e-7

double test_helmholtz_mpi::f(double x, double y) {
    return (2 * sin(M_PI * y) + k * k * (1 - x) * x * sin(M_PI * y) + M_PI * M_PI * (1 - x) * x * sin(M_PI * y));
}

double test_helmholtz_mpi::analytical(double x, double y) {
    return (1 - x) * x * sin(M_PI * y);
}

double test_helmholtz_mpi::norm(const std::vector<double>& lhs, const std::vector<double>& rhs) {
    double sum = 0.;
    for (size_t i = 0; i < n * n; ++i) {
        const double diff = lhs[i] - rhs[i];
        sum += diff*diff;
    }
    
    return sqrt(sum);
}

void test_helmholtz_mpi::run(int* argc, char*** argv, int i, std::vector<test_case> tests) {
    int size, rank;
    
    MPI_Init(argc, argv);
    
    for (auto& test: tests) {
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        MPI_Barrier(MPI_COMM_WORLD);
        test_helmholtz_mpi configuration(i);
        
        MPI_Barrier(MPI_COMM_WORLD);
        switch (test) {
            case test_helmholtz_mpi::jacobi_send_recv:
                if (rank == 0)
                    std::cout << "***** RUNNING MPI JACOBI SEND_RECV ALGORITHM *****\n";
                
                configuration.run_jacobi_send_recv();
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "@@@@@ FINISHED MPI JACOBI SEND_RECV ALGORITHM @@@@@\n\n";
                
                break;
                
            default:
                break;
                
//            case test_helmholtz_mpi::jacobi_sendrecv:
//                configuration.run_jacobi_sendrecv();
//                break;
//
//            case test_helmholtz_mpi::jacobi_isend_irecv:
//                configuration.run_jacobi_isend_irecv();
//                break;
//
//            case test_helmholtz_mpi::red_black_send_recv:
//                configuration.run_red_black_send_recv();
//                break;
//
//            case test_helmholtz_mpi::red_black_sendrecv:
//                configuration.run_red_black_sendrecv();
//                break;
//
//            case test_helmholtz_mpi::red_black_isend_irecv:
//                configuration.run_red_black_isend_irecv();
//                break;
        }
    }
    
    MPI_Finalize();
}

void test_helmholtz_mpi::prepare(
    std::vector<double>& matrix,
    const test_helmholtz_mpi& configuration
) {
    int size, rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    
}

void test_helmholtz_mpi::run_jacobi_send_recv() {
    int size, rank;
    int start;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        start = MPI_Wtime();
    }
}

test_helmholtz_mpi::test_helmholtz_mpi(int i) {
    const int t = 320 * i;
    n =  t * 10;
    h = 1.0 / double(n - 1);
    k = double(t) * 20.0;
}
