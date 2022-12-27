//
//  main.cpp
//  Parallels
//
//  Created by Арсений Токарев on 07.10.2022.
//

#include "../tests/helmholtz/mpi/test_helmholtz_mpi.hpp"

int main(int argc, char** argv) {
    test_helmholtz_mpi::run(&argc, &argv, 1, {
        test_helmholtz_mpi::jacobi_send_recv,
        test_helmholtz_mpi::jacobi_sendrecv,
        test_helmholtz_mpi::jacobi_isend_irecv,
        test_helmholtz_mpi::red_black_send_recv,
        test_helmholtz_mpi::red_black_sendrecv,
        test_helmholtz_mpi::red_black_isend_irecv
    });

    return 0;
}
