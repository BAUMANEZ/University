//
//  test_helmholtz_mpi.cpp
//  Parallels
//
//  Created by Арсений Токарев on 25.12.2022.
//

#include "test_helmholtz_mpi.hpp"

#define eps 1e-7

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
                
            case test_helmholtz_mpi::jacobi_sendrecv:
                if (rank == 0)
                    std::cout << "***** RUNNING MPI JACOBI SENDRECV ALGORITHM *****\n";
                
                configuration.run_jacobi_sendrecv();
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "@@@@@ FINISHED MPI JACOBI SENDRECV ALGORITHM @@@@@\n\n";
                
                break;
                
            case test_helmholtz_mpi::jacobi_isend_irecv:
                if (rank == 0)
                    std::cout << "***** RUNNING MPI JACOBI ISEND_IRECV ALGORITHM *****\n";
                
                configuration.run_jacobi_isend_irecv();
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "@@@@@ FINISHED MPI JACOBI ISEND_IRECV ALGORITHM @@@@@\n\n";
                
                break;
                
            case test_helmholtz_mpi::red_black_send_recv:
                if (rank == 0)
                    std::cout << "***** RUNNING MPI RED-BLACK SEND_RECV ALGORITHM *****\n";
                
                configuration.run_red_black_send_recv();
                
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "@@@@@ FINISHED MPI RED-BLACK SEND_RECV ALGORITHM @@@@@\n\n";
                
                break;
                
            case test_helmholtz_mpi::red_black_sendrecv:
                if (rank == 0)
                    std::cout << "***** RUNNING MPI RED-BLACK SENDRECV ALGORITHM *****\n";

                configuration.run_red_black_sendrecv();

                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "@@@@@ FINISHED MPI RED-BLACK SENDRECV ALGORITHM @@@@@\n\n";

                break;
                
            case test_helmholtz_mpi::red_black_isend_irecv:
                if (rank == 0)
                    std::cout << "***** RUNNING MPI RED-BLACK ISEND_IRECV ALGORITHM *****\n";

                configuration.run_red_black_isend_irecv();

                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "@@@@@ FINISHED MPI RED-BLACK ISEND_IRECV ALGORITHM @@@@@\n\n";

                break;
                
            default:
                break;
        }
    }
    
    MPI_Finalize();
}

void test_helmholtz_mpi::prepare(
    std::vector<double>& analytical,
    std::vector<int>& rows_per_process,
    std::vector<int>& shifts,
    int& recv_elements,
    std::vector<int>& elements_per_process,
    int& rows,
    int& rows_upto_rank
) {
    int size, rank;
    
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        make_analytical(analytical);
        
        for (int i = 0; i < size; ++i)
            rows_per_process[i] = (n / size);
        
        // MARK: неучтенные процессы из-за деления int
        for (int i = 0; i < n % size; ++i)
            ++rows_per_process[i];
        
        shifts[0] = 0;
        for (int i = 1; i < size; ++i)
            shifts[i] = shifts[i - 1] + rows_per_process[i - 1] * n;
        for (int i = 0; i < size; ++i)
            rows_per_process[i] += 2;
        
        rows_per_process[0] -= 1;
        rows_per_process[size - 1] -= 1;
    }
    
    MPI_Bcast(rows_per_process.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(shifts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    recv_elements = (rows_per_process[rank] - (is_first_or_last(rank, size) ? 1 : 2)) * n;
    elements_per_process = rows_per_process;
    for (int& element: elements_per_process)
        element *= n;
    
    rows = rows_per_process[rank];
    rows_upto_rank = -((rank == 0) ? 0 : rank * 2);
    for (int i = 0; i < rank; ++i)
        rows_upto_rank += rows_per_process[i];
}

void test_helmholtz_mpi::make_analytical(std::vector<double>& analytical) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            analytical[i * n + j] = this->analytical(j * h, i * h);
}

// MARK: - Jacobi Send then Receive

void test_helmholtz_mpi::run_jacobi_send_recv() {
    int size, rank;
    double start, finish;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows;
    int rows_upto_rank;
    int recv_elements;
    std::vector<int> elements_per_process;
    std::vector<double> y_per_processs;
    std::vector<double> y(n * n, 0.0);
    std::vector<double> analytical(n * n, 0.0);
    std::vector<int> rows_per_process(size, 0);
    std::vector<int> shifts(size, 0);
    
    prepare(analytical, rows_per_process, shifts, recv_elements, elements_per_process, rows, rows_upto_rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    y_per_processs.resize(rows_per_process[rank] * n, 0.0);
        
    if (rank == 0)
        start = MPI_Wtime();
        
    // главная часть
    int upper_elements = rank == (size - 1) ? 0 : n;
    int lower_elements = rank == 0 ? 0 : n;

    int upper_process = rank == (size - 1) ? 0 : rank + 1;
    int lower_process = rank == 0 ? size - 1 : rank - 1;
    
    double error;
    double error_in_process;
    int iterations = 0;
    std::vector<double> y_per_processs_prev(y_per_processs);
    
    double a = 1 / (4.0 + h * h * k * k);
    do {
        y_per_processs_prev.swap(y_per_processs);

        double* upper_row = y_per_processs_prev.data() + (rows - 1) * n;
        double* lower_row = y_per_processs_prev.data();

        MPI_Send(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 0, MPI_COMM_WORLD);
        MPI_Recv(lower_row, lower_elements, MPI_DOUBLE, lower_process, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        MPI_Send(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 0, MPI_COMM_WORLD);
        MPI_Recv(upper_row, upper_elements, MPI_DOUBLE, upper_process, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        for (int i = 1; i < rows - 1; ++i)
            for (int j = 1; j < n - 1; ++j)
                y_per_processs[i * n + j] = a * (
                    y_per_processs_prev[i * n + (j - 1)]
                    + y_per_processs_prev[(i - 1) * n + j]
                    + y_per_processs_prev[(i + 1) * n + j]
                    + y_per_processs_prev[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );

        error_in_process = 0.0;

        int i_begin = rank == 0 ? 0 : n;
        int i_end = rank == (size - 1) ? rows * n : (rows - 1) * n;

        for (int i = i_begin; i < i_end; ++i) {
            double local_error = y_per_processs[i] - y_per_processs_prev[i];
            error_in_process += local_error * local_error;
        }

        MPI_Allreduce(&error_in_process, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        error = sqrt(error);
        ++iterations;
    } while (error > eps);
     
    if (rank == 0)
        finish = MPI_Wtime();
        
    MPI_Gatherv((rank == 0) ? y_per_processs.data() : y_per_processs.data() + n, recv_elements, MPI_DOUBLE, y.data(), elements_per_process.data(), shifts.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
        std::cout << "Time=" << finish - start << ";  Iterations=" << iterations << ";  ||frobenius||=" << norm(analytical, y) << "\n";
}

// MARK: - Jacobi SendReceive

void test_helmholtz_mpi::run_jacobi_sendrecv() {
    int size, rank;
    double start, finish;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows;
    int rows_upto_rank;
    int recv_elements;
    std::vector<int> elements_per_process;
    std::vector<double> y_per_processs;
    std::vector<double> y(n * n, 0.0);
    std::vector<double> analytical(n * n, 0.0);
    std::vector<int> rows_per_process(size, 0);
    std::vector<int> shifts(size, 0);
    
    prepare(analytical, rows_per_process, shifts, recv_elements, elements_per_process, rows, rows_upto_rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    y_per_processs.resize(rows_per_process[rank] * n, 0.0);
        
    if (rank == 0)
        start = MPI_Wtime();
        
    // главная часть
    int upper_elements = rank == (size - 1) ? 0 : n;
    int lower_elements = rank == 0 ? 0 : n;

    int upper_process = rank == (size - 1) ? 0 : rank + 1;
    int lower_process = rank == 0 ? size - 1 : rank - 1;
    
    double error;
    double error_in_process;
    int iterations = 0;
    std::vector<double> y_per_processs_prev(y_per_processs);
    
    double a = 1 / (4.0 + h * h * k * k);
    do {
        y_per_processs_prev.swap(y_per_processs);

        double* upper_row = y_per_processs_prev.data() + (rows - 1) * n;
        double* lower_row = y_per_processs_prev.data();

        MPI_Sendrecv(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 1, y_per_processs_prev.data(), lower_elements, MPI_DOUBLE, lower_process, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Sendrecv(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 2, upper_row, upper_elements, MPI_DOUBLE, upper_process, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        for (int i = 1; i < rows - 1; ++i)
            for (int j = 1; j < n - 1; ++j)
                y_per_processs[i * n + j] = a * (
                    y_per_processs_prev[i * n + (j - 1)]
                    + y_per_processs_prev[(i - 1) * n + j]
                    + y_per_processs_prev[(i + 1) * n + j]
                    + y_per_processs_prev[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );

        error_in_process = 0.0;

        int i_begin = rank == 0 ? 0 : n;
        int i_end = rank == (size - 1) ? rows * n : (rows - 1) * n;

        for (int i = i_begin; i < i_end; ++i) {
            double local_error = y_per_processs[i] - y_per_processs_prev[i];
            error_in_process += local_error * local_error;
        }

        MPI_Allreduce(&error_in_process, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        error = sqrt(error);
        ++iterations;
    } while (error > eps);
     
    if (rank == 0)
        finish = MPI_Wtime();
        
    MPI_Gatherv((rank == 0) ? y_per_processs.data() : y_per_processs.data() + n, recv_elements, MPI_DOUBLE, y.data(), elements_per_process.data(), shifts.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
        std::cout << "Time=" << finish - start << ";  Iterations=" << iterations << ";  ||frobenius||=" << norm(analytical, y) << "\n";
}

// MARK: - Jacobi Isend Then IRecieve

void test_helmholtz_mpi::run_jacobi_isend_irecv() {
    int size, rank;
    double start, finish;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows;
    int rows_upto_rank;
    int recv_elements;
    std::vector<int> elements_per_process;
    std::vector<double> y_per_processs;
    std::vector<double> y(n * n, 0.0);
    std::vector<double> analytical(n * n, 0.0);
    std::vector<int> rows_per_process(size, 0);
    std::vector<int> shifts(size, 0);
    
    prepare(analytical, rows_per_process, shifts, recv_elements, elements_per_process, rows, rows_upto_rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    y_per_processs.resize(rows_per_process[rank] * n, 0.0);
        
    if (rank == 0)
        start = MPI_Wtime();
        
    // главная часть
    int upper_elements = rank == (size - 1) ? 0 : n;
    int lower_elements = rank == 0 ? 0 : n;

    int upper_process = rank == (size - 1) ? 0 : rank + 1;
    int lower_process = rank == 0 ? size - 1 : rank - 1;
    
    double error;
    double error_in_process;
    int iterations = 0;
    std::vector<double> y_per_processs_prev(y_per_processs);
    
    MPI_Request upper_send, upper_receive;
    MPI_Request lower_send, lower_receive;
    
    double* upper_row = y_per_processs.data() + (rows - 1) * n;
    double* lower_row = y_per_processs.data();
    
    MPI_Send_init(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 1, MPI_COMM_WORLD, &upper_send);
    MPI_Recv_init(lower_row, lower_elements, MPI_DOUBLE, lower_process, 1, MPI_COMM_WORLD, &lower_receive);

    MPI_Send_init(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 2, MPI_COMM_WORLD, &lower_send);
    MPI_Recv_init(upper_row, upper_elements, MPI_DOUBLE, upper_process, 2, MPI_COMM_WORLD, &upper_receive);
    
    MPI_Request upper_send_prev, upper_receive_prev;
    MPI_Request lower_send_prev, lower_receive_prev;
    
    upper_row = y_per_processs_prev.data() + (rows - 1) * n;
    lower_row = y_per_processs_prev.data();
    
    MPI_Send_init(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 3, MPI_COMM_WORLD, &upper_send_prev);
    MPI_Recv_init(lower_row, lower_elements, MPI_DOUBLE, lower_process, 3, MPI_COMM_WORLD, &lower_receive_prev);

    MPI_Send_init(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 4, MPI_COMM_WORLD, &lower_send_prev);
    MPI_Recv_init(upper_row, upper_elements, MPI_DOUBLE, upper_process, 4, MPI_COMM_WORLD, &upper_receive_prev);
    
    double a = 1 / (4.0 + h * h * k * k);
    do {
        y_per_processs_prev.swap(y_per_processs);
        
        if (!is_odd(iterations)) {
            MPI_Start(&upper_receive);
            MPI_Start(&lower_receive);
            MPI_Start(&upper_send);
            MPI_Start(&lower_send);
        } else {
            MPI_Start(&upper_receive_prev);
            MPI_Start(&lower_receive_prev);
            MPI_Start(&upper_send_prev);
            MPI_Start(&lower_send_prev);
        }
        
        int i = 2;
        for (; i < rows - 2; ++i)
            for (int j = 1; j < n - 1; ++j)
                y_per_processs[i * n + j] = a * (
                    y_per_processs_prev[i * n + (j - 1)]
                    + y_per_processs_prev[(i - 1) * n + j]
                    + y_per_processs_prev[(i + 1) * n + j]
                    + y_per_processs_prev[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );

        if (!is_odd(iterations)) {
            MPI_Wait(&upper_receive, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_receive, MPI_STATUSES_IGNORE);
            MPI_Wait(&upper_send, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_send, MPI_STATUSES_IGNORE);
        } else {
            MPI_Wait(&upper_receive_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_receive_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&upper_send_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_send_prev, MPI_STATUSES_IGNORE);
        }
        
        i = 1;
        for (int j = 1; j < n - 1; ++j)
            y_per_processs[i * n + j] = a * (
                y_per_processs_prev[i * n + (j - 1)]
                + y_per_processs_prev[(i - 1) * n + j]
                + y_per_processs_prev[(i + 1) * n + j]
                + y_per_processs_prev[i * n + (j + 1)]
                + h * h * f(j * h, (i + rows_upto_rank) * h)
            );
        
        i = rows - 2;
        for (int j = 1; j < n - 1; ++j)
            y_per_processs[i * n + j] = a * (
                y_per_processs_prev[i * n + (j - 1)]
                + y_per_processs_prev[(i - 1) * n + j]
                + y_per_processs_prev[(i + 1) * n + j]
                + y_per_processs_prev[i * n + (j + 1)]
                + h * h * f(j * h, (i + rows_upto_rank) * h)
            );
        
        error_in_process = 0.0;

        int i_begin = rank == 0 ? 0 : n;
        int i_end = rank == (size - 1) ? rows * n : (rows - 1) * n;

        for (int i = i_begin; i < i_end; ++i) {
            double local_error = y_per_processs[i] - y_per_processs_prev[i];
            error_in_process += local_error * local_error;
        }

        MPI_Allreduce(&error_in_process, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        error = sqrt(error);
        ++iterations;
    } while (error > eps);
    
    if (rank == 0)
        finish = MPI_Wtime();
        
    MPI_Gatherv((rank == 0) ? y_per_processs.data() : y_per_processs.data() + n, recv_elements, MPI_DOUBLE, y.data(), elements_per_process.data(), shifts.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
        std::cout << "Time=" << finish - start << ";  Iterations=" << iterations << ";  ||frobenius||=" << norm(analytical, y) << "\n";
}

// MARK: - Red-Black Send Then Receive
void test_helmholtz_mpi::run_red_black_send_recv() {
    int size, rank;
    double start, finish;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows;
    int rows_upto_rank;
    int recv_elements;
    std::vector<int> elements_per_process;
    std::vector<double> y_per_processs;
    std::vector<double> y(n * n, 0.0);
    std::vector<double> analytical(n * n, 0.0);
    std::vector<int> rows_per_process(size, 0);
    std::vector<int> shifts(size, 0);
    
    prepare(analytical, rows_per_process, shifts, recv_elements, elements_per_process, rows, rows_upto_rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    y_per_processs.resize(rows_per_process[rank] * n, 0.0);
        
    if (rank == 0)
        start = MPI_Wtime();
        
    // главная часть
    int upper_elements = rank == (size - 1) ? 0 : n;
    int lower_elements = rank == 0 ? 0 : n;

    int upper_process = rank == (size - 1) ? 0 : rank + 1;
    int lower_process = rank == 0 ? size - 1 : rank - 1;
    
    double error;
    double error_in_process;
    int iterations = 0;
    std::vector<double> y_per_processs_prev(y_per_processs);
    
    double a = 1 / (4.0 + h * h * k * k);
    do {
        y_per_processs_prev.swap(y_per_processs);

        double* upper_row = y_per_processs_prev.data() + (rows - 1) * n;
        double* lower_row = y_per_processs_prev.data();

        MPI_Send(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 0, MPI_COMM_WORLD);
        MPI_Recv(lower_row, lower_elements, MPI_DOUBLE, lower_process, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        MPI_Send(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 0, MPI_COMM_WORLD);
        MPI_Recv(upper_row, upper_elements, MPI_DOUBLE, upper_process, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        
        for (int i = 1; i < rows - 1; ++i)
            for (int j = (i + rows) % 2 + 1; j < n - 1; j += 2)
                y_per_processs[i * n + j] = a * (
                    y_per_processs_prev[i * n + (j - 1)]
                    + y_per_processs_prev[(i - 1) * n + j]
                    + y_per_processs_prev[(i + 1) * n + j]
                    + y_per_processs_prev[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );
        
        upper_row = y_per_processs.data() + (rows - 1) * n;
        lower_row = y_per_processs.data();
        
        MPI_Send(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 0, MPI_COMM_WORLD);
        MPI_Recv(lower_row, lower_elements, MPI_DOUBLE, lower_process, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

        MPI_Send(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 0, MPI_COMM_WORLD);
        MPI_Recv(upper_row, upper_elements, MPI_DOUBLE, upper_process, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        
        for (int i = 1; i < rows - 1; ++i)
            for (int j = (i + rows + 1) % 2 + 1; j < n - 1; j += 2)
                y_per_processs[i * n + j] = a * (
                    y_per_processs[i * n + (j - 1)]
                    + y_per_processs[(i - 1) * n + j]
                    + y_per_processs[(i + 1) * n + j]
                    + y_per_processs[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );
        
        error_in_process = 0.0;

        int i_begin = rank == 0 ? 0 : n;
        int i_end = rank == (size - 1) ? rows * n : (rows - 1) * n;

        for (int i = i_begin; i < i_end; ++i) {
            double local_error = y_per_processs[i] - y_per_processs_prev[i];
            error_in_process += local_error * local_error;
        }

        MPI_Allreduce(&error_in_process, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        error = sqrt(error);
        ++iterations;
    } while (error > eps);
     
    if (rank == 0)
        finish = MPI_Wtime();
        
    MPI_Gatherv((rank == 0) ? y_per_processs.data() : y_per_processs.data() + n, recv_elements, MPI_DOUBLE, y.data(), elements_per_process.data(), shifts.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
        std::cout << "Time=" << finish - start << ";  Iterations=" << iterations << ";  ||frobenius||=" << norm(analytical, y) << "\n";
}

// MARK: - Red-Black SendReceive
void test_helmholtz_mpi::run_red_black_sendrecv() {
    int size, rank;
    double start, finish;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows;
    int rows_upto_rank;
    int recv_elements;
    std::vector<int> elements_per_process;
    std::vector<double> y_per_processs;
    std::vector<double> y(n * n, 0.0);
    std::vector<double> analytical(n * n, 0.0);
    std::vector<int> rows_per_process(size, 0);
    std::vector<int> shifts(size, 0);
    
    prepare(analytical, rows_per_process, shifts, recv_elements, elements_per_process, rows, rows_upto_rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    y_per_processs.resize(rows_per_process[rank] * n, 0.0);
        
    if (rank == 0)
        start = MPI_Wtime();
        
    // главная часть
    int upper_elements = rank == (size - 1) ? 0 : n;
    int lower_elements = rank == 0 ? 0 : n;

    int upper_process = rank == (size - 1) ? 0 : rank + 1;
    int lower_process = rank == 0 ? size - 1 : rank - 1;
    
    double error;
    double error_in_process;
    int iterations = 0;
    std::vector<double> y_per_processs_prev(y_per_processs);
    
    double a = 1 / (4.0 + h * h * k * k);
    do {
        y_per_processs_prev.swap(y_per_processs);

        double* upper_row = y_per_processs_prev.data() + (rows - 1) * n;
        double* lower_row = y_per_processs_prev.data();

        MPI_Sendrecv(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 1, lower_row, lower_elements, MPI_DOUBLE, lower_process, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Sendrecv(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 2, upper_row, upper_elements, MPI_DOUBLE, upper_process, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        
        for (int i = 1; i < rows - 1; ++i)
            for (int j = (i + rows) % 2 + 1; j < n - 1; j += 2)
                y_per_processs[i * n + j] = a * (
                    y_per_processs_prev[i * n + (j - 1)]
                    + y_per_processs_prev[(i - 1) * n + j]
                    + y_per_processs_prev[(i + 1) * n + j]
                    + y_per_processs_prev[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );
        
        upper_row = y_per_processs.data() + (rows - 1) * n;
        lower_row = y_per_processs.data();
        
        MPI_Sendrecv(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 1, lower_row, lower_elements, MPI_DOUBLE, lower_process, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        MPI_Sendrecv(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 2, upper_row, upper_elements, MPI_DOUBLE, upper_process, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        
        for (int i = 1; i < rows - 1; ++i)
            for (int j = (i + rows + 1) % 2 + 1; j < n - 1; j += 2)
                y_per_processs[i * n + j] = a * (
                    y_per_processs[i * n + (j - 1)]
                    + y_per_processs[(i - 1) * n + j]
                    + y_per_processs[(i + 1) * n + j]
                    + y_per_processs[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );
        
        error_in_process = 0.0;

        int i_begin = rank == 0 ? 0 : n;
        int i_end = rank == (size - 1) ? rows * n : (rows - 1) * n;

        for (int i = i_begin; i < i_end; ++i) {
            double local_error = y_per_processs[i] - y_per_processs_prev[i];
            error_in_process += local_error * local_error;
        }

        MPI_Allreduce(&error_in_process, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        error = sqrt(error);
        ++iterations;
    } while (error > eps);
     
    if (rank == 0)
        finish = MPI_Wtime();
        
    MPI_Gatherv((rank == 0) ? y_per_processs.data() : y_per_processs.data() + n, recv_elements, MPI_DOUBLE, y.data(), elements_per_process.data(), shifts.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
        std::cout << "Time=" << finish - start << ";  Iterations=" << iterations << ";  ||frobenius||=" << norm(analytical, y) << "\n";
}

// MARK: - Red-Black Isend Then IRecieve

void test_helmholtz_mpi::run_red_black_isend_irecv() {
    int size, rank;
    double start, finish;
        
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows;
    int rows_upto_rank;
    int recv_elements;
    std::vector<int> elements_per_process;
    std::vector<double> y_per_processs;
    std::vector<double> y(n * n, 0.0);
    std::vector<double> analytical(n * n, 0.0);
    std::vector<int> rows_per_process(size, 0);
    std::vector<int> shifts(size, 0);
    
    prepare(analytical, rows_per_process, shifts, recv_elements, elements_per_process, rows, rows_upto_rank);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    y_per_processs.resize(rows_per_process[rank] * n, 0.0);
        
    if (rank == 0)
        start = MPI_Wtime();
        
    // главная часть
    int upper_elements = rank == (size - 1) ? 0 : n;
    int lower_elements = rank == 0 ? 0 : n;

    int upper_process = rank == (size - 1) ? 0 : rank + 1;
    int lower_process = rank == 0 ? size - 1 : rank - 1;
    
    double error;
    double error_in_process;
    int iterations = 0;
    std::vector<double> y_per_processs_prev(y_per_processs);
    
    MPI_Request upper_send, upper_receive;
    MPI_Request lower_send, lower_receive;
    
    double* upper_row = y_per_processs.data() + (rows - 1) * n;
    double* lower_row = y_per_processs.data();
    
    MPI_Send_init(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 1, MPI_COMM_WORLD, &upper_send);
    MPI_Recv_init(lower_row, lower_elements, MPI_DOUBLE, lower_process, 1, MPI_COMM_WORLD, &lower_receive);

    MPI_Send_init(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 2, MPI_COMM_WORLD, &lower_send);
    MPI_Recv_init(upper_row, upper_elements, MPI_DOUBLE, upper_process, 2, MPI_COMM_WORLD, &upper_receive);
    
    MPI_Request upper_send_prev, upper_receive_prev;
    MPI_Request lower_send_prev, lower_receive_prev;
    
    upper_row = y_per_processs_prev.data() + (rows - 1) * n;
    lower_row = y_per_processs_prev.data();
    
    MPI_Send_init(upper_row - n, upper_elements, MPI_DOUBLE, upper_process, 3, MPI_COMM_WORLD, &upper_send_prev);
    MPI_Recv_init(lower_row, lower_elements, MPI_DOUBLE, lower_process, 3, MPI_COMM_WORLD, &lower_receive_prev);

    MPI_Send_init(lower_row + n, lower_elements, MPI_DOUBLE, lower_process, 4, MPI_COMM_WORLD, &lower_send_prev);
    MPI_Recv_init(upper_row, upper_elements, MPI_DOUBLE, upper_process, 4, MPI_COMM_WORLD, &upper_receive_prev);
    
    double a = 1 / (4.0 + h * h * k * k);
    do {
        y_per_processs_prev.swap(y_per_processs);
        
        int i = 1;
        for (int j = (i + rows) % 2 + 1; j < n - 1; j += 2)
            y_per_processs[i * n + j] = a * (
                y_per_processs_prev[i * n + (j - 1)]
                + y_per_processs_prev[(i - 1) * n + j]
                + y_per_processs_prev[(i + 1) * n + j]
                + y_per_processs_prev[i * n + (j + 1)]
                + h * h * f(j * h, (i + rows_upto_rank) * h)
            );
        
        i = rows - 2;
        for (int j = (i + rows) % 2 + 1; j < n - 1; j += 2)
            y_per_processs[i * n + j] = a * (
                y_per_processs_prev[i * n + (j - 1)]
                + y_per_processs_prev[(i - 1) * n + j]
                + y_per_processs_prev[(i + 1) * n + j]
                + y_per_processs_prev[i * n + (j + 1)]
                + h * h * f(j * h, (i + rows_upto_rank) * h)
            );
        
        if (is_odd(iterations)) {
            MPI_Start(&upper_receive);
            MPI_Start(&lower_receive);
            MPI_Start(&upper_send);
            MPI_Start(&lower_send);
        } else {
            MPI_Start(&upper_receive_prev);
            MPI_Start(&lower_receive_prev);
            MPI_Start(&upper_send_prev);
            MPI_Start(&lower_send_prev);
        }
        
        i = 2;
        for (; i < rows - 2; ++i)
            for (int j = (i + rows) % 2 + 1; j < n - 1; j += 2)
                y_per_processs[i * n + j] = a * (
                    y_per_processs_prev[i * n + (j - 1)]
                    + y_per_processs_prev[(i - 1) * n + j]
                    + y_per_processs_prev[(i + 1) * n + j]
                    + y_per_processs_prev[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );

        if (is_odd(iterations)) {
            MPI_Wait(&upper_receive, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_receive, MPI_STATUSES_IGNORE);
            MPI_Wait(&upper_send, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_send, MPI_STATUSES_IGNORE);
        } else {
            MPI_Wait(&upper_receive_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_receive_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&upper_send_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_send_prev, MPI_STATUSES_IGNORE);
        }
        
        i = 1;
        for (int j = (i + rows + 1) % 2 + 1; j < n - 1; j += 2)
            y_per_processs[i * n + j] = a * (
                y_per_processs[i * n + (j - 1)]
                + y_per_processs[(i - 1) * n + j]
                + y_per_processs[(i + 1) * n + j]
                + y_per_processs[i * n + (j + 1)]
                + h * h * f(j * h, (i + rows_upto_rank) * h)
            );
        
        i = rows - 2;
        for (int j = (i + rows + 1) % 2 + 1; j < n - 1; j += 2)
            y_per_processs[i * n + j] = a * (
                y_per_processs[i * n + (j - 1)]
                + y_per_processs[(i - 1) * n + j]
                + y_per_processs[(i + 1) * n + j]
                + y_per_processs[i * n + (j + 1)]
                + h * h * f(j * h, (i + rows_upto_rank) * h)
            );
        
        if (is_odd(iterations)) {
            MPI_Start(&upper_receive);
            MPI_Start(&lower_receive);
            MPI_Start(&upper_send);
            MPI_Start(&lower_send);
        } else {
            MPI_Start(&upper_receive_prev);
            MPI_Start(&lower_receive_prev);
            MPI_Start(&upper_send_prev);
            MPI_Start(&lower_send_prev);
        }
        
        i = 2;
        for (; i < rows - 2; ++i)
            for (int j = (i + rows + 1) % 2 + 1; j < n - 1; j += 2)
                y_per_processs[i * n + j] = a * (
                    y_per_processs[i * n + (j - 1)]
                    + y_per_processs[(i - 1) * n + j]
                    + y_per_processs[(i + 1) * n + j]
                    + y_per_processs[i * n + (j + 1)]
                    + h * h * f(j * h, (i + rows_upto_rank) * h)
                );

        if (is_odd(iterations)) {
            MPI_Wait(&upper_receive, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_receive, MPI_STATUSES_IGNORE);
            MPI_Wait(&upper_send, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_send, MPI_STATUSES_IGNORE);
        } else {
            MPI_Wait(&upper_receive_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_receive_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&upper_send_prev, MPI_STATUSES_IGNORE);
            MPI_Wait(&lower_send_prev, MPI_STATUSES_IGNORE);
        }
        
        error_in_process = 0.0;

        int i_begin = rank == 0 ? 0 : n;
        int i_end = rank == (size - 1) ? rows * n : (rows - 1) * n;

        for (int i = i_begin; i < i_end; ++i) {
            double local_error = y_per_processs[i] - y_per_processs_prev[i];
            error_in_process += local_error * local_error;
        }

        MPI_Allreduce(&error_in_process, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        error = sqrt(error);
        ++iterations;
    } while (error > eps);
    
    if (rank == 0)
        finish = MPI_Wtime();
        
    MPI_Gatherv((rank == 0) ? y_per_processs.data() : y_per_processs.data() + n, recv_elements, MPI_DOUBLE, y.data(), elements_per_process.data(), shifts.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0)
        std::cout << "Time=" << finish - start << ";  Iterations=" << iterations << ";  ||frobenius||=" << norm(analytical, y) << "\n";
}

test_helmholtz_mpi::test_helmholtz_mpi(int i) {
    const int t = 320 * i;
    n =  t * 10;
    h = 1.0 / double(n - 1);
    k = double(t) * 20.0;
}
