#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <chrono>

#define PI atan(1)*4
#define ONE_BY_SQRT_2 1/sqrt(2)

int dboard(int N);

int main(int argc, char* argv[]){
    // validate input
    if (argc != 3) {
        std::cout << "input invalid!" << std::endl;
        return 0;
    }
    std::size_t pos;
    int N = std::stoi(argv[1], &pos);
    int R = std::stoi(argv[2], &pos);
    // set up MPI
    MPI_Init(&argc, &argv);
    // get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    srand(rank);

    int m, m_sum;
    double master_PI_sum = 0.0;
    auto t_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < R; i++) {
        m = dboard(N/p);
        MPI_Reduce(&m, &m_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            master_PI_sum += 2.0 * N / m_sum;
        }
    }
    if (rank == 0) {
        master_PI_sum /= (double)R;
        std::cout << "N=" << N << ", R=" << R << 
            ", P=" << p << ", PI=" << master_PI_sum << std::endl;
        auto t_end = std::chrono::high_resolution_clock::now();
        std::cout << "Time=" << 
            std::chrono::duration<double>(t_end-t_start).count() << "sec" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

int dboard(int N) {
    // Compute number of darts to be thrown locally
    //Throw darts (Generate random coords (x,y) )
    //check whether dart lands in square
    int m = 0;
    double a, theta;
    for (int i = 0; i < N; ++i) {
        a = (double)rand() / RAND_MAX;
        theta = (double)rand() / RAND_MAX * 2 * PI;
        if (abs(sqrt(a)*cos(theta)) <= ONE_BY_SQRT_2 &&
            abs(sqrt(a)*sin(theta)) <= ONE_BY_SQRT_2) {
            m++;
        }
    }
    return m;
}
