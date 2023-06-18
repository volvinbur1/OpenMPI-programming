#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <array>

std::array<std::string, 3> get_parameters(int curr_rank, int nodes_cnt) {
    
    return std::array<std::string, 3>();
}

double calculate_integral(int curr_rank, int nodes_cnt, std::array<std::string, 3> params) {
    return 0;
}

int main(int argc, char *argv[]) {
        MPI_Init(&argc, &argv);

        int curr_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

        int nodes_cnt;
        MPI_Comm_size(MPI_COMM_WORLD, &nodes_cnt);

        std::cout << "[MPI_INIT_INFO] Current node rank: " << curr_rank << "\tTotal nodes count: " << nodes_cnt << std::endl;

        auto params = get_parameters(curr_rank, nodes_cnt);
        auto result = calculate_integral(curr_rank, nodes_cnt, params);

        MPI_Finalize();
        return 0;
}
