#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <array>
#include <string>
#include <vector>
#include <chrono>

using params_array = std::array<double, 3>;

params_array read_params_from_file() {
    params_array params;
    try {
        std::ifstream file("input.txt");
        
        for (auto i = 0; i < 3; i++) {
            std::string line;
            std::getline(file, line);
            params.at(i) = std::stod(line);
        }

        file.close();
    } catch (std::exception& e) {
        std::cout << "[FATAL] " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    return params;
}

params_array get_parameters(int curr_rank) {
    params_array params;
    if (curr_rank == 0) {
        params = read_params_from_file();
        std::cout << "[HOST - 0] Input values. a: " << params.at(0) 
            << "\tb: " << params.at(1) 
            << "\taccuracy: " << params.at(2) 
            << std::endl;
    }
    
    MPI_Bcast(params.data(), 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return params;
}

bool validate_input(params_array params) {
    if (params.at(0) > params.at(1)) {
        std::cout << "[ERROR] Parameter A is greater then B. Exit..." << std::endl;
        return false;
    }

    if (params.at(2) < 0 || params.at(2) > 1) {
        std::cout << "[ERROR] Accuracy parameter in out of (0,1] range. Exit..." << std::endl;
        return false;
    }

    return true;
}

bool is_accurate_enough(double new_val, double prev_val, double accuracy) {
    return fabs(new_val - prev_val) / 3. < accuracy;
}

double integral_function(double x) {
    return log2(pow(x, 3));
}

double integrate_right_rectangle(double a_edge, double b_edge, double accuracy) {
    int itr_cnt = 1;
    double new_res = -1, prev_res = 0;

    while (!is_accurate_enough(new_res, prev_res, accuracy)) {
        itr_cnt *= 2;
        prev_res = new_res;
        new_res = 0;
        auto h = (b_edge - a_edge) / itr_cnt;
        for (auto i = 1; i < itr_cnt; i++) {
            new_res += integral_function(a_edge + i * h) * h;
        }
    }
    return new_res;
}

double calculate_integral(int curr_rank, int nodes_cnt, params_array params) {
    auto rank_step = (params.at(1) - params.at(0)) / nodes_cnt;
    auto rank_a_edge = params.at(0) + curr_rank * rank_step;
    auto rank_b_edge = params.at(0) + (curr_rank + 1) * rank_step;
    auto rank_result = integrate_right_rectangle(rank_a_edge, rank_b_edge, params.at(2));

    auto total_integrall_value = rank_result;
    std::vector<double> nodes_results(nodes_cnt - 1);
    if (curr_rank != 0) {
        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Isend(&rank_result, 1, MPI_DOUBLE, 0, curr_rank, MPI_COMM_WORLD, &request);
        MPI_Waitall(1, &request, nullptr);
    } else {
        MPI_Request requests[nodes_cnt - 1];
        for (auto i = 1; i < nodes_cnt; i++) {
            MPI_Irecv(nodes_results.data() + (i - 1), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &requests[i - 1]);
        }
        MPI_Waitall(nodes_cnt - 1, requests, nullptr);

        for (auto res : nodes_results) {
            total_integrall_value += res;
        }
    }

    return total_integrall_value;
}

void save_to_file(int curr_rank, double value) {
    if (curr_rank != 0) {
        return;
    }

    try {
        std::ofstream file("output.txt");
        file << std::fixed << std::setprecision(10) << value << std::endl;
        file.close();
    } catch (std::exception& e) {
        std::cout << "[FATAL] " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int curr_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &curr_rank);

    int nodes_cnt;
    MPI_Comm_size(MPI_COMM_WORLD, &nodes_cnt);

    std::cout << "[MPI_INIT_INFO] Current node rank: " << curr_rank << "\tTotal nodes count: " << nodes_cnt << std::endl;

    auto params = get_parameters(curr_rank);
    if (!validate_input(params)) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    auto start_time = std::chrono::steady_clock::now();
    auto result = calculate_integral(curr_rank, nodes_cnt, params);
    auto finish_time = std::chrono::steady_clock::now();

    save_to_file(curr_rank, result);

    MPI_Finalize();

    if (curr_rank == 0) {
        std::cout << "[HOST - 0] Integral calculation finished successfully: " << std::fixed << std::setprecision(10) << result << std::endl;  

        std::cout << std::endl;
        std::cout << "[HOST - 0] Execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time).count() << "[ms]" << std::endl;
        std::cout << "[HOST - 0] Execution time: " << std::chrono::duration_cast<std::chrono::microseconds>(finish_time - start_time).count() << "[µs]" << std::endl;
        std::cout << "[HOST - 0] Execution time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish_time - start_time).count() << "[ns]" << std::endl;
    }

    return 0;
}
