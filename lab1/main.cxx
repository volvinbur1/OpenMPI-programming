#include <mpi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <cmath>
#include <iomanip>
#include <chrono>

const double accuracy = 1e-8;
const int x_value_tag = 1;
const int seq_number_tag = 2;
const int series_element_value_tag = 3;
const int break_flag_tag = 4;

double read_x_from_file() {
        try {
                std::ifstream file("input.txt");
                std::stringstream input_str_buffer;
                input_str_buffer << file.rdbuf();

                return std::stod(input_str_buffer.str());
        } catch (std::exception& e) {
                std::cout << "[FATAL] " << e.what() << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
        }
        return 0;
}

double get_x_value(int curr_rank, int nodes_cnt) {
        double input_x_value;
        if (curr_rank == 0) {
                input_x_value = read_x_from_file();
                std::cout << "[HOST - 0] square root param value: " << std::fixed << std::setprecision(10) << input_x_value << std::endl;
                std::cout << "[HOST - 0] accuracy: " << std::defaultfloat << accuracy << std::endl;

                input_x_value -= 1;

                for (auto dest_rank = 1; dest_rank < nodes_cnt; dest_rank++) {
                        MPI_Send(&input_x_value, 1, MPI_DOUBLE, dest_rank, x_value_tag, MPI_COMM_WORLD);
                }
        } else {
                MPI_Recv(&input_x_value, 1, MPI_DOUBLE, 0, x_value_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        return input_x_value;
}

bool validate_input(double x_value) {
        if (fabs(x_value) > 1) {
                std::cout << "[ERROR] Input param value -1 is out of the range [-1, 1]. Exit..." << std::endl;
                return false;
        }

        return true;
}

double factorial(int value) {
        if (value < 0) {
                return 0;
        } else if (value == 0) {
                return 1;
        }
        return value * factorial(value - 1);
}

double calc_series_element(double x, int seq_number) {
        auto numerator = pow(-1, seq_number - 1) * factorial(2*seq_number);
        auto denominator = pow(4, seq_number) * pow(factorial(seq_number), 2) * (2 * seq_number - 1);
        return numerator / denominator * pow(x, seq_number);
}

int get_target_seq_number(int curr_rank, int nodes_cnt, int* last_seq_number) {
        int target_seq_number;
        if (curr_rank == 0) {
                target_seq_number = (*last_seq_number)++;
                for (auto dest_rank = 1; dest_rank < nodes_cnt; dest_rank++) {
                        MPI_Send(last_seq_number, 1, MPI_INT, dest_rank, seq_number_tag, MPI_COMM_WORLD);
                        (*last_seq_number)++;
                }
        } else {
                MPI_Recv(&target_seq_number, 1, MPI_INT, 0, seq_number_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        return target_seq_number;
}

bool get_break_flag(int curr_rank, int nodes_cnt, double* series_elements_sum, double series_element_value) {
        auto break_calculation = false;
        if (curr_rank == 0) {
                (*series_elements_sum) += series_element_value;
                if (fabs(series_element_value) < accuracy) {
                        break_calculation = true;
                        std::cout << "[HOST - 0] break" << std::endl;
                } else {
                        for (auto dest_rank = 1; dest_rank < nodes_cnt; dest_rank++) {
                                MPI_Recv(&series_element_value, 1, MPI_DOUBLE, dest_rank, series_element_value_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                (*series_elements_sum) += series_element_value;
                                if (fabs(series_element_value) < accuracy) {
                                        break_calculation = true;
                                        std::cout << "[NODE - " << dest_rank << "] break" << std::endl;
                                        break;
                                }
                        }
                }

                for (auto dest_rank = 1; dest_rank < nodes_cnt; dest_rank++) {
                        MPI_Send(&break_calculation, 1, MPI_BYTE, dest_rank, break_flag_tag, MPI_COMM_WORLD);
                }
        } else {
                MPI_Send(&series_element_value, 1, MPI_DOUBLE, 0, series_element_value_tag, MPI_COMM_WORLD);
                MPI_Recv(&break_calculation, 1, MPI_BYTE, 0, break_flag_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        return break_calculation;
}

double calculate_square_root(int curr_rank, int nodes_cnt, double x) {
        if (curr_rank == 0) {
                std::cout << "[HOST - 0] square root calculation started." << std::endl;
        }

        auto last_seq_number = 0;
        auto series_elements_sum = 0.0;

        while (true) {
                auto target_seq_number = get_target_seq_number(curr_rank, nodes_cnt, &last_seq_number);
                auto series_element_value = calc_series_element(x, target_seq_number);
                
                // std::cout << "[CalcNode - "<< curr_rank << "][SeqNumber - " << target_seq_number 
                //         << "] calculate_square_root::series_element_value: " << series_element_value << std::endl;
                
                if (get_break_flag(curr_rank, nodes_cnt, &series_elements_sum, series_element_value)) {
                        break;
                }
        }

        return series_elements_sum;
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

        auto sqrt_param = get_x_value(curr_rank, nodes_cnt);
        if (!validate_input(sqrt_param)) {
                MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        auto start_time = std::chrono::steady_clock::now();
        auto result = calculate_square_root(curr_rank, nodes_cnt, sqrt_param);
        auto finish_time = std::chrono::steady_clock::now();

        save_to_file(curr_rank, result);

        MPI_Finalize();

        if (curr_rank == 0) {
                std::cout << "[HOST - 0] Calculation finished. Square root of '" <<
                        std::fixed << std::setprecision(10) << 1 + sqrt_param << "' is equal " << result << std::endl;

                std::cout << std::endl;
                std::cout << "[HOST - 0] Execution time: " << std::chrono::duration_cast<std::chrono::microseconds>(finish_time - start_time).count() << "[µs]" << std::endl;
                std::cout << "[HOST - 0] Execution time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish_time - start_time).count() << "[ns]" << std::endl;
        }

        return 0;
}
