#ifndef MATRIX_OP_H
#define MATRIX_OP_H



// =================================================
// UTILITY  AND MATRIX OPERATION FUNCTIONS        ||
// =================================================

std::vector<std::vector<float>> make_matrix(int rows, int cols, float value);

void print_matrix(std::vector<std::vector<float>>& mat, int precision = 4, int width = 10);

void print_n_row(std::vector<std::vector<float>>& mat, int n, int precision = 4, int width = 10);


void num_print(float num, int precision = 4, int width = 10);

void print_vector(std::vector<float>& vec, int precision = 4, int width = 10);

std::vector<std::vector<float>> csv_to_matrix(const std::string& filename, char delimiter = ',');

bool write_matrix_to_CSV(const std::vector<std::vector<float>>& matrix, const std::string& filename, char delimiter = ',');

std::vector<std::vector<float>> normalize_matrix(std::vector<std::vector<float>>& matrix);

std::vector<std::vector<float>> add_mat(std::vector<std::vector<float>>& mat_1, std::vector<std::vector<float>>& mat_2);

std::vector<std::vector<float>> subtract_mat(std::vector<std::vector<float>> &mat_1, std::vector<std::vector<float>> &mat_2);

std::vector<std::vector<float>> matrix_transpose(std::vector<std::vector<float>>& mat);

std::vector<std::vector<float>> matrix_multiply(const std::vector<std::vector<float>>& a, const std::vector<std::vector<float>>& b);

std::vector<std::vector<float>> matrix_scaler_multiply(const std::vector<std::vector<float>> mat, float value);

float dot_product(std::vector<float>& vec_1, std::vector<float>& vec_2);

std::vector<std::vector<float>> generate_linear_regression_dataset(int num_samples, float m, float c, float noise_level);

std::vector<float> add_vec(std::vector<float>& vec_1, std::vector<float>& vec_2);

std::vector<float> subtract_vec(std::vector<float>& vec_1, std::vector<float>& vec_2);

std::vector<float> vector_scaler_multiply(std::vector<float>& vec, float m);

std::vector<float> vec_mul(std::vector<float>& vec_1, std::vector<float>& vec_2);

float vec_sum(std::vector<float>& vec);

float mat_sum(std::vector<std::vector<float>>& mat);

std::vector<std::vector<float>> random_matrix(int rows, int cols, float upper, float lower);

void mat_dim(std::vector<std::vector<float>>& mat);

std::vector<std::vector<float>> broad_cast(const std::vector<std::vector<float>>& a, const std::vector<std::vector<float>>& b);

float relu(float x);

std::vector<std::vector<float>> ReLU(const std::vector<std::vector<float>> &matrix);

std::vector<std::vector<float>> mat_element_square(std::vector<std::vector<float>> &mat);

std::vector<std::vector<float>> sum_rows(const std::vector<std::vector<float>>& matrix);

std::vector<std::vector<float>> elementwise_multiply(const std::vector<std::vector<float>>& a, const std::vector<std::vector<float>>& b);


float relu_derivative_single(float x);

std::vector<std::vector<float>> relu_derivative(const std::vector<std::vector<float>>& z);


#endif