#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <random>
#include "matrix_op.h"

std::vector<std::vector<float>> make_matrix(int rows, int cols, float value)
{
    return std::vector<std::vector<float>>(rows, std::vector<float>(cols, value));
}

void print_matrix(std::vector<std::vector<float>> &mat, int precision, int width)
{
    if (mat.empty() || mat[0].empty())
    {
        std::cout << "(empty matrix)\n";
        return;
    }
    std::cout << std::fixed << std::setprecision(precision);
    for (const auto &row : mat)
    {
        for (float val : row)
        {
            std::cout << std::setw(width) << val << " ";
        }
        std::cout << std::endl;
        ;
    }
    std::cout << "\n"
              << std::defaultfloat;
}

void print_n_row(std::vector<std::vector<float>> &mat, int n, int precision, int width)
{
    if (mat.empty() || mat[0].empty())
    {
        std::cout << "(empty matrix)\n";
        return;
    }

    std::cout << std::fixed << std::setprecision(precision);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < mat[0].size(); j++)
        {
            std::cout << std::setw(width) << mat[i][j] << " ";
        }
        std::cout << std::endl;
        ;
    }
    std::cout << "...." << std::endl;
    std::cout << "....." << std::endl;
    std::cout << "\n"
              << std::defaultfloat;
}

void num_print(float num, int precision, int width)
{
    std::cout << std::fixed << std::setprecision(precision);
    std::cout << std::setw(width) << num << std::endl;
    std::cout << "\n"
              << std::endl;
}

void print_vector(std::vector<float> &vec, int precision, int width)
{
    int size = vec.size();

    std::cout << std::fixed << std::setprecision(precision);
    for (int i = 0; i < size; i++)
    {
        std::cout << std::setw(width) << vec[i] << " ";
    }
    std::cout << "\n"
              << std::endl;
}

void mat_dim(std::vector<std::vector<float>> &mat)
{
    int rows = mat.size();
    int cols = mat[0].size();
    std::cout << rows << "x" << cols << std::endl;
}

std::vector<std::vector<float>> csv_to_matrix(const std::string &filename, char delimiter)
{

    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::vector<std::vector<float>> matrix;
    std::string line;

    while (std::getline(file, line))
    {
        std::vector<float> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, delimiter))
        {
            try
            {
                row.push_back(std::stof(cell));
            }
            catch (const std::invalid_argument &)
            {
                throw std::runtime_error("Invalid number in CSV: " + cell);
            }
        }

        // Skip empty rows
        if (!row.empty())
        {
            // Validate consistent column count
            if (!matrix.empty() && row.size() != matrix[0].size())
            {
                throw std::runtime_error("Inconsistent column count in CSV");
            }
            matrix.push_back(row);
        }
    }

    if (matrix.empty())
    {
        throw std::runtime_error("CSV file is empty or contains no valid data");
    }

    return matrix;
}

bool write_matrix_to_CSV(const std::vector<std::vector<float>> &matrix, const std::string &filename, char delimiter)
{
    // Open the output file
    std::ofstream outputFile(filename);

    // Check if file opened successfully
    if (!outputFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    // Write each row of the matrix
    for (const auto &row : matrix)
    {
        for (size_t i = 0; i < row.size(); ++i)
        {
            outputFile << row[i];

            // Add delimiter unless it's the last element
            if (i != row.size() - 1)
            {
                outputFile << delimiter;
            }
        }
        outputFile << "\n";
    }

    // Close the file
    outputFile.close();
    return true;
}

std::vector<std::vector<float>> normalize_matrix(std::vector<std::vector<float>> &matrix)
{
    if (matrix.empty())
        return {};

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    // Initialize min and max values for each column
    std::vector<float> min_values(cols, std::numeric_limits<float>::max());
    std::vector<float> max_values(cols, std::numeric_limits<float>::lowest());

    // Find min and max for column
    for (size_t j = 0; j < cols; ++j)
    {
        for (size_t i = 0; i < rows; ++i)
        {
            min_values[j] = std::min(min_values[j], matrix[i][j]);
            max_values[j] = std::max(max_values[j], matrix[i][j]);
        }
    }

    // Normalize the matrix
    std::vector<std::vector<float>> normalized(rows, std::vector<float>(cols));
    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            if (max_values[j] != min_values[j])
            { // Avoid division by zero
                normalized[i][j] = (matrix[i][j] - min_values[j]) / (max_values[j] - min_values[j]);
            }
            else
            {
                normalized[i][j] = 0.000; // if all values are same
            }
        }
    }
    return normalized;
}

std::vector<std::vector<float>> add_mat(std::vector<std::vector<float>> &mat_1, std::vector<std::vector<float>> &mat_2)
{

    std::vector<std::vector<float>> sum(mat_1.size(), std::vector<float>(mat_1[0].size()));

    if ((mat_1.size() != mat_2.size()) && (mat_1[0].size() != mat_2[0].size()))
    {
        std::cout << "Matrix Size Don't Match.. ERROR" << std::endl;
    }
    else
    {
        for (int i = 0; i < mat_1.size(); i++)
        {
            for (int j = 0; j < mat_1[0].size(); j++)
            {
                sum[i][j] = mat_1[i][j] + mat_2[i][j];
            }
        }
    }
    return sum;
}

std::vector<std::vector<float>> subtract_mat(std::vector<std::vector<float>> &mat_1, std::vector<std::vector<float>> &mat_2)
{
    std::vector<std::vector<float>> dif(mat_1.size(), std::vector<float>(mat_1[0].size()));

    if ((mat_1.size() != mat_2.size()) && (mat_1[0].size() != mat_2[0].size()))
    {
        std::cout << "Matrix Size Don't Match.. ERROR" << std::endl;
    }
    else
    {
        for (int i = 0; i < mat_1.size(); i++)
        {
            for (int j = 0; j < mat_1[0].size(); j++)
            {
                dif[i][j] = mat_1[i][j] - mat_2[i][j];
            }
        }
    }
    return dif;
}

std::vector<std::vector<float>> matrix_transpose(std::vector<std::vector<float>> &mat)
{
    std::vector<std::vector<float>> transposed(mat[0].size(), std::vector<float>(mat.size()));
    if (mat.empty())
    {
        std::cout << "ERROR... Empty matrix..." << std::endl;
    }
    else
    {
        for (int i = 0; i < mat.size(); i++)
        {
            for (int j = 0; j < mat[0].size(); j++)
            {
                transposed[j][i] = mat[i][j];
            }
        }
    }
    return transposed;
}

std::vector<std::vector<float>> matrix_multiply(const std::vector<std::vector<float>> &a, const std::vector<std::vector<float>> &b)
{
    // Check dimensions
    if (a.empty() || b.empty() || a[0].size() != b.size())
    {
        std::cerr << "Error: Invalid matrix dimensions for multiplication\n";
        return {};
    }

    size_t a_rows = a.size();
    size_t a_cols = a[0].size();
    size_t b_cols = b[0].size();

    std::vector<std::vector<float>> result(a_rows, std::vector<float>(b_cols, 0.0f));

    for (size_t i = 0; i < a_rows; ++i)
    {
        for (size_t j = 0; j < b_cols; ++j)
        {
            for (size_t k = 0; k < a_cols; ++k)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return result;
}

std::vector<std::vector<float>> matrix_scaler_multiply(const std::vector<std::vector<float>> mat, float value)
{
    std::vector<std::vector<float>> new_matrix(mat.size(), std::vector<float>(mat[0].size()));
    for (int i = 0; i < mat.size(); i++)
    {
        for (int j = 0; j < mat[0].size(); j++)
        {
            new_matrix[i][j] = value * mat[i][j];
        }
    }
    return new_matrix;
}

float dot_product(std::vector<float> &vec_1, std::vector<float> &vec_2)
{
    if (vec_1.size() != vec_2.size())
    {
        std::cerr << "Error: Vector dimensions don't match (" << vec_1.size() << " vs " << vec_2.size() << ")\n";
        return 0.0f; // Or throw an exception
    }

    float result = 0.0f;
    for (size_t i = 0; i < vec_1.size(); ++i)
    {
        result += vec_1[i] * vec_2[i]; // Accumulate the sum
    }
    return result;
}

std::vector<std::vector<float>> generate_linear_regression_dataset(int num_samples, float m, float c, float noise_level)
{
    std::vector<std::vector<float>> dataset;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> x_dist(-10, 10);      // X values between -10 and 10
    std::normal_distribution<float> noise_dist(0, noise_level); // Gaussian noise

    for (int i = 0; i < num_samples; i++)
    {
        float x = x_dist(gen);
        float noise = noise_dist(gen);
        float y = m * x + c + noise;
        dataset.push_back({x, y});
    }

    return dataset;
}

std::vector<float> add_vec(std::vector<float> &vec_1, std::vector<float> &vec_2)
{
    std::vector<float> sum(vec_1.size());
    if ((vec_1.size() == vec_2.size()))
    {
        for (int i = 0; i < vec_1.size(); i++)
        {
            sum[i] = vec_1[i] + vec_2[i];
        }
    }
    else
    {
        std::cout << "ERROR--->>>Dimension Mismatch." << std::endl;
    }
    return sum;
}

std::vector<float> subtract_vec(std::vector<float> &vec_1, std::vector<float> &vec_2)
{
    std::vector<float> sub(vec_1.size());
    if ((vec_1.size() == vec_2.size()))
    {
        for (int i = 0; i < vec_1.size(); i++)
        {
            sub[i] = vec_1[i] - vec_2[i];
        }
    }
    else
    {
        std::cout << "ERROR--->>>Dimension Mismatch." << std::endl;
    }
    return sub;
}

std::vector<float> vector_scaler_multiply(std::vector<float> &vec, float m)
{
    std::vector<float> vec_mul(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        vec_mul[i] = vec[i] * m;
    }
    return vec_mul;
}

std::vector<float> vec_mul(std::vector<float> &vec_1, std::vector<float> &vec_2)
{
    std::vector<float> mul(vec_1.size());
    if (vec_1.size() == vec_2.size())
    {
        for (int i = 0; i < vec_1.size(); i++)
        {
            mul[i] = vec_1[i] * vec_2[i];
        }
    }
    else
    {
        std::cout << "ERROR--->>>Dimension Mismatch." << std::endl;
    }
    return mul;
}

float vec_sum(std::vector<float> &vec)
{
    float sum = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        sum += vec[i];
    }
    return sum;
}

float mat_sum(std::vector<std::vector<float>> &mat)
{
    float sum = 0;
    for (int i = 0; i < mat.size(); i++)
    {
        for (int j = 0; j < mat[0].size(); j++)
        {
            sum += mat[i][j];
        }
    }
    return sum;
}

std::vector<std::vector<float>> random_matrix(int rows, int cols, float upper, float lower)
{
    std::vector<std::vector<float>> rand(rows, std::vector<float>(cols, 1));
    std::random_device rd;  // seed
    std::mt19937 gen(rd()); // Mersenne twister Engine
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            std::uniform_real_distribution<float> dis(lower, upper);
            rand[i][j] = dis(gen);
        }
    }
    return rand;
}

std::vector<std::vector<float>> broad_cast(const std::vector<std::vector<float>> &a, const std::vector<std::vector<float>> &b)
{
    // Check for empty matrices
    if (a.empty() || b.empty() || a[0].empty() || b[0].empty())
    {
        throw std::invalid_argument("Matrices must not be empty");
    }

    const size_t a_rows = a.size();
    const size_t a_cols = a[0].size();
    const size_t b_rows = b.size();
    const size_t b_cols = b[0].size();

    // Determine output dimensions
    const size_t out_rows = std::max(a_rows, b_rows);
    const size_t out_cols = std::max(a_cols, b_cols);

    // Create result matrix
    std::vector<std::vector<float>> result(out_rows, std::vector<float>(out_cols, 0));

    for (size_t i = 0; i < out_rows; ++i)
    {
        for (size_t j = 0; j < out_cols; ++j)
        {
            // Calculate indices with broadcasting
            const size_t a_i = (a_rows == 1) ? 0 : i;
            const size_t a_j = (a_cols == 1) ? 0 : j;
            const size_t b_i = (b_rows == 1) ? 0 : i;
            const size_t b_j = (b_cols == 1) ? 0 : j;

            // Perform bounds checking
            if (a_i >= a_rows || a_j >= a_cols || b_i >= b_rows || b_j >= b_cols)
            {
                throw std::out_of_range("Matrix index out of bounds");
            }

            result[i][j] = a[a_i][a_j] + b[b_i][b_j];
        }
    }
    return result;
}

float relu(float x) {
    return std::max(0.0f, x);
}

std::vector<std::vector<float>> ReLU(const std::vector<std::vector<float>> &matrix)
{
    std::vector<std::vector<float>> result = matrix; // Copy the input matrix

    for (auto &row : result)
    {
        for (auto &element : row)
        {
            element = relu(element); // Apply ReLU to each element
        }
    }

    return result;
}

std::vector<std::vector<float>> mat_element_square(std::vector<std::vector<float>> &mat){
    std::vector<std::vector<float>> square(mat.size(), std::vector<float>(mat[0].size(),0));
    for(int i = 0; i<mat.size(); i++){
        for(int j = 0; j<mat[0].size(); j++){
            square[i][j] = mat[i][j] * mat[i][j];
        }
    }
    return square;
}

std::vector<std::vector<float>> sum_rows(const std::vector<std::vector<float>>& matrix) {
    if(matrix.empty()) return {};
    
    // Create a 1xN matrix to store column sums
    std::vector<std::vector<float>> result(1, std::vector<float>(matrix[0].size(), 0.0f));
    
    for(const auto& row : matrix) {         // For each training example
        for(size_t j = 0; j < row.size(); j++) { // For each neuron/column
            result[0][j] += row[j];         // Accumulate column values
        }
    }
    
    return result;
}

std::vector<std::vector<float>> elementwise_multiply(const std::vector<std::vector<float>>& a,  const std::vector<std::vector<float>>& b) {
if(a.size() != b.size() || a[0].size() != b[0].size()) {
throw std::invalid_argument("Matrices must have same dimensions");
}

std::vector<std::vector<float>> result(a.size(), std::vector<float>(a[0].size()));
for(size_t i = 0; i < a.size(); ++i) {
for(size_t j = 0; j < a[0].size(); ++j) {
result[i][j] = a[i][j] * b[i][j];
}
}
return result;
}

float relu_derivative_single(float x) {
    return (x > 0.0f) ? 1.0f : 0.0f;
}

std::vector<std::vector<float>> relu_derivative(const std::vector<std::vector<float>>& z) {
    std::vector<std::vector<float>> result(z.size(), std::vector<float>(z[0].size()));
    
    for (size_t i = 0; i < z.size(); ++i) {
        for (size_t j = 0; j < z[0].size(); ++j) {
            result[i][j] = relu_derivative_single(z[i][j]);
        }
    }
    
    return result;
}



