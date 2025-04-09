#include <iostream>
#include <vector>
#include <cmath>
#include "matrix_op.h"

int main()
{
    std::vector<std::vector<float>> data_train = csv_to_matrix("data_train");
    std::vector<std::vector<float>> data_valid = csv_to_matrix("data_valid");

    // separating the input and output
    std::vector<std::vector<float>> x_train = make_matrix(18000, 4, 0.00);
    std::vector<std::vector<float>> y_train = make_matrix(18000, 1, 0.00);
    for (int i = 0; i < data_train.size(); i++)
    {
        x_train[i][0] = data_train[i][0];
        x_train[i][1] = data_train[i][1];
        x_train[i][2] = data_train[i][2];
        x_train[i][3] = data_train[i][3];

        y_train[i][0] = data_train[i][4];
    }

    std::vector<std::vector<float>> x_valid = make_matrix(2000, 4, 0.00);
    std::vector<std::vector<float>> y_valid = make_matrix(2000, 4, 0.00);

    for (int i = 0; i < data_valid.size(); i++)
    {
        x_valid[i][0] = data_valid[i][0];
        x_valid[i][1] = data_valid[i][1];
        x_valid[i][2] = data_valid[i][2];
        x_valid[i][3] = data_valid[i][3];

        y_valid[i][0] = data_valid[i][4];
    }

    // ====================================================================================================
    // Making it a 2 layer 10, 10 node Neural Network with no activation function so that it can predict  ||
    // the exact value which is no bound with any range                                                   ||                        
    // ====================================================================================================

    std::vector<std::vector<float>> w1 = csv_to_matrix("w1.1");
    std::vector<std::vector<float>> w2 = csv_to_matrix("w2.1");
    std::vector<std::vector<float>> b1 = csv_to_matrix("b1.1");
    std::vector<std::vector<float>> b2 = csv_to_matrix("b2.1");

    int batch_size = 2000;
    int avg_valid = (mat_sum(y_valid) / 2000);

    // predictions
    std::vector<std::vector<float>> z_ = matrix_multiply(x_valid, w1);
    std::vector<std::vector<float>> z1 = broad_cast(z_, b1);
    std::vector<std::vector<float>> a1 = ReLU(z1);
    std::vector<std::vector<float>> z2 = matrix_multiply(a1, w2);
    std::vector<std::vector<float>> a2 = broad_cast(z2, b2);

    // Calculate validation error (Mean Squared Error)
    std::vector<std::vector<float>> error = subtract_mat(a2, y_valid);
    std::vector<std::vector<float>> err_2 = mat_element_square(error);
    float validation_loss = mat_sum(err_2) / (2 * 2000.0f); // Corrected denominator

    std::cout << "Validation Loss (MSE): " << validation_loss << std::endl;

    return 0;
}
