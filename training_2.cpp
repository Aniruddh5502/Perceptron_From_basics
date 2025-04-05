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

    // ====================================================================================================
    // Making it a 2 layer 10, 10 node Neural Network with no activation function so that it can predict  ||
    // the exact value which is no bound with any range                                                   ||                        ||
    // ====================================================================================================

    std::vector<std::vector<float>> w1 = csv_to_matrix("w1");
    std::vector<std::vector<float>> w2 = csv_to_matrix("w2");
    std::vector<std::vector<float>> b1 = csv_to_matrix("b1");
    std::vector<std::vector<float>> b2 = csv_to_matrix("b2");

    float learning_rate = 0.002;
    int epochs = 10000;
    int batch_size = 18000;

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        // predictions
        std::vector<std::vector<float>> z_ = matrix_multiply(x_train, w1);
        std::vector<std::vector<float>> z1 = broad_cast(z_, b1);
        std::vector<std::vector<float>> a1 = ReLU(z1);
        std::vector<std::vector<float>> z2 = matrix_multiply(a1, w2);
        std::vector<std::vector<float>> a2 = broad_cast(z2, b2);

        // Error calculation
        std::vector<std::vector<float>> error = subtract_mat(a2, y_train);
        std::vector<std::vector<float>> err_2 = mat_element_square(error); // element square
        float ERR = mat_sum(err_2) / (2 * 18000);

        // ===============================
        // Backpropagation              ||
        // ===============================

        // Output layer gradients
        std::vector<std::vector<float>> dz2 = matrix_scaler_multiply(error, 1.0f / 18000);
        std::vector<std::vector<float>> a1t = matrix_transpose(a1);
        std::vector<std::vector<float>> dw2 = matrix_multiply(a1t, dz2);
        std::vector<std::vector<float>> db2 = sum_rows(dz2);

        // Hidden layer gradients
        std::vector<std::vector<float>> dz1 = elementwise_multiply(matrix_multiply(dz2, matrix_transpose(w2)), relu_derivative(z1));
        std::vector<std::vector<float>> dw1 = matrix_multiply(matrix_transpose(x_train), dz1);
        std::vector<std::vector<float>> db1 = sum_rows(dz1);

        // =======================================
        // UPDATE PARAMETERS
        // =======================================
        std::vector<std::vector<float>> msm = matrix_scaler_multiply(dw2, learning_rate);
        w2 = subtract_mat(w2, msm);
        std::vector<std::vector<float>> ms2 = matrix_scaler_multiply(db2, learning_rate);
        b2 = subtract_mat(b2, ms2);
        std::vector<std::vector<float>> ms3 = matrix_scaler_multiply(dw1, learning_rate);
        w1 = subtract_mat(w1, ms3);
        std::vector<std::vector<float>> ms4 = matrix_scaler_multiply(db1, learning_rate);
        b1 = subtract_mat(b1, ms4);

        // Print progress
        if (epoch % 100 == 0)
        {
            std::cout << "Epoch "        << epoch << " |            Loss: " << ERR << std::endl;
        }
    }

    //=================================================================
    //  AFTER TRAINING SAVE THE WEIGHTS AND BIASES
    //=================================================================

    write_matrix_to_CSV(w1,"w1.1");
    write_matrix_to_CSV(w2,"w2.1");
    write_matrix_to_CSV(b1,"b1.1");
    write_matrix_to_CSV(b2,"b2.1");

    return 0;
}