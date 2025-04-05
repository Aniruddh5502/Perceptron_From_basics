#include <iostream>
#include <vector>
#include "matrix_op.h"

int main(){
    // data_logistic_regression has 5 columns and 20,000 rows
    std::vector<std::vector<float>> data = csv_to_matrix("dataset_Perceptron");
    print_n_row(data, 6);

    /*
    separating dataset for training and testing
    */

    std::vector<std::vector<float>> data_train = make_matrix(18000, 5, 0.00);
    std::vector<std::vector<float>> data_valid = make_matrix( 2000, 5, 0.00);


    
    for(int i = 0; i<18000; i++){
        data_train[i][0] = data[i][0];
        data_train[i][1] = data[i][1];
        data_train[i][2] = data[i][2];
        data_train[i][3] = data[i][3];
        data_train[i][4] = data[i][4];
    }
    for(int i = 18000; i<20000; i++){
        data_valid[i-18000][0] = data[i][0];
        data_valid[i-18000][1] = data[i][1];
        data_valid[i-18000][2] = data[i][2];
        data_valid[i-18000][3] = data[i][3];
        data_valid[i-18000][4] = data[i][4];
    }

    
    std::cout << "Data Train: " << std::endl;
    print_n_row(data_train, 5);
    std::cout << "Data train size: " << data_train.size() << " rows." <<"\n\n" <<std::endl;
    std::cout << "Data Valid: " << std::endl;
    print_n_row(data_valid, 5);
    std::cout << "Data Valid size: " << data_valid.size() << " rows." <<"\n\n" <<std::endl;

    /*
    Saving the training and validation data as a csv file
    */
    write_matrix_to_CSV(data_train, "data_train");
    write_matrix_to_CSV(data_valid, "data_valid");
    std::cout << "Data Saved.. For training and validation.." << std::endl;

    return 0;
}