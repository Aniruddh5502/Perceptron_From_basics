#include <iostream>
#include <vector>
#include <random>
#include "matrix_op.h"

// =====================================================================================================
// MAKING THE DATASET FOR THE LOGISTIC REGRESSION MODEL WITH PRECISE CALCULATION TO APPROXIMATE        ||
// AND ADDING UNCERTAINTIES FOR RANDOMIZATION AND USING THE DATA WITH LITTLE AMOUNT OF DISTURBANCES    ||
// =====================================================================================================

int main(){
    std::vector<std::vector<float>> data = random_matrix(20000, 5, -1.00, 1.00);
    print_n_row(data, 6);


    // generating random input(first 5 columns) and then calculating 
    // the y_true (last column) from the   GOVERNING FORMULA 
    // ===================================================================
    //                                                                   ||
    // x_5  =  x_1*1.2131 + x_2*2.3313 + x_3*x_3*5.673 + x_4*3.1232;     ||
    //                                                                   ||
    // ===================================================================

    for(int i = 0; i<data.size(); i++){
        data[i][4] = data[i][0]*1.2132 + data[i][1]*2.3313 + data[i][2]*data[i][2]*5.673 + data[i][3]*3.1232;
    }

    write_matrix_to_CSV(data, "dataset_Perceptron");
    print_n_row(data,6);
    std::cout << "Data generation Done.......\n\n" << std::endl;



    return 0;
}
