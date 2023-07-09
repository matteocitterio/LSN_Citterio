#include "lib.h"

double Error(double averages, double squared, int n){

    /*
    This computes the blocking average error
    */

    if (n ==0){                                                     // If it is the first block return 0
        return 0;
    }

    else{
        return sqrt((squared - pow(averages, 2)) / (double) n);
    }
}
