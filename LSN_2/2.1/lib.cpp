#include "lib.h"

double Error(double averages, double squared, int n){

    if (n ==0){
        return 0;
    }

    else{
        return sqrt((squared - pow(averages, 2)) / (double) n);
    }
}
