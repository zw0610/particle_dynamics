#ifndef kernel_helper_cuh
#define kernel_helper_cuh

#include <iostream>

template<class T>
void print_dev_array(T * list, size_t length, size_t num_col) {

    T * ha = new T[length];
    cudaMemcpy( ha , list, length*sizeof(T) , cudaMemcpyDeviceToHost);

    for (size_t i = 0; i<length; i++) {
        std::cout << ha[i] << " ";
        if ( i%num_col==0 ) {std::cout << std::endl;};
    }

    delete [] ha;
}

#endif
