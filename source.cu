#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>


#include <cublas_v2.h>
#include <curand.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "basic_operations.h"
#include "matrix_operations.h"
#include "defs01.h"
#include "tests01.h"




int main() {
	tester_01();
	return 0;
}







