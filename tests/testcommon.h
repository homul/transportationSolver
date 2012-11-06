#ifndef TESTCOMMON_H
#define TESTCOMMON_H
#include "common.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace ublas = boost::numeric::ublas;
typedef ublas::compressed_vector<DD::floatType> SparseArray;
typedef ublas::matrix<DD::floatType> Matrix;
#endif
