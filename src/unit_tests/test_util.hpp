#ifndef TARQUIN__TEST_UTIL__
#define TARQUIN__TEST_UTIL__

#include "common.hpp"
#include <Eigen/Dense>

cvm::cvector read_ref_data(std::string filename);
void write_ref_data(cvm::cvector cvec, std::string filename);
Eigen::VectorXcd read_ref_data_eigen(std::string filename);
void write_ref_data_eigen(Eigen::VectorXcd cvec, std::string filename);

#endif
