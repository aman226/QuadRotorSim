#ifndef _QUATERNION_OP_HPP
#define _QUATERNION_OP_HPP

#include <armadillo>

// Hat function for skew-symmetric matrix
arma::mat hat(const arma::vec &v);

// L function for quaternion
arma::mat L(const arma::vec &q);

// Quaternion to rotation matrix
arma::mat qtoQ(const arma::vec &q);

// G function for quaternion
arma::mat G(const arma::vec &q);

// rotation vector to quaternion
arma::vec rptoq(const arma::vec &phi);

// quaternion to rotation vector
arma::vec qtorp(const arma::vec &q);

// E function for quaternion
arma::mat E(const arma::vec &q);

#endif