#include <QuaternionUtils/quaternion_op.hpp>

// Hat function for skew-symmetric matrix
arma::mat hat(const arma::vec &v) {
    arma::mat hat_matrix = arma::zeros<arma::mat>(3,3);
    hat_matrix(0,1) = -v(2);
    hat_matrix(0,2) = v(1);
    hat_matrix(1,0) = v(2);
    hat_matrix(1,2) = -v(0);
    hat_matrix(2,0) = -v(1);
    hat_matrix(2,1) = v(0);
    return hat_matrix;
}

// L function for quaternion
arma::mat L(const arma::vec &q) {
    double s = q(0);
    arma::vec v = q.subvec(1, 3);
    arma::mat L = arma::join_vert(arma::join_horiz(arma::vec({s}), -v.t()),
                                    arma::join_horiz(v, s*arma::eye<arma::mat>(3,3) + hat(v)));
    return L;
}

// Quaternion to rotation matrix
arma::mat qtoQ(const arma::vec &q) {
    arma::mat T = arma::diagmat(arma::vec({1, -1, -1, -1}));
    arma::mat H = arma::join_vert(arma::zeros<arma::rowvec>(3), arma::eye<arma::mat>(3,3));
    return H.t() * T * L(q) * T * L(q) * H;
}

// G function for quaternion
arma::mat G(const arma::vec &q) {
    return L(q) * arma::join_vert(arma::zeros<arma::rowvec>(3), arma::eye<arma::mat>(3,3));
}

// rotation vector to quaternion
arma::vec rptoq(const arma::vec &phi) {
    return (1 / sqrt(1 + arma::dot(phi, phi))) * arma::join_vert(arma::vec({1}), phi);
}

// quaternion to rotation vector
arma::vec qtorp(const arma::vec &q) {
    return q.subvec(1, 3) / q(0);
}

// E function for quaternion
arma::mat E(const arma::vec &q) {
    arma::mat Gq = G(q);
    arma::mat E = arma::zeros<arma::mat>(13, 12);
    E(arma::span(0, 2), arma::span(0, 2)) = arma::eye<arma::mat>(3, 3);
    E(arma::span(3, 6), arma::span(3, 5)) = Gq;
    E(arma::span(7, 12), arma::span(6, 11)) = arma::eye<arma::mat>(6, 6);
    return E;
}
