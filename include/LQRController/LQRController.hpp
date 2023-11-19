#ifndef LQRCONTROLLER_HPP
#define LQRCONTROLLER_HPP

#include <armadillo>
class LQRController {
public:
    LQRController(const arma::mat& A,
                  const arma::mat& B,
                  const arma::mat& Q,
                  const arma::mat& R);

    arma::mat computeOptimalGain();

private:
    arma::mat A_; // System matrix
    arma::mat B_; // Input matrix
    arma::mat Q_; // State cost matrix
    arma::mat R_; // Control cost matrix

    arma::mat solveDARE(const arma::mat& A,
                        const arma::mat& B,
                        const arma::mat& Q,
                        const arma::mat& R);
};

#endif