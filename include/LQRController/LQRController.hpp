#ifndef LQRCONTROLLER_HPP
#define LQRCONTROLLER_HPP

#include <armadillo>
class LQRController {
public:
    arma::mat computeOptimalGain(const arma::mat& A,
                                 const arma::mat& B,
                                 const arma::mat& Q,
                                 const arma::mat& R);
                                 
    arma::mat controller(const arma::vec &x, const arma::vec &uhover, const arma::vec &x0, const arma::mat &K);

private:


    arma::mat solveDARE(const arma::mat& A,
                        const arma::mat& B,
                        const arma::mat& Q,
                        const arma::mat& R);
    
    
};

#endif