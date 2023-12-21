#ifndef _SMC_HPP
#define _SMC_HPP

#include <armadillo>

class SlidingModeController{
    private:
    public:
        SlidingModeController();
        ~SlidingModeController();
        arma::vec get_s(const arma::vec*, const arma::mat *);
        arma::vec generate_control_input(const arma::vec*);

};
#endif