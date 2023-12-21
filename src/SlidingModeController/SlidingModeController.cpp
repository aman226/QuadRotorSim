#include <SlidingModeController/SlidingModeController.hpp>

SlidingModeController::SlidingModeController(){

}

SlidingModeController::~SlidingModeController(){
    return;
}

arma::vec SlidingModeController::get_s(const arma::vec* state, const arma::mat* J){
    
}
arma::vec SlidingModeController::generate_control_input(const arma::vec* state){
    arma::vec control_input, omega, p, p_dot;
    return control_input;
}