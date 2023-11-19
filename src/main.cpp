#include <armadillo>
#include <iostream>
#include <chrono>
#include <LQRController/LQRController.hpp>
#include <QuadrotorDynamics/QuadrotorDynamics.hpp>

int main() {
    double Tfinal = 10.0;
    double h = 0.05; 
    int Nt = (int) (Tfinal / h) + 1;

    // Quadrotor parameters
    double m = 0.5;
    double l = 0.1750;
    arma::mat J = arma::diagmat(arma::vec({0.0023, 0.0023, 0.004}));
    double g = 9.81;
    double kt = 1.0;
    double km = 0.0245;

    // create an instance of QuadrotorDynamics
    QuadrotorDynamics quadDynamics(m, l, J, g, kt, km);
    
    // initial state vector
    arma::vec r0 = {0.0, 0.0, 1.0}; // Initial position
    arma::vec q0 = {1.0, 0.0, 0.0, 0.0}; // Initial orientation (quaternion)
    arma::vec v0 = arma::zeros<arma::vec>(3); // Initial linear velocity
    arma::vec omega0 = arma::zeros<arma::vec>(3); // Initial angular velocity
    arma::vec uhover = (m * g / 4.0) * arma::ones<arma::vec>(4); // Hover control input
    arma::vec x0 = arma::join_cols(arma::join_cols(arma::join_cols(r0, q0), v0), omega0); // Combined initial state

    arma::mat A, B;
    quadDynamics.jacobian(&A, &B, x0, uhover, h, 1e-2);

    arma::mat A_red = quadDynamics.E(q0).t() * A * quadDynamics.E(q0);
    arma::mat B_red = quadDynamics.E(q0).t() * B;

    // define cost matrices
    arma::mat Q = arma::diagmat(arma::vec({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));
    arma::mat R = arma::diagmat(arma::vec({1, 1, 1, 1})) * 1e-1;

    // define LQR controller
    LQRController lqr(A_red, B_red, Q, R);

    // compute optimal gain matrix
    arma::mat K = lqr.computeOptimalGain();

    // initialize control input and state vectors
    arma::mat uhist = arma::zeros<arma::mat>(4, Nt);
    arma::mat xhist = arma::zeros<arma::mat>(13, Nt);
    arma::vec thist = arma::linspace<arma::vec>(0, h * (Nt - 1), Nt); // Time vector

    xhist.col(0) = arma::join_vert(arma::join_vert(arma::join_vert(r0 + 2*arma::randn(3), quadDynamics.rptoq(arma::vec({1.0, 0.0, 0.0}))), v0), omega0);
    
    // simulate the system
    for (int k = 0; k < Nt - 1; ++k) {
        uhist.col(k) = quadDynamics.controller(xhist.col(k), uhover, x0, K);
        xhist.col(k + 1) = quadDynamics.quad_dynamics_rk4(xhist.col(k), uhist.col(k), h);
    }
    
    quadDynamics.plot_state_vs_time(xhist, thist);
    quadDynamics.plot_control_vs_time(uhist, thist);
    return 0;
}
