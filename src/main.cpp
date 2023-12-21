#include <armadillo>
#include <iostream>
#include <chrono>
#include <LQRController/LQRController.hpp>
#include <QuadrotorDynamics/QuadrotorDynamics.hpp>
#include <QuaternionUtils/quaternion_op.hpp>

int main() {
    unsigned int microsecond = 1000000;
    double Tfinal = 100.0;
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
    
    // state vector for linearizing dynamics
    arma::vec r0 = {1.0, 1.0, 1.0}; // Initial position
    arma::vec q0 = {1.0, 0.0, 0.0, 0.0}; // Initial orientation (quaternion)
    arma::vec v0 = arma::zeros<arma::vec>(3); // Initial linear velocity
    arma::vec omega0 = arma::zeros<arma::vec>(3); // Initial angular velocity
    arma::vec x0 = arma::join_cols(arma::join_cols(arma::join_cols(r0, q0), v0), omega0); // Combined initial state
    arma::vec uhover = (m * g / 4.0) * arma::ones<arma::vec>(4); // Hover control input
    arma::mat A, B;

    // define cost matrices
    arma::mat Q = arma::diagmat(arma::vec({1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));
    arma::mat R = arma::diagmat(arma::vec({1, 1, 1, 1})) * 1e-1;

    // define LQR controller
    LQRController lqr;

    arma::mat x_goal = x0;
    arma::mat x_start =  arma::join_vert(arma::join_vert(arma::join_vert(r0 + 4*arma::randn(3), rptoq(arma::vec({1.0, 0.0, 0.0}))), v0), omega0);
    // initialize control input and state vectors
    arma::mat uhist = arma::zeros<arma::mat>(4, Nt);
    arma::mat xhist = arma::zeros<arma::mat>(13, Nt);

    arma::vec thist = arma::linspace<arma::vec>(0, h * (Nt - 1), Nt); // Time vector

    auto start = std::chrono::high_resolution_clock::now();
    
    // Linearize
    quadDynamics.jacobian(&A, &B, x0, uhover, h, 1e-2);
    
    arma::mat C({{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}});

    arma::mat A_red = E(q0).t() * A * E(q0);
    arma::mat B_red = E(q0).t() * B;

    // Trajectory
    arma::vec theta = arma::linspace<arma::vec>(0, 2*arma::datum::pi, 50);
    double radius = 2.0;
    arma::vec x_circle = radius*arma::cos(theta);
    arma::vec y_circle = radius*arma::sin(theta);
    arma::vec z_circle = arma::vec(arma::ones(50));


    // compute optimal gain matrix
    arma::mat K = lqr.computeOptimalGain(A_red, B_red, Q, R);
    xhist.col(0) = x_start;
    
    quadDynamics.draw(xhist.col(0),x_circle,y_circle,z_circle);
    arma::vec ref;
    double traj_cnt = 0;
    // simulate the system
    for (int k = 0; k < Nt - 1; ++k) {
        ref = arma::vec({x_circle(traj_cnt), y_circle(traj_cnt), 1.0});
        uhist.col(k) = lqr.controller(xhist.col(k), uhover, x_goal, K) + K(arma::span(0, 3), arma::span(0, 2)) * (ref - C*x_goal);
        xhist.col(k + 1) = quadDynamics.quad_dynamics_rk4(xhist.col(k), uhist.col(k), h);
        quadDynamics.draw(xhist.col(k + 1),x_circle,y_circle,z_circle);
        double dist_to_ref = arma::norm(xhist(arma::span(0, 2), k) - ref);

        if (dist_to_ref < 0.1) {
            traj_cnt += 1;
            if (traj_cnt > 49) {
            traj_cnt = 0;
            }
        }
        
        usleep(0.05 * microsecond);
    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout << "Time taken by Sim: " << duration.count() << " microseconds" << std::endl;
    quadDynamics.plot_state_vs_time(xhist, thist);
    quadDynamics.plot_control_vs_time(uhist, thist);
    matplot::show();
    return 0;
}
