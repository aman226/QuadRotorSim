#ifndef QUADROTOR_DYNAMICS_HPP
#define QUADROTOR_DYNAMICS_HPP
#include <armadillo>
#include <memory>
#include <matplot/matplot.h>
class QuadrotorDynamics {
private:
    // Quadrotor parameters
    double m,l,g,kt,km;
    arma::mat J;
    std::shared_ptr<matplot::figure_type> fig = matplot::figure("3DAnim");
    std::vector<double> x_,y_,z_;

public:
    // Constructor
    QuadrotorDynamics(double m, double l, arma::mat J, double g, double kt, double km);
    ~QuadrotorDynamics();

    // Controller function
    arma::mat controller(const arma::vec &x, const arma::vec &uhover, const arma::vec &x0, const arma::mat &K);

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

    // quad dynamics function
    arma::vec quad_dynamics(const arma::vec &x, const arma::vec &u);

    // RK4 Integration Method
    arma::vec quad_dynamics_rk4(const arma::vec &x, const arma::vec &u, double h);
    arma::vec rk4_ss_integrator(const arma::vec &x, const arma::vec &x_int, const arma::vec &x0, const arma::vec &x_goal, const arma::mat &C, const arma::vec &r, double h);

    // Calculate Jacobian
    void jacobian(arma::mat* A, arma::mat* B, arma::colvec x, arma::colvec u, double h, double jac_eps);

    void plot_state_vs_time(const arma::mat &xhist, const arma::vec &thist);
    void plot_control_vs_time(const arma::mat &uhist, const arma::vec &thist);
    void draw(const arma::vec &x);
};

#endif