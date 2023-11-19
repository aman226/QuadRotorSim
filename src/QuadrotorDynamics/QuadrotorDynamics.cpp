#include "QuadrotorDynamics/QuadrotorDynamics.hpp"
#include "matplot/matplot.h"
QuadrotorDynamics::QuadrotorDynamics(double m, double l, arma::mat J, double g, double kt, double km) {
    // Quadrotor parameters
    this->m = m;
    this->l = l;
    this->J = J;
    this->g = g;
    this->kt = kt;
    this->km = km;
}

QuadrotorDynamics::~QuadrotorDynamics() {
}

// Controller function
arma::mat QuadrotorDynamics::controller(const arma::vec &x, const arma::vec &uhover, const arma::vec &x0, const arma::mat &K ) {
    arma::vec q0 = x0.subvec(3, 6);
    arma::vec q = x.subvec(3, 6);
    arma::vec phi = qtorp(L(q0).t() * q);
    arma::vec dx = arma::join_vert(arma::join_vert(arma::join_vert(x.subvec(0, 2) - x0.subvec(0, 2), phi), x.subvec(7, 9) - x0.subvec(7,9)), x.subvec(10, 12) - x0.subvec(10, 12));
    arma::vec u = uhover - K * dx;
    return u;
}

// Hat function for skew-symmetric matrix
arma::mat QuadrotorDynamics::hat(const arma::vec &v) {
    arma::mat Omega = arma::zeros<arma::mat>(3,3);
    Omega(0,1) = -v(2);
    Omega(0,2) = v(1);
    Omega(1,0) = v(2);
    Omega(1,2) = -v(0);
    Omega(2,0) = -v(1);
    Omega(2,1) = v(0);
    return Omega;
}

// L function for quaternion
arma::mat QuadrotorDynamics::L(const arma::vec &q) {
    double s = q(0);
    arma::vec v = q.subvec(1, 3);
    arma::mat L = arma::join_vert(arma::join_horiz(arma::vec({s}), -v.t()),
                                    arma::join_horiz(v, s*arma::eye<arma::mat>(3,3) + hat(v)));
    return L;
}

// Quaternion to rotation matrix
arma::mat QuadrotorDynamics::qtoQ(const arma::vec &q) {
    arma::mat T = arma::diagmat(arma::vec({1, -1, -1, -1}));
    arma::mat H = arma::join_vert(arma::zeros<arma::rowvec>(3), arma::eye<arma::mat>(3,3));
    return H.t() * T * L(q) * T * L(q) * H;
}

// G function for quaternion
arma::mat QuadrotorDynamics::G(const arma::vec &q) {
    return L(q) * arma::join_vert(arma::zeros<arma::rowvec>(3), arma::eye<arma::mat>(3,3));
}

// rotation vector to quaternion
arma::vec QuadrotorDynamics::rptoq(const arma::vec &phi) {
    return (1 / sqrt(1 + arma::dot(phi, phi))) * arma::join_vert(arma::vec({1}), phi);
}

// quaternion to rotation vector
arma::vec QuadrotorDynamics::qtorp(const arma::vec &q) {
    return q.subvec(1, 3) / q(0);
}

// E function for quaternion
arma::mat QuadrotorDynamics::E(const arma::vec &q) {
    arma::mat Gq = G(q);
    arma::mat E = arma::zeros<arma::mat>(13, 12);
    E(arma::span(0, 2), arma::span(0, 2)) = arma::eye<arma::mat>(3, 3);
    E(arma::span(3, 6), arma::span(3, 5)) = Gq;
    E(arma::span(7, 12), arma::span(6, 11)) = arma::eye<arma::mat>(6, 6);
    return E;
}


// quad dynamics function
arma::vec QuadrotorDynamics::quad_dynamics(const arma::vec &x, const arma::vec &u) {
    // get components from state vector
    arma::vec r = x.subvec(0, 2);
    arma::vec q = x.subvec(3, 6) / arma::norm(x.subvec(3, 6));
    arma::vec v = x.subvec(7, 9);
    arma::vec omega = x.subvec(10, 12);

    // calc dynamics
    arma::mat Q = qtoQ(q);
    arma::vec r_dot = Q * v;
    arma::mat H = arma::join_vert(arma::zeros<arma::rowvec>(3), arma::eye<arma::mat>(3,3));
    arma::vec q_dot = 0.5 * L(q) * H * omega;
    arma::vec v_dot = Q.t() * arma::vec({0, 0, -g}) + (1/m) * arma::join_vert(arma::zeros<arma::mat>(2, 4), arma::repmat(arma::rowvec({kt}), 1, 4)) * u - hat(omega) * v;
    arma::vec omega_dot = arma::solve(J, -hat(omega) * J * omega + arma::mat({{0, l*kt, 0, -l*kt}, {-l*kt, 0, l*kt, 0}, {km, -km, km, -km}}) * u);

    // combine and return the result
    return arma::join_vert(arma::join_vert(arma::join_vert(r_dot, q_dot), v_dot), omega_dot);
    //return arma::vec({0.0});
}

// RK4 Integration Method
arma::vec QuadrotorDynamics::quad_dynamics_rk4(const arma::vec &x, const arma::vec &u, double h) {
    arma::vec f1 = quad_dynamics(x, u);
    arma::vec f2 = quad_dynamics(x + 0.5 * h * f1, u);
    arma::vec f3 = quad_dynamics(x + 0.5 * h * f2, u);
    arma::vec f4 = quad_dynamics(x + h * f3, u);
    arma::vec xn = x + (h / 6.0) * (f1 + 2 * f2 + 2 * f3 + f4);

    // normalize quaternion
    xn.subvec(3, 6) /= arma::norm(xn.subvec(3, 6));
    return xn;
}

void QuadrotorDynamics::jacobian(arma::mat* A, arma::mat* B, arma::colvec x, arma::colvec u, double h, double jac_eps){
arma::colvec u_eps(u.n_elem);
arma::colvec x_eps(x.n_elem);

for(int i = 0; i < u.n_elem; i++){
    u_eps = u;
    u_eps(i) += jac_eps;
    if(i == 0)
        *B = (quad_dynamics_rk4(x, u_eps, h) - quad_dynamics_rk4(x, u, h)) / jac_eps;
    else
        *B = arma::join_horiz(*B, (quad_dynamics_rk4(x, u_eps, h) - quad_dynamics_rk4(x, u, h)) / jac_eps);
} 

for(int i = 0; i < x.n_elem; i++){
    x_eps = x;
    x_eps(i) += jac_eps;
    if(i == 0)
        *A = (quad_dynamics_rk4(x_eps, u, h) - quad_dynamics_rk4(x, u, h)) / jac_eps;
    else
        *A = arma::join_horiz(*A, (quad_dynamics_rk4(x_eps, u, h) - quad_dynamics_rk4(x, u, h)) / jac_eps);
}   
}

void QuadrotorDynamics::plot_state_vs_time(const arma::mat &xhist, const arma::vec &thist) {
    // Plot the state vs time
    matplot::figure("Position");
    matplot::hold(matplot::on);
    matplot::vector_1d t = arma::conv_to<matplot::vector_1d>::from(thist);
    
    matplot::vector_1d x = arma::conv_to<matplot::vector_1d>::from(xhist.row(0).t());
    matplot::vector_1d y = arma::conv_to<matplot::vector_1d>::from(xhist.row(1).t());
    matplot::vector_1d z = arma::conv_to<matplot::vector_1d>::from(xhist.row(2).t());
    matplot::plot(t, x, "-b");
    matplot::plot(t, y, "-r");
    matplot::plot(t, z, "-g");
    matplot::xlabel("Time (s)");
    matplot::ylabel("Position (m)");
    matplot::title("Position vs Time");
    matplot::legend({"x", "y", "z"});
    matplot::grid(matplot::on);
    matplot::hold(matplot::off);

    matplot::show();

    matplot::figure("3D Position");
    matplot::hold(matplot::on);
    matplot::plot3(x, y, z);
    matplot::xlabel("x (m)");
    matplot::ylabel("y (m)");
    matplot::zlabel("z (m)");
    matplot::title("3D Position");
    matplot::grid(matplot::on);
    matplot::hold(matplot::off);

    
    matplot::show();
}

void QuadrotorDynamics::plot_control_vs_time(const arma::mat &uhist, const arma::vec &thist){
    // Plot the control vs time
    matplot::figure("Control");
    matplot::hold(matplot::on);
    matplot::vector_1d t = arma::conv_to<matplot::vector_1d>::from(thist);
    
    matplot::vector_1d u1 = arma::conv_to<matplot::vector_1d>::from(uhist.row(0).t());
    matplot::vector_1d u2 = arma::conv_to<matplot::vector_1d>::from(uhist.row(1).t());
    matplot::vector_1d u3 = arma::conv_to<matplot::vector_1d>::from(uhist.row(2).t());
    matplot::vector_1d u4 = arma::conv_to<matplot::vector_1d>::from(uhist.row(3).t());
    matplot::plot(t, u1, "-b");
    matplot::plot(t, u2, "-r");
    matplot::plot(t, u3, "-g");
    matplot::plot(t, u4, "-k");
    matplot::xlabel("Time (s)");
    matplot::ylabel("Control Input (N)");
    matplot::title("Control Input vs Time");
    matplot::legend({"u1", "u2", "u3", "u4"});
    matplot::grid(matplot::on);
    matplot::hold(matplot::off);

    matplot::show();
}