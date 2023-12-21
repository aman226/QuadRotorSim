#include "QuadrotorDynamics/QuadrotorDynamics.hpp"
#include "QuaternionUtils/quaternion_op.hpp"
#include <iostream>
#include <omp.h>
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

arma::vec QuadrotorDynamics::rk4_ss_integrator(const arma::vec &x, const arma::vec &x_int, const arma::vec &x0, const arma::vec &x_goal, const arma::mat &C, const arma::vec &r, double h){
    return x_int + h*(-C*x + r - x_goal(arma::span(0,0)) - C*x0);
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
    std::shared_ptr<matplot::figure_type> pos_fig = matplot::figure("Position");
    matplot::hold(matplot::on);
    matplot::vector_1d t = arma::conv_to<matplot::vector_1d>::from(thist);
    
    matplot::vector_1d x = arma::conv_to<matplot::vector_1d>::from(xhist.row(0).t());
    matplot::vector_1d y = arma::conv_to<matplot::vector_1d>::from(xhist.row(1).t());
    matplot::vector_1d z = arma::conv_to<matplot::vector_1d>::from(xhist.row(2).t());
    matplot::plot(t, x, "-b");
    matplot::plot(t, y, "-r");
    matplot::plot(t, z, "-g");
    matplot::xlim({-4.0,4.0});
    matplot::ylim({-4.0,4.0});
    matplot::zlim({0.0,2.5});
    matplot::xlabel("Time (s)");
    matplot::ylabel("Position (m)");
    matplot::title("Position vs Time");
    matplot::legend({"x", "y", "z"});
    matplot::grid(matplot::on);
    matplot::hold(matplot::off);

    pos_fig->draw();

    std::shared_ptr<matplot::figure_type> pos_3d_fig = matplot::figure("3D Position");
    matplot::hold(matplot::on);
    matplot::plot3(x, y, z);
    matplot::xlabel("x (m)");
    matplot::ylabel("y (m)");
    matplot::zlabel("z (m)");
    matplot::title("3D Position");
    matplot::grid(matplot::on);
    matplot::hold(matplot::off);
    pos_3d_fig->draw();
}

void QuadrotorDynamics::plot_control_vs_time(const arma::mat &uhist, const arma::vec &thist){
    // Plot the control vs time
    std::shared_ptr<matplot::figure_type> control_fig = matplot::figure("Control");
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
    control_fig->draw();
}

void QuadrotorDynamics::draw(const arma::vec &x, const arma::vec &x_traj_, const arma::vec &y_traj_, const arma::vec &z_traj_){
    this->x_.push_back(x(0));
    this->y_.push_back(x(1));
    this->z_.push_back(x(2));
    matplot::vector_1d  x_traj,y_traj,z_traj;
    if (!x_traj_.is_empty()){
        x_traj = arma::conv_to<matplot::vector_1d>::from(x_traj_);
        y_traj = arma::conv_to<matplot::vector_1d>::from(y_traj_);
        z_traj = arma::conv_to<matplot::vector_1d>::from(z_traj_);
    }

    // Get the rotation matrix
    arma::vec q = x.subvec(3, 6);
    arma::mat Q = qtoQ(q);
    Q = Q.t();
    arma::vec r1 = Q.col(0)/arma::norm(Q.col(0),2)*0.5;
    arma::vec r2 = Q.col(1)/arma::norm(Q.col(1),2)*0.5;
    arma::vec r3 = Q.col(2)/arma::norm(Q.col(2),2)*0.5;

    std::vector<double> X = {x(0),x(0),x(0)};
    std::vector<double> Y = {x(1),x(1),x(1)};
    std::vector<double> Z = {x(2),x(2),x(2)};

    std::vector<double> U = {r1(0),r2(0),r3(0)};
    std::vector<double> V = {r1(1),r2(1),r3(1)};
    std::vector<double> W = {r1(2),r2(2),r3(2)};
    
    matplot::hold(matplot::on);
    matplot::cla();
    matplot::quiver3(X, Y, Z, U, V, W)->line_width(3);

    if (!x_traj_.is_empty())
        matplot::plot3(x_traj,y_traj,z_traj,"g--")->line_width(3);

    matplot::plot3(this->x_,this->y_,this->z_,"k-")->line_width(3);
    matplot::hold(matplot::off);
    matplot::xlabel("x (m)");
    matplot::ylabel("y (m)");
    matplot::zlabel("z (m)");
    fig->width(1920);
    fig->height(1080);
    matplot::title("3D Position");
    matplot::xlim({-4.0,4.0});
    matplot::ylim({-4.0,4.0});
    matplot::zlim({0.0,2.5});
    fig->draw();
}