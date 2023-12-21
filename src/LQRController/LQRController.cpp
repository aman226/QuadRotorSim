#include <LQRController/LQRController.hpp>
#include <QuaternionUtils/quaternion_op.hpp>

arma::mat LQRController::computeOptimalGain(const arma::mat& A,
                                            const arma::mat& B,
                                            const arma::mat& Q,
                                            const arma::mat& R) {
    return solveDARE(A, B, Q, R);
}

arma::mat LQRController::solveDARE(const arma::mat& A,
                                   const arma::mat& B,
                                   const arma::mat& Q,
                                   const arma::mat& R) {
    // Initial guess for P
    arma::mat P = Q;

    // Convergence tolerance
    double tol = 1e-5;
    int maxIter = 1000;
    for (int i = 0; i < maxIter; ++i) {
        arma::mat P_next = A.t() * P * A - (A.t() * P * B) * arma::inv(R + B.t() * P * B) * (B.t() * P * A) + Q;
        
        // Check for convergence
        if (arma::norm(P_next - P, "fro") < tol) {
            P = P_next;
            break;
        }

        P = P_next;
    }

    // Calculate and return the optimal gain matrix
    return arma::inv(R + B.t() * P * B) * (B.t() * P * A);
}

// Controller function
arma::mat LQRController::controller(const arma::vec &x, const arma::vec &uhover, const arma::vec &x0, const arma::mat &K) {
    arma::vec q0 = x0.subvec(3, 6);
    arma::vec q = x.subvec(3, 6);
    arma::vec phi = qtorp(L(q0).t() * q);
    arma::vec dx = arma::join_vert(arma::join_vert(arma::join_vert(x.subvec(0, 2) - x0.subvec(0, 2), phi), x.subvec(7, 9) - x0.subvec(7,9)), x.subvec(10, 12) - x0.subvec(10, 12));
    arma::vec u = uhover - K * dx;
    return u;
}
