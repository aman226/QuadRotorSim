#include <LQRController/LQRController.hpp>

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
