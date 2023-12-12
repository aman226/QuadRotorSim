#pragma once

#include <armadillo>
#include "glob_opts.hpp"

#ifdef __cplusplus
extern "C"
{
#endif

    typedef arma::dmat::fixed<NSTATES, 1ULL> tiny_VectorNx;
    typedef arma::dmat::fixed<NINPUTS, 1> tiny_VectorNu;
    typedef arma::dmat::fixed<NSTATES, NSTATES> tiny_MatrixNxNx;
    typedef arma::dmat::fixed<NSTATES, NINPUTS> tiny_MatrixNxNu;
    typedef arma::dmat::fixed<NINPUTS, NSTATES> tiny_MatrixNuNx;
    typedef arma::dmat::fixed<NINPUTS, NINPUTS> tiny_MatrixNuNu;

    typedef arma::dmat::fixed<NSTATES, NHORIZON> tiny_MatrixNxNh;       // Nu x Nh
    typedef arma::dmat::fixed<NINPUTS, NHORIZON - 1> tiny_MatrixNuNhm1; // Nu x Nh-1

    /**
     * Matrices that must be recomputed with changes in time step, rho
     */
    typedef struct
    {
        tinytype rho;
        tiny_MatrixNuNx Kinf;
        tiny_MatrixNxNx Pinf;
        tiny_MatrixNuNu Quu_inv;
        tiny_MatrixNxNx AmBKt;
        tiny_MatrixNxNu coeff_d2p;
    } TinyCache;

    /**
     * User settings
     */
    typedef struct
    {
        tinytype abs_pri_tol;
        tinytype abs_dua_tol;
        int max_iter;
        int check_termination;
        int en_state_bound;
        int en_input_bound;
    } TinySettings;

    /**
     * Problem variables
     */
    typedef struct
    {
        // State and input
        tiny_MatrixNxNh x;
        tiny_MatrixNuNhm1 u;

        // Linear control cost terms
        tiny_MatrixNxNh q;
        tiny_MatrixNuNhm1 r;

        // Linear Riccati backward pass terms
        tiny_MatrixNxNh p;
        tiny_MatrixNuNhm1 d;

        // Auxiliary variables
        tiny_MatrixNxNh v;
        tiny_MatrixNxNh vnew;
        tiny_MatrixNuNhm1 z;
        tiny_MatrixNuNhm1 znew;

        // Dual variables
        tiny_MatrixNxNh g;
        tiny_MatrixNuNhm1 y;

        tinytype primal_residual_state;
        tinytype primal_residual_input;
        tinytype dual_residual_state;
        tinytype dual_residual_input;
        int status;
        int iter;

        tiny_VectorNx Q;
        tiny_VectorNx Qf;
        tiny_VectorNu R;
        tiny_MatrixNxNx Adyn;
        tiny_MatrixNxNu Bdyn;

        tiny_MatrixNuNhm1 u_min;
        tiny_MatrixNuNhm1 u_max;
        tiny_MatrixNxNh x_min;
        tiny_MatrixNxNh x_max;
        tiny_MatrixNxNh Xref;   // Nx x Nh
        tiny_MatrixNuNhm1 Uref; // Nu x Nh-1

        // Temporaries
        tiny_VectorNu Qu;
    } TinyWorkspace;

    /**
     * Main TinyMPC solver structure that holds all information.
     */
    typedef struct
    {
        TinySettings *settings; // Problem settings
        TinyCache *cache;       // Problem cache
        TinyWorkspace *work;    // Solver workspace
    } TinySolver;

#ifdef __cplusplus
}
#endif
