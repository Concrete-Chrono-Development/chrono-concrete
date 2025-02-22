﻿// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Alessandro Tasora
// =============================================================================

#include <numeric>
#include <iomanip>

#include "chrono_modal/ChEigenvalueSolver.h"
#include "chrono_modal/ChKrylovSchurEig.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono/solver/ChDirectSolverLScomplex.h"
#include "chrono/utils/ChConstants.h"
#include "chrono/core/ChMatrix.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>

#include <Spectra/KrylovSchurGEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseRegularInverse.h>
#include <Spectra/GenEigsBase.h>

#include <unsupported/Eigen/SparseExtra>  //TODO: remove after debug

using namespace Spectra;
using namespace Eigen;

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using SpMatrix = Eigen::SparseMatrix<double>;

namespace chrono {

/// TODO: this functions returns a sparse matrix in column major order as needed by Eigen;
Eigen::Map<SpMatrix> getColMajorSparseMatrix(const ChSparseMatrix& mat) {
    return Eigen::Map<SpMatrix>(
        const_cast<ChSparseMatrix&>(mat).rows(), const_cast<ChSparseMatrix&>(mat).cols(),
        const_cast<ChSparseMatrix&>(mat).nonZeros(), const_cast<ChSparseMatrix&>(mat).outerIndexPtr(),
        const_cast<ChSparseMatrix&>(mat).innerIndexPtr(), const_cast<ChSparseMatrix&>(mat).valuePtr());
}

namespace modal {

// This is an helper class for using Krylov-Schur eigen solver also with the shift&invert mode,
// because at the moment it is not yet available in Spectra.
template <typename OpType, typename BOpType>
class KrylovSchurGEigsShiftInvert : public KrylovSchurGEigsBase<SymGEigsShiftInvertOp<OpType, BOpType>, BOpType> {
  private:
    using Scalar = typename OpType::Scalar;
    using Index = Eigen::Index;
    using Array = Eigen::Array<Scalar, Eigen::Dynamic, 1>;

    using ModeMatOp = SymGEigsShiftInvertOp<OpType, BOpType>;
    using Base = KrylovSchurGEigsBase<ModeMatOp, BOpType>;
    using Base::m_nev;
    using Base::m_ritz_val;

    const Scalar m_sigma;

    // Set shift and forward
    static ModeMatOp set_shift_and_move(ModeMatOp&& op, const Scalar& sigma) {
        op.set_shift(sigma);
        return std::move(op);
    }

    // First transform back the Ritz values, and then sort
    void sort_ritzpair(SortRule sort_rule) override {
        // The eigenvalues we get from the iteration is nu = 1 / (lambda - sigma)
        // So the eigenvalues of the original problem is lambda = 1 / nu + sigma
        m_ritz_val.head(m_nev).array() = Scalar(1) / m_ritz_val.head(m_nev).array() + m_sigma;
        Base::sort_ritzpair(sort_rule);
    }

  public:
    KrylovSchurGEigsShiftInvert(OpType& op, BOpType& Bop, Index nev, Index ncv, const Scalar& sigma)
        : Base(set_shift_and_move(ModeMatOp(op, Bop), sigma), Bop, nev, ncv), m_sigma(sigma) {}
};

void placeMatrix(Eigen::SparseMatrix<double, Eigen::ColMajor, int>& HCQ,
                 const ChSparseMatrix& H,
                 int row_start,
                 int col_start) {
    for (int k = 0; k < H.outerSize(); ++k)
        for (ChSparseMatrix::InnerIterator it(H, k); it; ++it) {
            HCQ.coeffRef(it.row() + row_start, it.col() + col_start) = it.value();
        }
}

bool ChGeneralizedEigenvalueSolverKrylovSchur::Solve(
    const ChSparseMatrix& M,   ///< input M matrix, n_v x n_v
    const ChSparseMatrix& K,   ///< input K matrix, n_v x n_v
    const ChSparseMatrix& Cq,  ///< input Cq matrix of constraint jacobians, n_c x n_v
    ChMatrixDynamic<std::complex<double>>&
        eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
    ChVectorDynamic<std::complex<double>>& eigvals,  ///< output vector with n eigenvalues, will be resized.
    ChVectorDynamic<double>& freq,       ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ChEigenvalueSolverSettings settings  ///< optional: settings for the solver, or n. of desired lower eigenvalues. If
                                         ///< =0, return all eigenvalues.
) const {
    m_timer_matrix_assembly.start();
    // Assembly the A and B for the generalized constrained eigenvalue problem.
    // Note that those sparse matrices must be column-major for better compatibility with Spectra.
    int n_vars = M.rows();
    int n_constr = Cq.rows();

    // Scale constraints matrix
    double scaling = 1.0;
    if (settings.scaleCq) {
        // std::cout << "Scaling Cq\n";
        scaling = K.diagonal().mean();
        // for (int k = 0; k < Cq.outerSize(); ++k)
        //     for (ChSparseMatrix::InnerIterator it(Cq, k); it; ++it) {
        //         it.valueRef() *= scaling;
        //     }
    }

    // A  =  [ -K   -Cq' ]
    //       [ -Cq    0  ]
    Eigen::SparseMatrix<double> A(n_vars + n_constr, n_vars + n_constr);
    A.setZero();
    placeMatrix(A, -K, 0, 0);
    placeMatrix(A, -Cq.transpose() * scaling, 0, n_vars);
    placeMatrix(A, -Cq * scaling, n_vars, 0);
    A.makeCompressed();

    // B  =  [  M     0  ]
    //       [  0     0  ]
    Eigen::SparseMatrix<double> B(n_vars + n_constr, n_vars + n_constr);
    B.setZero();
    placeMatrix(B, M, 0, 0);
    B.makeCompressed();

    m_timer_matrix_assembly.stop();
    m_timer_eigen_setup.start();

    int m = 2 * settings.n_modes >= 30 ? 2 * settings.n_modes
                                       : 30;  // minimum subspace size   //**TO DO*** make parametric?
    if (m > n_vars + n_constr - 1)
        m = n_vars + n_constr - 1;
    if (m <= settings.n_modes)
        m = settings.n_modes + 1;

    // Construct matrix operation objects using the wrapper classes
    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;
    OpType op(A, B);
    BOpType Bop(B);

    // Eigen::saveMarket(A, "C:/workspace/_temp/ChronoDump/generalized_splitmatrix_A.dat");
    // Eigen::saveMarket(B, "C:/workspace/_temp/ChronoDump/generalized_splitmatrix_B.dat");

    // The Krylov-Schur solver, using the shift and invert mode:
    KrylovSchurGEigsShiftInvert<OpType, BOpType> eigen_solver(
        op, Bop, settings.n_modes, m,
        settings.sigma
            .real());  //// TODO: OK EIGVECTS, WRONG EIGVALS REQUIRE eigen_values(i) = (1.0 / eigen_values(i)) + sigma;

    eigen_solver.init();

    m_timer_eigen_setup.stop();

    m_timer_eigen_solver.start();
    int nconv = eigen_solver.compute(SortRule::LargestMagn, settings.max_iterations, settings.tolerance);
    m_timer_eigen_solver.stop();

    if (settings.verbose) {
        if (eigen_solver.info() != CompInfo::Successful) {
            std::cout << "KrylovSchurGEigsSolver FAILED." << std::endl;
            if (eigen_solver.info() == CompInfo::NotComputed)
                std::cout << " Error: not computed." << std::endl;
            if (eigen_solver.info() == CompInfo::NotConverging)
                std::cout << " Error: not converging." << std::endl;
            if (eigen_solver.info() == CompInfo::NumericalIssue)
                std::cout << " Error: numerical issue." << std::endl;
            std::cout << " nconv  = " << nconv << std::endl;
            std::cout << " niter  = " << eigen_solver.num_iterations() << std::endl;
            std::cout << " nops   = " << eigen_solver.num_operations() << std::endl;
            return false;
        } else {
            std::cout << "KrylovSchurGEigsSolver successful." << std::endl;
            std::cout << " nconv   = " << nconv << std::endl;
            std::cout << " niter   = " << eigen_solver.num_iterations() << std::endl;
            std::cout << " nops    = " << eigen_solver.num_operations() << std::endl;
            std::cout << " n_modes = " << settings.n_modes << std::endl;
            std::cout << " n_vars  = " << n_vars << std::endl;
            std::cout << " n_constr= " << n_constr << std::endl;
        }
    }

    m_timer_solution_postprocessing.start();

    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues();
    Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();

    // ***HACK***
    // Correct eigenvals for shift-invert because KrylovSchurGEigsShiftInvert does not take care of it.
    // This should be automatically done by KrylovSchurGEigsShiftInvert::sort_ritz_pairs() at the end of compute(),
    // but at the moment such sort_ritz_pairs() is not called by the base KrylovSchurGEigsBase, differently from
    // SymGEigsShiftSolver, for example.
    for (int i = 0; i < eigen_values.rows(); ++i) {
        eigen_values(i) = (1.0 / eigen_values(i)) + settings.sigma;
    }

    // Return values
    eigvects.setZero(M.rows(), settings.n_modes);
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);

    for (int i = 0; i < settings.n_modes; i++) {
        eigvects.col(i) =
            eigen_vectors.col(i).head(n_vars);  // store only displacement part of eigenvector, no constraint part

        // normalize w.r.t. mass matrix
        double gen_mass = eigvects.col(i).real().transpose() * M * eigvects.col(i).real();
        if (gen_mass > 0)
            eigvects.col(i) *= pow(1.0 / gen_mass, 0.5);
        else
            eigvects.col(i).normalize();

        eigvals(i) = eigen_values(i);
        freq(i) = (1.0 / CH_2PI) * sqrt(-eigvals(i).real());
    }

    m_timer_solution_postprocessing.stop();

    return true;
}

bool ChGeneralizedEigenvalueSolverKrylovSchur::Solve(ChAssembly& assembly,
                                                     ChMatrixDynamic<std::complex<double>>& eigvects,
                                                     ChVectorDynamic<std::complex<double>>& eigvals,
                                                     ChVectorDynamic<double>& freq,
                                                     ChEigenvalueSolverSettings settings) const {
    // Assembly the A and B for the generalized constrained eigenvalue problem.
    // Note that those sparse matrices must be column-major for better compatibility with Spectra.

    m_timer_matrix_assembly.start();

    ChSystemDescriptor sysd;
    ChSystemDescriptor temp_descriptor;

    assembly.InjectVariables(temp_descriptor);
    assembly.InjectKRMMatrices(temp_descriptor);
    assembly.InjectConstraints(temp_descriptor);

    temp_descriptor.UpdateCountsAndOffsets();

    // Generate the A and B in state space
    int n_vars = temp_descriptor.CountActiveVariables();
    int n_constr = temp_descriptor.CountActiveConstraints();

    eigvects.setZero(n_vars, settings.n_modes);  // TODO: check row size if correct
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);

    // Eigen::SparseMatrix<double> A(n_vars + n_constr, n_vars + n_constr);
    // Eigen::SparseMatrix<double> B(n_vars + n_constr, n_vars + n_constr);

    // A  =  [ -K   -Cq' ]
    //       [ -Cq    0  ]

    // B  =  [  M     0  ]
    //       [  0     0  ]

    ChSparseMatrix A(n_vars + n_constr, n_vars + n_constr);
    ChSparseMatrix B(n_vars + n_constr, n_vars + n_constr);

    A.setZeroValues();
    B.setZeroValues();

    // Stiffness matrix
    assembly.LoadKRMMatrices(+1.0, 0.0, 0.0);
    temp_descriptor.SetMassFactor(0.0);
    temp_descriptor.PasteMassKRMMatrixInto(A, 0, 0);

    // Mass matrix
    assembly.LoadKRMMatrices(0.0, 0.0, 1.0);
    temp_descriptor.SetMassFactor(1.0);
    temp_descriptor.PasteMassKRMMatrixInto(B, 0, 0);

    // Constraint Jacobian
    assembly.LoadConstraintJacobians();
    temp_descriptor.PasteConstraintsJacobianMatrixInto(A, n_vars, 0);
    temp_descriptor.PasteConstraintsJacobianMatrixTransposedInto(A, 0, n_vars);

    // A scaling
    for (unsigned int k = 0; k < n_vars + n_constr; ++k) {
        for (ChSparseMatrix::InnerIterator it(A, k); it; ++it) {
            it.valueRef() *= -1.0;  // TODO: apply scaling here
        }
    }

    A.makeCompressed();
    B.makeCompressed();

    m_timer_matrix_assembly.stop();
    m_timer_eigen_setup.start();

    int m = 2 * settings.n_modes >= 30 ? 2 * settings.n_modes
                                       : 30;  // minimum subspace size   //**TO DO*** make parametric?
    if (m > n_vars + n_constr - 1)
        m = n_vars + n_constr - 1;
    if (m <= settings.n_modes)
        m = settings.n_modes + 1;

    // Construct matrix operation objects using the wrapper classes
    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;
    OpType op(getColMajorSparseMatrix(A), getColMajorSparseMatrix(B));
    BOpType Bop(getColMajorSparseMatrix(B));

    // Dump data for test. ***TODO*** remove when well tested
    Eigen::saveMarket(A, "C:/workspace/_temp/ChronoDump/generalized_onematrix_A.dat");
    Eigen::saveMarket(B, "C:/workspace/_temp/ChronoDump/generalized_onematrix_B.dat");
    std::ofstream fileSigma("C:/workspace/_temp/ChronoDump/generalized_onematrix_sigma.dat");
    fileSigma << settings.sigma.real();

    // The Krylov-Schur solver, using the shift and invert mode:
    KrylovSchurGEigsShiftInvert<OpType, BOpType> eigen_solver(
        op, Bop, settings.n_modes, m,
        settings.sigma
            .real());  //// TODO: OK EIGVECTS, WRONG EIGVALS REQUIRE eigen_values(i) = (1.0 / eigen_values(i)) + sigma;

    eigen_solver.init();

    m_timer_eigen_setup.start();
    m_timer_eigen_solver.start();

    int nconv = eigen_solver.compute(SortRule::LargestMagn, settings.max_iterations, settings.tolerance);
    m_timer_eigen_solver.stop();

    if (settings.verbose) {
        if (eigen_solver.info() != CompInfo::Successful) {
            std::cout << "KrylovSchurGEigsSolver FAILED." << std::endl;
            if (eigen_solver.info() == CompInfo::NotComputed)
                std::cout << " Error: not computed." << std::endl;
            if (eigen_solver.info() == CompInfo::NotConverging)
                std::cout << " Error: not converging." << std::endl;
            if (eigen_solver.info() == CompInfo::NumericalIssue)
                std::cout << " Error: numerical issue." << std::endl;
            std::cout << " nconv  = " << nconv << std::endl;
            std::cout << " niter  = " << eigen_solver.num_iterations() << std::endl;
            std::cout << " nops   = " << eigen_solver.num_operations() << std::endl;
            return false;
        } else {
            std::cout << "KrylovSchurGEigsSolver successful." << std::endl;
            std::cout << " nconv   = " << nconv << std::endl;
            std::cout << " niter   = " << eigen_solver.num_iterations() << std::endl;
            std::cout << " nops    = " << eigen_solver.num_operations() << std::endl;
            std::cout << " n_modes = " << settings.n_modes << std::endl;
            std::cout << " n_vars  = " << n_vars << std::endl;
            std::cout << " n_constr= " << n_constr << std::endl;
        }
    }

    m_timer_solution_postprocessing.start();

    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues();
    Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();

    // ***HACK***
    // Correct eigenvals for shift-invert because KrylovSchurGEigsShiftInvert does not take care of it.
    // This should be automatically done by KrylovSchurGEigsShiftInvert::sort_ritz_pairs() at the end of compute(),
    // but at the moment such sort_ritz_pairs() is not called by the base KrylovSchurGEigsBase, differently from
    // SymGEigsShiftSolver, for example.
    for (int i = 0; i < eigen_values.rows(); ++i) {
        eigen_values(i) = (1.0 / eigen_values(i)) + settings.sigma;
    }

    // Return values
    eigvects.setZero(n_vars, settings.n_modes);
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);

    for (int i = 0; i < settings.n_modes; i++) {
        eigvects.col(i) =
            eigen_vectors.col(i).head(n_vars);  // store only displacement part of eigenvector, no constraint part

        // TODO: re-enable this
        //// normalize w.r.t. mass matrix
        // double gen_mass = eigvects.col(i).real().transpose() * M * eigvects.col(i).real();
        // if (gen_mass > 0)
        //     eigvects.col(i) *= pow(1.0 / gen_mass, 0.5);
        // else
        //     eigvects.col(i).normalize();

        eigvals(i) = eigen_values(i);
        freq(i) = (1.0 / CH_2PI) * sqrt(-eigvals(i).real());
    }

    m_timer_solution_postprocessing.stop();

    return true;
}

bool ChGeneralizedEigenvalueSolverLanczos::Solve(
    const ChSparseMatrix& M,   ///< input M matrix, n_v x n_v
    const ChSparseMatrix& K,   ///< input K matrix, n_v x n_v
    const ChSparseMatrix& Cq,  ///< input Cq matrix of constraint jacobians, n_c x n_v
    ChMatrixDynamic<std::complex<double>>&
        eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
    ChVectorDynamic<std::complex<double>>& eigvals,  ///< output vector with n eigenvalues, will be resized.
    ChVectorDynamic<double>& freq,       ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ChEigenvalueSolverSettings settings  ///< optional: settings for the solver, or n. of desired lower eigenvalues. If
                                         ///< =0, return all eigenvalues.)
) const {
    m_timer_matrix_assembly.start();

    // Assembly the A and B for the generalized constrained eigenvalue problem.
    // Note that those sparse matrices must be column-major for better compatibility with Spectra.
    int n_vars = M.rows();
    int n_constr = Cq.rows();

    // Scale constraints matrix
    double scaling = 1.0;
    if (settings.scaleCq) {
        // std::cout << "Scaling Cq\n";
        scaling = K.diagonal().mean();
        // for (int k = 0; k < Cq.outerSize(); ++k)
        //     for (ChSparseMatrix::InnerIterator it(Cq, k); it; ++it) {
        //         it.valueRef() *= scaling;
        //     }
    }

    // A  =  [ -K   -Cq' ]
    //       [ -Cq    0  ]
    Eigen::SparseMatrix<double> A(n_vars + n_constr, n_vars + n_constr);
    A.setZero();
    placeMatrix(A, -K, 0, 0);
    placeMatrix(A, -Cq.transpose() * scaling, 0, n_vars);
    placeMatrix(A, -Cq * scaling, n_vars, 0);
    A.makeCompressed();

    // B  =  [  M     0  ]
    //       [  0     0  ]
    Eigen::SparseMatrix<double> B(n_vars + n_constr, n_vars + n_constr);
    B.setZero();
    placeMatrix(B, M, 0, 0);
    B.makeCompressed();

    m_timer_matrix_assembly.stop();
    m_timer_eigen_setup.start();

    int m = 2 * settings.n_modes >= 20 ? 2 * settings.n_modes : 20;  // minimum subspace size
    if (m > n_vars + n_constr - 1)
        m = n_vars + n_constr - 1;
    if (m <= settings.n_modes)
        m = settings.n_modes + 1;

    // Construct matrix operation objects using the wrapper classes
    using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = SparseSymMatProd<double>;
    OpType op(A, B);
    BOpType Bop(B);

    // The Lanczos solver, using the shift and invert mode
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eigen_solver(op, Bop, settings.n_modes, m,
                                                                              settings.sigma.real());

    eigen_solver.init();
    m_timer_eigen_setup.stop();

    m_timer_eigen_solver.start();
    int nconv = eigen_solver.compute(SortRule::LargestMagn, settings.max_iterations, settings.tolerance);
    m_timer_eigen_solver.start();

    if (settings.verbose) {
        if (eigen_solver.info() != CompInfo::Successful) {
            std::cout << "Lanczos eigenvalue solver FAILED." << std::endl;
            if (eigen_solver.info() == CompInfo::NotComputed)
                std::cout << " Error: not computed." << std::endl;
            if (eigen_solver.info() == CompInfo::NotConverging)
                std::cout << " Error: not converging." << std::endl;
            if (eigen_solver.info() == CompInfo::NumericalIssue)
                std::cout << " Error: numerical issue." << std::endl;
            std::cout << " nconv  = " << nconv << std::endl;
            std::cout << " niter  = " << eigen_solver.num_iterations() << std::endl;
            std::cout << " nops   = " << eigen_solver.num_operations() << std::endl;
            return false;
        } else {
            std::cout << "Lanczos eigenvalue solver successfull." << std::endl;
            std::cout << " nconv   = " << nconv << std::endl;
            std::cout << " niter   = " << eigen_solver.num_iterations() << std::endl;
            std::cout << " nops    = " << eigen_solver.num_operations() << std::endl;
            std::cout << " n_modes = " << settings.n_modes << std::endl;
            std::cout << " n_vars  = " << n_vars << std::endl;
            std::cout << " n_constr= " << n_constr << std::endl;
        }
    }

    m_timer_solution_postprocessing.start();

    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues();
    Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();

    // Return values
    eigvects.setZero(M.rows(), settings.n_modes);
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);

    for (int i = 0; i < settings.n_modes; i++) {
        eigvects.col(i) =
            eigen_vectors.col(i).head(n_vars);  // store only displacement part of eigenvector, no constraint part

        // normalize w.r.t. mass matrix
        double gen_mass = eigvects.col(i).real().transpose() * M * eigvects.col(i).real();
        if (gen_mass > 0)
            eigvects.col(i) *= pow(1.0 / gen_mass, 0.5);
        else
            eigvects.col(i).normalize();

        eigvals(i) = eigen_values(i);
        freq(i) = (1.0 / CH_2PI) * sqrt(-eigvals(i).real());
    }

    m_timer_solution_postprocessing.stop();

    return true;
}

int ChModalSolveUndamped::Solve(
    const ChSparseMatrix& M,   ///< input M matrix, n_v x n_v
    const ChSparseMatrix& K,   ///< input K matrix, n_v x n_v
    const ChSparseMatrix& Cq,  ///< input Cq matrix of constraint jacobians, n_c x n_v
    ChMatrixDynamic<std::complex<double>>&
        eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
    ChVectorDynamic<std::complex<double>>& eigvals,  ///< output vector with n eigenvalues, will be resized.
    ChVectorDynamic<double>& freq  ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
) const {
    int found_eigs = 0;
    eigvects.resize(0, 0);
    eigvals.resize(0);
    freq.resize(0);

    // for each freq_spans finds the closest modes to i-th input frequency:
    for (int i = 0; i < this->freq_spans.size(); ++i) {
        int nmodes_goal_i = this->freq_spans[i].nmodes;
        double sigma_i =
            -pow(this->freq_spans[i].freq * CH_2PI, 2);  // sigma for shift&invert, as lowest eigenvalue, from Hz info

        ChMatrixDynamic<std::complex<double>> eigvects_i;
        ChVectorDynamic<std::complex<double>> eigvals_i;
        ChVectorDynamic<double> freq_i;

        eigvects_i.setZero(M.rows(), nmodes_goal_i);
        eigvals_i.setZero(nmodes_goal_i);
        freq_i.setZero(nmodes_goal_i);

        ChEigenvalueSolverSettings settings_i(nmodes_goal_i, this->max_iterations, this->tolerance, this->verbose,
                                              sigma_i);

        if (!this->msolver.Solve(M, K, Cq, eigvects_i, eigvals_i, freq_i, settings_i))
            return found_eigs;

        // append to list of results

        int nmodes_out_i = eigvals_i.size();

        // Sort modes by frequencies if not exactly in increasing order. Some solver sometime fail at this.
        std::vector<int> order(nmodes_out_i);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) { return freq_i[a] < freq_i[b]; });
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
        perm.indices() = Eigen::Map<Eigen::ArrayXi>(order.data(), order.size());
        eigvects_i = eigvects_i * perm;
        eigvals_i = perm * eigvals_i;
        freq_i = perm * freq_i;

        // avoid overlap when multiple shifts were used, and too close.. If is it may happen that the lowest eigvals of
        // some shift are smaller than the highest of the previous shift.
        int i_nodes_notoverlap = nmodes_out_i;
        if (freq.size() > 0) {
            double upper_freq = freq[freq.size() - 1];
            for (int j = 0; j < nmodes_out_i; ++j)
                if (freq_i[j] < upper_freq)
                    i_nodes_notoverlap--;
        }

        if (i_nodes_notoverlap) {
            eigvects.conservativeResize(M.rows(), eigvects.cols() + i_nodes_notoverlap);
            eigvals.conservativeResize(eigvals.size() + i_nodes_notoverlap);
            freq.conservativeResize(freq.size() + i_nodes_notoverlap);
            // TODO: seems a bug in below three lines
            // because we might have i_nodes_notoverlap < nmodes_out_i when there are overlap modes
            eigvects.rightCols(i_nodes_notoverlap) = eigvects_i;
            eigvals.tail(i_nodes_notoverlap) = eigvals_i;
            freq.tail(i_nodes_notoverlap) = freq_i;
            found_eigs = eigvals.size();
        }
    }

    return found_eigs;
}

int ChModalSolveUndamped::Solve(ChAssembly& assembly,
                                ChMatrixDynamic<std::complex<double>>& eigvects,
                                ChVectorDynamic<std::complex<double>>& eigvals,
                                ChVectorDynamic<double>& freq) const {
    int found_eigs = 0;
    eigvects.resize(0, 0);
    eigvals.resize(0);
    freq.resize(0);

    // for each freq_spans finds the closest modes to i-th input frequency:
    for (int i = 0; i < this->freq_spans.size(); ++i) {
        int nmodes_goal_i = this->freq_spans[i].nmodes;
        double sigma_i =
            -pow(this->freq_spans[i].freq * CH_2PI, 2);  // sigma for shift&invert, as lowest eigenvalue, from Hz info

        ChMatrixDynamic<std::complex<double>> eigvects_i;
        ChVectorDynamic<std::complex<double>> eigvals_i;
        ChVectorDynamic<double> freq_i;

        ChEigenvalueSolverSettings settings_i(nmodes_goal_i, this->max_iterations, this->tolerance, this->verbose,
                                              sigma_i);

        if (!this->msolver.Solve(assembly, eigvects_i, eigvals_i, freq_i, settings_i))
            return found_eigs;

        // append to list of results

        int nmodes_out_i = eigvals_i.size();

        // Sort modes by frequencies if not exactly in increasing order. Some solver sometime fail at this.
        std::vector<int> order(nmodes_out_i);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) { return freq_i[a] < freq_i[b]; });
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
        perm.indices() = Eigen::Map<Eigen::ArrayXi>(order.data(), order.size());
        eigvects_i = eigvects_i * perm;
        eigvals_i = perm * eigvals_i;
        freq_i = perm * freq_i;

        // avoid overlap when multiple shifts were used, and too close.. If is it may happen that the lowest eigvals of
        // some shift are smaller than the highest of the previous shift.
        int i_nodes_notoverlap = nmodes_out_i;
        if (freq.size() > 0) {
            double upper_freq = freq[freq.size() - 1];
            for (int j = 0; j < nmodes_out_i; ++j)
                if (freq_i[j] < upper_freq)
                    i_nodes_notoverlap--;
        }

        if (i_nodes_notoverlap) {
            eigvects.conservativeResize(assembly.GetNumCoordsVelLevel(), eigvects.cols() + i_nodes_notoverlap);
            eigvals.conservativeResize(eigvals.size() + i_nodes_notoverlap);
            freq.conservativeResize(freq.size() + i_nodes_notoverlap);
            // TODO: seems a bug in below three lines
            // because we might have i_nodes_notoverlap < nmodes_out_i when there are overlap modes
            eigvects.rightCols(i_nodes_notoverlap) = eigvects_i;
            eigvals.tail(i_nodes_notoverlap) = eigvals_i;
            freq.tail(i_nodes_notoverlap) = freq_i;
            found_eigs = eigvals.size();
        }
    }
    return found_eigs;
}

//-------------------------------------------------------------------------------------------------------------------

bool ChQuadraticEigenvalueSolverNullspaceDirect::Solve(
    const ChSparseMatrix& M,
    const ChSparseMatrix& R,
    const ChSparseMatrix& K,
    const ChSparseMatrix& Cq,
    ChMatrixDynamic<std::complex<double>>& eigvects,  ///< output matrix with eigenvectors as columns, will be resized
    ChVectorDynamic<std::complex<double>>&
        eigvals,  ///< output vector with eigenvalues (real part not zero if some damping), will be resized
    ChVectorDynamic<double>& freq,           ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ChVectorDynamic<double>& damping_ratio,  ///< output vector with n damping rations r=damping/critical_damping.
    ChEigenvalueSolverSettings settings) const {
    // The folowing is adapted from the former implementation of Peng Chao

    // Compute the null space of the Cq matrix
    ChMatrixDynamic<> Cq_full = Cq;  // because the fullPivLu is not available for sparse matrix in Eigen
    Eigen::MatrixXd Cq_null_space = Cq_full.fullPivLu().kernel();
    Eigen::MatrixXd M_hat = Cq_null_space.transpose() * M * Cq_null_space;
    Eigen::MatrixXd K_hat = Cq_null_space.transpose() * K * Cq_null_space;
    Eigen::MatrixXd R_hat = Cq_null_space.transpose() * R * Cq_null_space;

    // frequency-shift，Solve the matrix singular problem - it is not needed now, you can set it to 0
    double freq_shift = 0;  // shift value, can take any value
    Eigen::MatrixXd M_bar = M_hat;
    Eigen::MatrixXd K_bar = pow(freq_shift, 2) * M_hat + freq_shift * R_hat + K_hat;
    Eigen::MatrixXd R_bar = 2 * freq_shift * M_hat + R_hat;
    Eigen::MatrixXd M_bar_inv =
        M_bar.inverse();  // performance warning! dense inverse matrix! Only small sizes should be used.

    // Generate the A matrix of the state equation, whose eigenvalues ​​are the modal frequencies
    int dim = M_bar.rows();
    Eigen::MatrixXd A_tilde(2 * dim, 2 * dim);
    A_tilde << Eigen::MatrixXd::Zero(dim, dim), Eigen::MatrixXd::Identity(dim, dim), -M_bar_inv * K_bar,
        -M_bar_inv * R_bar;

    // Call EIGEN3, dense matrix to directly solve the eigenvalues ​​and eigenvectors
    // NOTE：EIGEN3 has a very fast calculation speed in release mode, and a very slow calculation speed in debug mode
    Eigen::EigenSolver<Eigen::MatrixXd> eigen_solver(A_tilde);

    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues() + freq_shift;
    Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();

    class eig_vect_and_val {
      public:
        Eigen::VectorXcd eigen_vect;
        std::complex<double> eigen_val;

        bool operator<(const eig_vect_and_val& str) const { return (eigen_val.imag() < str.eigen_val.imag()); }
    };

    std::vector<eig_vect_and_val> all_eigen_values_and_vectors;
    for (int i = 0; i < eigen_values.size(); i++) {
        eig_vect_and_val vec_and_val{eigen_vectors.col(i), eigen_values(i)};
        all_eigen_values_and_vectors.push_back(vec_and_val);
    }

    // sort
    std::sort(all_eigen_values_and_vectors.begin(), all_eigen_values_and_vectors.end());  // sort by imaginary part

    // organize and return results
    int middle_number = (int)(all_eigen_values_and_vectors.size() /
                              2);  // The eigenvalues ​​filtered out are conjugate complex roots, just take half
    int DOF_counts = (int)((all_eigen_values_and_vectors.at(0)).eigen_vect.rows() /
                           2);  // The number of degrees of freedom in the model, used when extracting the mode shape.

    // Return values

    int nmodes = settings.n_modes;

    // if nmodes==0, then compute all eigs:
    if (nmodes == 0)
        nmodes = M.rows() - Cq.rows();

    // cannot use more modes than n. of dofs, if so, clamp
    nmodes = std::min(nmodes, (int)(M.rows() - Cq.rows()));

    eigvects.setZero(M.rows(), nmodes);
    eigvals.setZero(nmodes);
    freq.setZero(nmodes);
    damping_ratio.setZero(nmodes);

    for (int i = 0; i < nmodes; i++) {
        int i_half =
            middle_number +
            i;  // because the n.of eigenvalues is double (conjugate pairs), so just use the 2nd half after sorting
        auto mv = all_eigen_values_and_vectors.at(i_half).eigen_vect.head(DOF_counts);
        auto mvhns = Cq_null_space * mv;

        eigvects.col(i) = mvhns;

        // normalize w.r.t. mass matrix
        double gen_mass = eigvects.col(i).real().transpose() * M * eigvects.col(i).real();
        if (gen_mass > 0)
            eigvects.col(i) *= pow(1.0 / gen_mass, 0.5);
        else
            eigvects.col(i).normalize();

        eigvals(i) = all_eigen_values_and_vectors.at(i_half).eigen_val;
        freq(i) = (1.0 / CH_2PI) * (std::abs(eigvals(i)));  // undamped freq.
        damping_ratio(i) = -eigvals(i).real() / std::abs(eigvals(i));
    }

    return true;
}

//
//-------------------------------------------------------------------------------------------------------------------

ChQuadraticEigenvalueSolverKrylovSchur::ChQuadraticEigenvalueSolverKrylovSchur(
    ChDirectSolverLScomplex* mlinear_solver) {
    linear_solver = mlinear_solver;
}

bool ChQuadraticEigenvalueSolverKrylovSchur::Solve(
    const ChSparseMatrix& M,
    const ChSparseMatrix& R,
    const ChSparseMatrix& K,
    const ChSparseMatrix& Cq,
    ChMatrixDynamic<std::complex<double>>& eigvects,  ///< output matrix with eigenvectors as columns, will be resized
    ChVectorDynamic<std::complex<double>>&
        eigvals,  ///< output vector with eigenvalues (real part not zero if some damping), will be resized
    ChVectorDynamic<double>& freq,           ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ChVectorDynamic<double>& damping_ratio,  ///< output vector with n damping rations r=damping/critical_damping.
    ChEigenvalueSolverSettings settings) const {
    m_timer_matrix_assembly.start();

    // Generate the A and B in state space
    int n_vars = M.rows();
    int n_constr = Cq.rows();

    // Scale constraints matrix
    double scaling = 0;
    if (settings.scaleCq) {
        // std::cout << "Scaling Cq" << std::endl;
        scaling = K.diagonal().mean();
        for (int k = 0; k < Cq.outerSize(); ++k)
            for (ChSparseMatrix::InnerIterator it(Cq, k); it; ++it) {
                it.valueRef() *= scaling;
            }
    }

    Eigen::SparseMatrix<double, Eigen::ColMajor> As(2 * n_vars + n_constr, 2 * n_vars + n_constr);
    Eigen::SparseMatrix<double, Eigen::ColMajor> Bs(2 * n_vars + n_constr, 2 * n_vars + n_constr);
    ChSparseMatrix identity_n_vars(n_vars, n_vars);
    identity_n_vars.setIdentity();

    Eigen::VectorXi As_resSize;
    As_resSize.resize(As.cols());
    As_resSize.setZero();

    for (int k = 0; k < K.outerSize(); ++k)
        for (ChSparseMatrix::InnerIterator it(K, k); it; ++it)
            As_resSize[it.col()]++;
    for (int k = 0; k < Cq.outerSize(); ++k)
        for (ChSparseMatrix::InnerIterator it(Cq, k); it; ++it)
            As_resSize[it.col()]++;
    for (int k = 0; k < R.outerSize(); ++k)
        for (ChSparseMatrix::InnerIterator it(R, k); it; ++it)
            As_resSize[it.col() + n_vars]++;
    for (auto col_i = 0; col_i < Cq.outerSize(); col_i++) {
        As_resSize[col_i + 2 * n_vars] +=
            Cq.isCompressed() ? Cq.outerIndexPtr()[col_i + 1] - Cq.outerIndexPtr()[col_i] : Cq.innerNonZeroPtr()[col_i];
    }
    As_resSize.segment(n_vars, n_vars) = As_resSize.segment(n_vars, n_vars) + 1;
    As.reserve(As_resSize);

    // A  =  [  0     I     0 ]
    //       [ -K    -R  -Cq' ]
    //       [ -Cq    0     0 ]

    As.setZero();
    placeMatrix(As, -K, n_vars, 0);
    placeMatrix(As, -Cq, 2 * n_vars, 0);
    placeMatrix(As, identity_n_vars, 0, n_vars);
    placeMatrix(As, -R, n_vars, n_vars);
    placeMatrix(As, -Cq.transpose(), n_vars, 2 * n_vars);
    As.makeCompressed();

    Eigen::VectorXi Bs_resSize;
    Bs_resSize.resize(Bs.cols());
    Bs_resSize.setZero();
    for (int k = 0; k < M.outerSize(); ++k)
        for (ChSparseMatrix::InnerIterator it(M, k); it; ++it)
            Bs_resSize[it.col() + n_vars]++;
    Bs_resSize.segment(0, n_vars) = Bs_resSize.segment(0, n_vars) + 1;
    Bs.reserve(Bs_resSize);

    // B  =  [  I     0     0 ]
    //       [  0     M     0 ]
    //       [  0     0     0 ]
    Bs.setZero();
    placeMatrix(Bs, identity_n_vars, 0, 0);
    placeMatrix(Bs, M, n_vars, n_vars);
    Bs.makeCompressed();

    m_timer_matrix_assembly.stop();

    int n_computed_eigs = 2 * settings.n_modes;
    int m =
        2 * n_computed_eigs >= 30 ? 2 * n_computed_eigs : 30;  // minimum subspace size   //**TO DO*** make parametric
    if (m > 2 * n_vars + n_constr)
        m = 2 * n_vars + n_constr;

    // Setup the Krylov Schur solver:
    m_timer_eigen_setup.start();
    ChVectorDynamic<std::complex<double>> eigen_values;
    ChMatrixDynamic<std::complex<double>> eigen_vectors;
    ChVectorDynamic<std::complex<double>> v1;
    v1.setRandom(As.cols());  // note: to make deterministic may be preceded by something like  std::srand((unsigned
                              // int)1234567);

    // Setup the callback for matrix * vector
    callback_Ax_sparse_complexshiftinvert Ax_function3(As, Bs, settings.sigma, this->linear_solver);

    m_timer_eigen_setup.stop();

    m_timer_eigen_solver.start();

    bool isC, flag;
    int nconv, niter;
    ChKrylovSchurEig eigen_solver(
        eigen_vectors,    ///< output matrix with eigenvectors as columns, will be resized
        eigen_values,     ///< output vector with eigenvalues (real part not zero if some damping), will be resized
        isC,              ///< 0 = k-th eigenvalue is real, 1= k-th and k-th+1 are complex conjugate pairs
        flag,             ///< 0 = has converged, 1 = hasn't converged
        nconv,            ///< number of converged eigenvalues
        niter,            ///< number of used iterations
        &Ax_function3,    ///< callback for A*v
        v1,               ///< initial approx of eigenvector, or random
        As.cols(),        ///< size of A
        n_computed_eigs,  ///< number of needed eigenvalues
        m,                ///< Krylov restart threshold (largest dimension of krylov subspace)
        settings.max_iterations,  ///< max iteration number
        settings.tolerance        ///< tolerance
    );

    // Eigen::saveMarket(M, "D:/workspace/KrylovSchur-master/ChronoDump/M.dat");
    // Eigen::saveMarket(R, "D:/workspace/KrylovSchur-master/ChronoDump/R.dat");
    // Eigen::saveMarket(K, "D:/workspace/KrylovSchur-master/ChronoDump/K.dat");
    // Eigen::saveMarket(Cq, "D:/workspace/KrylovSchur-master/ChronoDump/Cq.dat");
    m_timer_eigen_solver.stop();

    // Eigen::saveMarket(As, "D:/workspace/KrylovSchur-master/ChronoDump/As.dat");
    // Eigen::saveMarket(Bs, "D:/workspace/KrylovSchur-master/ChronoDump/Bs.dat");

    // Eigen::saveMarket(ChSparseMatrix(eigen_values.real().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_values_real.dat");
    // Eigen::saveMarket(ChSparseMatrix(eigen_values.imag().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_values_imag.dat");
    // Eigen::saveMarket(ChSparseMatrix(eigen_vectors.real().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_vectors_real.dat");
    // Eigen::saveMarket(ChSparseMatrix(eigen_vectors.imag().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_vectors_imag.dat");

    if (settings.verbose) {
        if (flag == 1) {
            std::cout << "KrylovSchurEig FAILED." << std::endl;
            std::cout << " shift   = (" << settings.sigma.real() << "," << settings.sigma.imag() << ")" << std::endl;
            std::cout << " nconv = " << nconv << std::endl;
            std::cout << " niter = " << niter << std::endl;
            return false;
        } else {
            std::cout << "KrylovSchurEig successfull." << std::endl;
            std::cout << " shift   = (" << settings.sigma.real() << "," << settings.sigma.imag() << ")" << std::endl;
            std::cout << " nconv   = " << nconv << std::endl;
            std::cout << " niter   = " << niter << std::endl;
        }
    }

    m_timer_solution_postprocessing.start();

    // Restore eigenvals, trasform back  from  shift-inverted problem to original problem:
    for (int i = 0; i < eigen_values.rows(); ++i) {
        eigen_values(i) = (1.0 / eigen_values(i)) + settings.sigma;
    }

    class eig_vect_and_val {
      public:
        Eigen::VectorXcd eigen_vect;
        std::complex<double> eigen_val;

        bool operator<(const eig_vect_and_val& str) const { return (eigen_val.imag() < str.eigen_val.imag()); }
    };

    std::vector<eig_vect_and_val> all_eigen_values_and_vectors;
    for (int i = 0; i < eigen_values.size(); i++) {
        eig_vect_and_val vec_and_val{eigen_vectors.col(i), eigen_values(i)};
        all_eigen_values_and_vectors.push_back(vec_and_val);
    }

    // sort
    std::sort(all_eigen_values_and_vectors.begin(), all_eigen_values_and_vectors.end());  // sort by imaginary part

    if (settings.verbose) {
        ChVectorDynamic<std::complex<double>> resCallback(2 * n_vars + n_constr);
        for (int i = 0; i < all_eigen_values_and_vectors.size(); i++) {
            ChVectorDynamic<std::complex<double>> temp;
            Ax_function3.compute(temp, all_eigen_values_and_vectors[i].eigen_vect);
            resCallback = temp - all_eigen_values_and_vectors[i].eigen_vect * 1.0 /
                                     (all_eigen_values_and_vectors[i].eigen_val - settings.sigma);
            std::cout << "   Sorted Eig " << i << "= " << all_eigen_values_and_vectors[i].eigen_val.real() << ", "
                      << all_eigen_values_and_vectors[i].eigen_val.imag() << "i "
                      << "   freq= " << (1.0 / CH_2PI) * std::abs(all_eigen_values_and_vectors[i].eigen_val)
                      << "; res: " << resCallback.norm() << std::endl;
        }
    }

    // organize and return results
    int middle_number = (int)(all_eigen_values_and_vectors.size() /
                              2);  // The eigenvalues ​​filtered out are conjugate complex roots, just take half

    // Return values
    eigvects.setZero(M.rows(), settings.n_modes);
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);
    damping_ratio.setZero(settings.n_modes);

    for (int i = 0; i < settings.n_modes; i++) {
        // because the n.of eigenvalues is double (conjugate pairs), so just use the 2nd half after sorting
        int i_half = middle_number + i;

        // store only displacement part of eigenvector, no speed part, no constraint part
        eigvects.col(i) = all_eigen_values_and_vectors.at(i_half).eigen_vect.head(n_vars);

        // normalize w.r.t. mass matrix
        double gen_mass = eigvects.col(i).real().transpose() * M * eigvects.col(i).real();
        if (gen_mass > 0)
            eigvects.col(i) *= pow(1.0 / gen_mass, 0.5);
        else
            eigvects.col(i).normalize();

        eigvals(i) = all_eigen_values_and_vectors.at(i_half).eigen_val;
        freq(i) = (1.0 / CH_2PI) * (std::abs(eigvals(i)));  // undamped freq.
        damping_ratio(i) = -eigvals(i).real() / std::abs(eigvals(i));
    }

    return true;
}

bool ChQuadraticEigenvalueSolverKrylovSchur::Solve(
    ChAssembly& assembly,  //// TODO: cannot make const since Inject___ aren't
    ChMatrixDynamic<std::complex<double>>& eigvects,
    ChVectorDynamic<std::complex<double>>& eigvals,
    ChVectorDynamic<double>& freq,
    ChVectorDynamic<double>& damping_ratio,
    ChEigenvalueSolverSettings settings) const {
    m_timer_matrix_assembly.start();

    ChSystemDescriptor sysd;
    ChSystemDescriptor temp_descriptor;

    assembly.InjectVariables(temp_descriptor);
    assembly.InjectKRMMatrices(temp_descriptor);
    assembly.InjectConstraints(temp_descriptor);

    temp_descriptor.UpdateCountsAndOffsets();

    // Generate the A and B in state space
    int n_vars = temp_descriptor.CountActiveVariables();
    int n_constr = temp_descriptor.CountActiveConstraints();

    eigvects.setZero(n_vars, settings.n_modes);  // TODO: check row size if correct
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);
    damping_ratio.setZero(settings.n_modes);

    // Eigen::SparseMatrix<double, Eigen::ColMajor> As(2 * n_vars + n_constr, 2 * n_vars + n_constr);
    // Eigen::SparseMatrix<double, Eigen::ColMajor> Bs(2 * n_vars + n_constr, 2 * n_vars + n_constr);

    ChSparseMatrix As(2 * n_vars + n_constr, 2 * n_vars + n_constr);
    ChSparseMatrix Bs(2 * n_vars + n_constr, 2 * n_vars + n_constr);

    // A  =  [  0     I     0 ]
    //       [ -K    -R  -Cq' ]
    //       [ -Cq    0     0 ]

    // B  =  [  I     0     0 ]
    //       [  0     M     0 ]
    //       [  0     0     0 ]

    As.setZeroValues();
    Bs.setZeroValues();

    // Stiffness matrix
    assembly.LoadKRMMatrices(-1.0, 0.0, 0.0);
    temp_descriptor.SetMassFactor(0.0);
    temp_descriptor.PasteMassKRMMatrixInto(As, n_vars, 0);

    // Damping matrix
    assembly.LoadKRMMatrices(0.0, -1.0, 0.0);
    temp_descriptor.SetMassFactor(0.0);
    temp_descriptor.PasteMassKRMMatrixInto(As, n_vars, n_vars);

    // Mass matrix
    assembly.LoadKRMMatrices(0.0, 0.0, 1.0);
    temp_descriptor.SetMassFactor(1.0);
    temp_descriptor.PasteMassKRMMatrixInto(Bs, n_vars, n_vars);

    // Constraint Jacobian
    assembly.LoadConstraintJacobians();
    temp_descriptor.PasteConstraintsJacobianMatrixInto(As, 2 * n_vars, 0);
    temp_descriptor.PasteConstraintsJacobianMatrixTransposedInto(As, n_vars, 2 * n_vars);

    // Identity matrix
    for (unsigned int id_sel = 0; id_sel < n_vars; ++id_sel) {
        As.SetElement(id_sel, id_sel + n_vars, 1.0);
        Bs.SetElement(id_sel + n_vars, id_sel, 1.0);
    }

    // Cq scaling
    for (unsigned int k = 0; k < n_vars; ++k) {
        for (ChSparseMatrix::InnerIterator it(As, 2 * n_vars + k); it; ++it) {
            it.valueRef() *= -1.0;  // TODO: apply scaling here
        }
    }

    // CqT scaling
    for (unsigned int k = 0; k < n_vars; ++k) {
        for (ChSparseMatrix::InnerIterator it(As, n_vars + k); it; ++it) {
            if (it.index() >= n_vars && it.index() < 2 * n_vars) {
                it.valueRef() *= -1.0;  // TODO: apply scaling here
            }
        }
    }

    As.makeCompressed();
    Bs.makeCompressed();

    m_timer_matrix_assembly.stop();

    int n_computed_eigs = 2 * settings.n_modes;
    int m =
        2 * n_computed_eigs >= 30 ? 2 * n_computed_eigs : 30;  // minimum subspace size   //**TO DO*** make parametric
    if (m > 2 * n_vars + n_constr)
        m = 2 * n_vars + n_constr;

    // Setup the Krylov Schur solver:
    m_timer_eigen_setup.start();
    ChVectorDynamic<std::complex<double>> eigen_values;
    ChMatrixDynamic<std::complex<double>> eigen_vectors;
    ChVectorDynamic<std::complex<double>> v1;
    v1.setRandom(As.cols());  // note: to make deterministic may be preceded by something like  std::srand((unsigned
                              // int)1234567);

    // Setup the callback for matrix * vector
    callback_Ax_sparse_complexshiftinvert Ax_function3(getColMajorSparseMatrix(As), getColMajorSparseMatrix(Bs),
                                                       settings.sigma, this->linear_solver);

    m_timer_eigen_setup.stop();

    m_timer_eigen_solver.start();

    bool isC, flag;
    int nconv, niter;
    ChKrylovSchurEig eigen_solver(
        eigen_vectors,    ///< output matrix with eigenvectors as columns, will be resized
        eigen_values,     ///< output vector with eigenvalues (real part not zero if some damping), will be resized
        isC,              ///< 0 = k-th eigenvalue is real, 1= k-th and k-th+1 are complex conjugate pairs
        flag,             ///< 0 = has converged, 1 = hasn't converged
        nconv,            ///< number of converged eigenvalues
        niter,            ///< number of used iterations
        &Ax_function3,    ///< callback for A*v
        v1,               ///< initial approx of eigenvector, or random
        As.cols(),        ///< size of A
        n_computed_eigs,  ///< number of needed eigenvalues
        m,                ///< Krylov restart threshold (largest dimension of krylov subspace)
        settings.max_iterations,  ///< max iteration number
        settings.tolerance        ///< tolerance
    );

    // Eigen::saveMarket(M, "D:/workspace/KrylovSchur-master/ChronoDump/M.dat");
    // Eigen::saveMarket(R, "D:/workspace/KrylovSchur-master/ChronoDump/R.dat");
    // Eigen::saveMarket(K, "D:/workspace/KrylovSchur-master/ChronoDump/K.dat");
    // Eigen::saveMarket(Cq, "D:/workspace/KrylovSchur-master/ChronoDump/Cq.dat");
    m_timer_eigen_solver.stop();

    // Eigen::saveMarket(As, "D:/workspace/KrylovSchur-master/ChronoDump/As.dat");
    // Eigen::saveMarket(Bs, "D:/workspace/KrylovSchur-master/ChronoDump/Bs.dat");

    // Eigen::saveMarket(ChSparseMatrix(eigen_values.real().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_values_real.dat");
    // Eigen::saveMarket(ChSparseMatrix(eigen_values.imag().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_values_imag.dat");
    // Eigen::saveMarket(ChSparseMatrix(eigen_vectors.real().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_vectors_real.dat");
    // Eigen::saveMarket(ChSparseMatrix(eigen_vectors.imag().sparseView()),
    // "D:/workspace/KrylovSchur-master/ChronoDump/eigen_vectors_imag.dat");

    if (settings.verbose) {
        if (flag == 1) {
            std::cout << "KrylovSchurEig FAILED." << std::endl;
            std::cout << " shift   = (" << settings.sigma.real() << "," << settings.sigma.imag() << ")" << std::endl;
            std::cout << " nconv = " << nconv << std::endl;
            std::cout << " niter = " << niter << std::endl;
            return false;
        } else {
            std::cout << "KrylovSchurEig successfull." << std::endl;
            std::cout << " shift   = (" << settings.sigma.real() << "," << settings.sigma.imag() << ")" << std::endl;
            std::cout << " nconv   = " << nconv << std::endl;
            std::cout << " niter   = " << niter << std::endl;
        }
    }

    m_timer_solution_postprocessing.start();

    // Restore eigenvals, trasform back  from  shift-inverted problem to original problem:
    for (int i = 0; i < eigen_values.rows(); ++i) {
        eigen_values(i) = (1.0 / eigen_values(i)) + settings.sigma;
    }

    class eig_vect_and_val {
      public:
        Eigen::VectorXcd eigen_vect;
        std::complex<double> eigen_val;

        bool operator<(const eig_vect_and_val& str) const { return (eigen_val.imag() < str.eigen_val.imag()); }
    };

    std::vector<eig_vect_and_val> all_eigen_values_and_vectors;
    for (int i = 0; i < eigen_values.size(); i++) {
        eig_vect_and_val vec_and_val{eigen_vectors.col(i), eigen_values(i)};
        all_eigen_values_and_vectors.push_back(vec_and_val);
    }

    // sort
    std::sort(all_eigen_values_and_vectors.begin(), all_eigen_values_and_vectors.end());  // sort by imaginary part

    if (settings.verbose) {
        ChVectorDynamic<std::complex<double>> resCallback(2 * n_vars + n_constr);
        for (int i = 0; i < all_eigen_values_and_vectors.size(); i++) {
            ChVectorDynamic<std::complex<double>> temp;
            Ax_function3.compute(temp, all_eigen_values_and_vectors[i].eigen_vect);
            resCallback = temp - all_eigen_values_and_vectors[i].eigen_vect * 1.0 /
                                     (all_eigen_values_and_vectors[i].eigen_val - settings.sigma);
            std::cout << "   Sorted Eig " << i << "= " << all_eigen_values_and_vectors[i].eigen_val.real() << ", "
                      << all_eigen_values_and_vectors[i].eigen_val.imag() << "i "
                      << "   freq= " << (1.0 / CH_2PI) * std::abs(all_eigen_values_and_vectors[i].eigen_val)
                      << "; res: " << resCallback.norm() << std::endl;
        }
    }

    // organize and return results
    int middle_number = (int)(all_eigen_values_and_vectors.size() /
                              2);  // The eigenvalues ​​filtered out are conjugate complex roots, just take half

    // Return values
    eigvects.setZero(n_vars, settings.n_modes);
    eigvals.setZero(settings.n_modes);
    freq.setZero(settings.n_modes);
    damping_ratio.setZero(settings.n_modes);

    for (int i = 0; i < settings.n_modes; i++) {
        int i_half =
            middle_number +
            i;  // because the n.of eigenvalues is double (conjugate pairs), so just use the 2nd half after sorting

        eigvects.col(i) = all_eigen_values_and_vectors.at(i_half).eigen_vect.head(
            n_vars);  // store only displacement part of eigenvector, no speed part, no constraint part

        //// TODO: re-eneable this part
        // normalize w.r.t. mass matrix
        // double gen_mass = eigvects.col(i).real().transpose() * M * eigvects.col(i).real();
        // if (gen_mass > 0)
        //    eigvects.col(i) *= pow(1.0 / gen_mass, 0.5);
        // else
        //    eigvects.col(i).normalize();

        eigvals(i) = all_eigen_values_and_vectors.at(i_half).eigen_val;
        freq(i) = (1.0 / CH_2PI) * (std::abs(eigvals(i)));  // undamped freq.
        damping_ratio(i) = -eigvals(i).real() / std::abs(eigvals(i));
    }

    return true;
}

//-------------------------------------------------------------------------------------

int ChModalSolveDamped::Solve(
    const ChSparseMatrix& M,   ///< input M matrix, n_v x n_v
    const ChSparseMatrix& R,   ///< input R matrix, n_v x n_v
    const ChSparseMatrix& K,   ///< input K matrix, n_v x n_v
    const ChSparseMatrix& Cq,  ///< input Cq matrix of constraint jacobians, n_c x n_v
    ChMatrixDynamic<std::complex<double>>&
        eigvects,  ///< output matrix n x n_v with eigenvectors as columns, will be resized
    ChVectorDynamic<std::complex<double>>& eigvals,  ///< output vector with n eigenvalues, will be resized.
    ChVectorDynamic<double>& freq,        ///< output vector with n frequencies [Hz], as f=w/(2*PI), will be resized.
    ChVectorDynamic<double>& damp_ratios  ///< output vector with n damping ratios, will be resized.
) const {
    int found_eigs = 0;
    eigvects.resize(0, 0);
    eigvals.resize(0);
    freq.resize(0);
    damp_ratios.resize(0);

    // for each freq_spans finds the closest modes to i-th input frequency:
    for (int i = 0; i < this->freq_spans.size(); ++i) {
        int nmodes_goal_i = this->freq_spans[i].nmodes;
        double sigma_i =
            (this->freq_spans[i].freq * CH_2PI);  // sigma for shift&invert, as lowest eigenvalue, from Hz info
        // Note, for the damped case, the sigma_i assumed as a *complex* shift, as w are on the imaginary axis.

        ChMatrixDynamic<std::complex<double>> eigvects_i;
        ChVectorDynamic<std::complex<double>> eigvals_i;
        ChVectorDynamic<double> freq_i;
        ChVectorDynamic<double> damp_ratios_i;

        eigvects_i.setZero(M.rows(), nmodes_goal_i);
        eigvals_i.setZero(nmodes_goal_i);
        freq_i.setZero(nmodes_goal_i);
        damp_ratios_i.setZero(nmodes_goal_i);

        ChEigenvalueSolverSettings settings_i(nmodes_goal_i, this->max_iterations, this->tolerance, this->verbose,
                                              sigma_i);

        if (!this->msolver.Solve(M, R, K, Cq, eigvects_i, eigvals_i, freq_i, damp_ratios_i, settings_i))
            return found_eigs;

        // append to list of results

        int nmodes_out_i = eigvals_i.size();

        // Sort modes by eigenvalue imag part (or modulus(?)) if not exactly in increasing order. Some solver sometime
        // fail at this.
        std::vector<int> order(nmodes_out_i);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            // return std::abs(eigvals_i[a]) < std::abs(eigvals_i[b]);
            return eigvals_i[a].imag() < eigvals_i[b].imag();
        });
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
        perm.indices() = Eigen::Map<Eigen::ArrayXi>(order.data(), order.size());
        eigvects_i = eigvects_i * perm;
        eigvals_i = perm * eigvals_i;
        freq_i = perm * freq_i;
        damp_ratios_i = perm * damp_ratios_i;

        // avoid overlap when multiple shifts were used, and too close.. If is it may happen that the lowest eigvals of
        // some shift are smaller than the highest of the previous shift.
        int i_nodes_notoverlap = nmodes_out_i;
        if (freq.size() > 0) {
            double upper_freq = freq[freq.size() - 1];
            for (int j = 0; j < nmodes_out_i; ++j)
                if (freq_i[j] < upper_freq)
                    i_nodes_notoverlap--;
        }

        if (i_nodes_notoverlap) {
            eigvects.conservativeResize(M.rows(), eigvects.cols() + i_nodes_notoverlap);
            eigvals.conservativeResize(eigvals.size() + i_nodes_notoverlap);
            freq.conservativeResize(freq.size() + i_nodes_notoverlap);
            damp_ratios.conservativeResize(damp_ratios.size() + i_nodes_notoverlap);
            // TODO: seems a bug in below three lines
            // because we might have i_nodes_notoverlap < nmodes_out_i when there are overlap modes
            eigvects.rightCols(i_nodes_notoverlap) = eigvects_i;
            eigvals.tail(i_nodes_notoverlap) = eigvals_i;
            freq.tail(i_nodes_notoverlap) = freq_i;
            damp_ratios.tail(i_nodes_notoverlap) = damp_ratios_i;
            found_eigs = eigvals.size();
        }
    }
    return found_eigs;
}

int ChModalSolveDamped::Solve(ChAssembly& assembly,
                              ChMatrixDynamic<std::complex<double>>& eigvects,
                              ChVectorDynamic<std::complex<double>>& eigvals,
                              ChVectorDynamic<double>& freq,
                              ChVectorDynamic<double>& damp_ratios) const {
    int found_eigs = 0;
    eigvects.resize(0, 0);
    eigvals.resize(0);
    freq.resize(0);
    damp_ratios.resize(0);

    // for each freq_spans finds the closest modes to i-th input frequency:
    for (int i = 0; i < this->freq_spans.size(); ++i) {
        int nmodes_goal_i = this->freq_spans[i].nmodes;
        double sigma_i =
            (this->freq_spans[i].freq * CH_2PI);  // sigma for shift&invert, as lowest eigenvalue, from Hz info
        // Note, for the damped case, the sigma_i assumed as a *complex* shift, as w are on the imaginary axis.

        ChMatrixDynamic<std::complex<double>> eigvects_i;
        ChVectorDynamic<std::complex<double>> eigvals_i;
        ChVectorDynamic<double> freq_i;
        ChVectorDynamic<double> damp_ratios_i;

        ChEigenvalueSolverSettings settings_i(nmodes_goal_i, this->max_iterations, this->tolerance, this->verbose,
                                              sigma_i);

        if (!this->msolver.Solve(assembly, eigvects_i, eigvals_i, freq_i, damp_ratios_i, settings_i))
            return found_eigs;

        // append to list of results

        int nmodes_out_i = eigvals_i.size();

        // Sort modes by eigenvalue imag part (or modulus(?)) if not exactly in increasing order. Some solver sometime
        // fail at this.
        std::vector<int> order(nmodes_out_i);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            // return std::abs(eigvals_i[a]) < std::abs(eigvals_i[b]);
            return eigvals_i[a].imag() < eigvals_i[b].imag();
        });
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
        perm.indices() = Eigen::Map<Eigen::ArrayXi>(order.data(), order.size());
        eigvects_i = eigvects_i * perm;
        eigvals_i = perm * eigvals_i;
        freq_i = perm * freq_i;
        damp_ratios_i = perm * damp_ratios_i;

        // avoid overlap when multiple shifts were used, and too close.. If is it may happen that the lowest eigvals of
        // some shift are smaller than the highest of the previous shift.
        int i_nodes_notoverlap = nmodes_out_i;
        if (freq.size() > 0) {
            double upper_freq = freq[freq.size() - 1];
            for (int j = 0; j < nmodes_out_i; ++j)
                if (freq_i[j] < upper_freq)
                    i_nodes_notoverlap--;
        }

        if (i_nodes_notoverlap) {
            eigvects.conservativeResize(assembly.GetNumCoordsVelLevel(),
                                        eigvects.cols() + i_nodes_notoverlap);  // TODO: check row size
            eigvals.conservativeResize(eigvals.size() + i_nodes_notoverlap);
            freq.conservativeResize(freq.size() + i_nodes_notoverlap);
            damp_ratios.conservativeResize(damp_ratios.size() + i_nodes_notoverlap);
            // TODO: seems a bug in below three lines
            // because we might have i_nodes_notoverlap < nmodes_out_i when there are overlap modes
            eigvects.rightCols(i_nodes_notoverlap) = eigvects_i;
            eigvals.tail(i_nodes_notoverlap) = eigvals_i;
            freq.tail(i_nodes_notoverlap) = freq_i;
            damp_ratios.tail(i_nodes_notoverlap) = damp_ratios_i;
            found_eigs = eigvals.size();
        }
    }
    return found_eigs;
}

}  // end namespace modal

}  // end namespace chrono
