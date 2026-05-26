#pragma once
// =============================================================================
// Matrix.hpp
// Classe Matrix para operações de álgebra linear usadas nos métodos de
// autovalores/autovetores. Inclui decomposição LU para o método da potência
// inverso e operações básicas (produto, transposta, norma).
// =============================================================================
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

namespace nm {

class Matrix {
public:
    // ---------- Construtores -------------------------------------------------
    Matrix() : rows_(0), cols_(0) {}
    Matrix(int rows, int cols, double fill = 0.0);
    Matrix(int n);                                  ///< Matriz quadrada n×n
    static Matrix identity(int n);                  ///< Matriz identidade n×n

    // ---------- Acesso a elementos -------------------------------------------
    double&       at(int i, int j);
    const double& at(int i, int j) const;
    double&       operator()(int i, int j)       { return at(i, j); }
    const double& operator()(int i, int j) const { return at(i, j); }

    int rows() const { return rows_; }
    int cols() const { return cols_; }
    bool isSquare() const { return rows_ == cols_; }

    // ---------- Operações matriciais -----------------------------------------
    Matrix  operator+(const Matrix& rhs) const;
    Matrix  operator-(const Matrix& rhs) const;
    Matrix  operator*(const Matrix& rhs) const;
    Matrix  operator*(double scalar)     const;

    /// @brief Multiplica matriz por vetor coluna: y = A * x
    std::vector<double> multiply(const std::vector<double>& x) const;

    Matrix  transpose() const;

    // ---------- Decomposição LU com pivotamento parcial ---------------------
    /// @brief Decompõe A = P*L*U (in-place). Preenche piv com permutações.
    /// @return false se a matriz for singular.
    bool luDecompose(std::vector<int>& piv);

    /// @brief Resolve Ax = b após luDecompose ter sido chamado.
    std::vector<double> luSolve(const std::vector<double>& b) const;

    /// @brief Cria cópia e resolve via LU: A * x = b.
    static std::vector<double> solve(Matrix A, const std::vector<double>& b);

    // ---------- Utilitários --------------------------------------------------
    void print(const std::string& label = "", int width = 10, int prec = 4) const;

    /// @brief Norma de Frobenius
    double frobeniusNorm() const;

    /// @brief Norma do vetor
    static double vecNorm(const std::vector<double>& v);

    /// @brief Produto interno de dois vetores
    static double dot(const std::vector<double>& a, const std::vector<double>& b);

    /// @brief Normaliza vetor (retorna vetor unitário)
    static std::vector<double> normalize(const std::vector<double>& v);

private:
    int rows_, cols_;
    std::vector<double> data_;        ///< Armazenamento row-major
    std::vector<int>    pivots_;      ///< Permutações da decomposição LU

    int idx(int i, int j) const { return i * cols_ + j; }
};

} // namespace nm