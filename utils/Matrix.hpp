#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

namespace nm
{

    class Matrix
    {
    public:
        Matrix() : rows_(0), cols_(0) {}
        Matrix(int rows, int cols, double fill = 0.0);
        Matrix(int n);
        static Matrix identity(int n);

        double &at(int i, int j);
        const double &at(int i, int j) const;
        double &operator()(int i, int j) { return at(i, j); }
        const double &operator()(int i, int j) const { return at(i, j); }

        int rows() const { return rows_; }
        int cols() const { return cols_; }
        bool isSquare() const { return rows_ == cols_; }

        Matrix operator+(const Matrix &rhs) const;
        Matrix operator-(const Matrix &rhs) const;
        Matrix operator*(const Matrix &rhs) const;
        Matrix operator*(double scalar) const;

        std::vector<double> multiply(const std::vector<double> &x) const;

        Matrix transpose() const;

        bool luDecompose(std::vector<int> &piv);

        std::vector<double> luSolve(const std::vector<double> &b) const;

        static std::vector<double> solve(Matrix A, const std::vector<double> &b);

        void print(const std::string &label = "", int width = 10, int prec = 4) const;

        double frobeniusNorm() const;

        static double vecNorm(const std::vector<double> &v);

        static double dot(const std::vector<double> &a, const std::vector<double> &b);

        static std::vector<double> normalize(const std::vector<double> &v);

    private:
        int rows_, cols_;
        std::vector<double> data_;
        std::vector<int> pivots_;

        int idx(int i, int j) const { return i * cols_ + j; }
    };

}