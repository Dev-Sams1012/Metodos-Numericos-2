#include "Matrix.hpp"
#include <cassert>
#include <numeric>
#include <sstream>
#include <algorithm>

namespace nm
{

    Matrix::Matrix(int rows, int cols, double fill)
        : rows_(rows), cols_(cols), data_(rows * cols, fill) {}

    Matrix::Matrix(int n) : Matrix(n, n) {}

    Matrix Matrix::identity(int n)
    {
        Matrix m(n, n, 0.0);
        for (int i = 0; i < n; ++i)
            m.at(i, i) = 1.0;
        return m;
    }

    double &Matrix::at(int i, int j)
    {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
            throw std::out_of_range("Matrix::at — índice fora dos limites.");
        return data_[idx(i, j)];
    }

    const double &Matrix::at(int i, int j) const
    {
        if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
            throw std::out_of_range("Matrix::at — índice fora dos limites.");
        return data_[idx(i, j)];
    }

    Matrix Matrix::operator+(const Matrix &rhs) const
    {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_)
            throw std::invalid_argument("Matrix::+ — dimensões incompatíveis.");
        Matrix result(rows_, cols_);
        for (int k = 0; k < rows_ * cols_; ++k)
            result.data_[k] = data_[k] + rhs.data_[k];
        return result;
    }

    Matrix Matrix::operator-(const Matrix &rhs) const
    {
        if (rows_ != rhs.rows_ || cols_ != rhs.cols_)
            throw std::invalid_argument("Matrix::- — dimensões incompatíveis.");
        Matrix result(rows_, cols_);
        for (int k = 0; k < rows_ * cols_; ++k)
            result.data_[k] = data_[k] - rhs.data_[k];
        return result;
    }

    Matrix Matrix::operator*(const Matrix &rhs) const
    {
        if (cols_ != rhs.rows_)
            throw std::invalid_argument("Matrix::* — dimensões incompatíveis.");
        Matrix result(rows_, rhs.cols_, 0.0);
        for (int i = 0; i < rows_; ++i)
            for (int k = 0; k < cols_; ++k)
                for (int j = 0; j < rhs.cols_; ++j)
                    result.data_[i * rhs.cols_ + j] += data_[idx(i, k)] * rhs.data_[k * rhs.cols_ + j];
        return result;
    }

    Matrix Matrix::operator*(double scalar) const
    {
        Matrix result(rows_, cols_);
        for (int k = 0; k < rows_ * cols_; ++k)
            result.data_[k] = data_[k] * scalar;
        return result;
    }

    std::vector<double> Matrix::multiply(const std::vector<double> &x) const
    {
        if ((int)x.size() != cols_)
            throw std::invalid_argument("Matrix::multiply — tamanho do vetor incompatível.");
        std::vector<double> y(rows_, 0.0);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < cols_; ++j)
                y[i] += at(i, j) * x[j];
        return y;
    }

    Matrix Matrix::transpose() const
    {
        Matrix T(cols_, rows_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < cols_; ++j)
                T.at(j, i) = at(i, j);
        return T;
    }

    bool Matrix::luDecompose(std::vector<int> &piv)
    {
        if (!isSquare())
            throw std::invalid_argument("luDecompose — matriz deve ser quadrada.");
        int n = rows_;
        piv.resize(n);
        std::iota(piv.begin(), piv.end(), 0);
        pivots_ = piv;

        for (int k = 0; k < n; ++k)
        {

            int maxRow = k;
            double maxVal = std::abs(at(k, k));
            for (int i = k + 1; i < n; ++i)
            {
                if (std::abs(at(i, k)) > maxVal)
                {
                    maxVal = std::abs(at(i, k));
                    maxRow = i;
                }
            }
            if (maxVal < 1e-14)
                return false;

            if (maxRow != k)
            {
                for (int j = 0; j < n; ++j)
                    std::swap(data_[idx(k, j)], data_[idx(maxRow, j)]);
                std::swap(piv[k], piv[maxRow]);
            }

            for (int i = k + 1; i < n; ++i)
            {
                at(i, k) /= at(k, k);
                for (int j = k + 1; j < n; ++j)
                    at(i, j) -= at(i, k) * at(k, j);
            }
        }
        pivots_ = piv;
        return true;
    }

    std::vector<double> Matrix::luSolve(const std::vector<double> &b) const
    {
        int n = rows_;

        std::vector<double> x(n);
        for (int i = 0; i < n; ++i)
            x[i] = b[pivots_[i]];

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                x[i] -= at(i, j) * x[j];

        for (int i = n - 1; i >= 0; --i)
        {
            for (int j = i + 1; j < n; ++j)
                x[i] -= at(i, j) * x[j];
            x[i] /= at(i, i);
        }
        return x;
    }

    std::vector<double> Matrix::solve(Matrix A, const std::vector<double> &b)
    {
        std::vector<int> piv;
        if (!A.luDecompose(piv))
            throw std::runtime_error("Matrix::solve — matriz singular.");
        return A.luSolve(b);
    }

    void Matrix::print(const std::string &label, int width, int prec) const
    {
        if (!label.empty())
            std::cout << label << "\n";
        for (int i = 0; i < rows_; ++i)
        {
            std::cout << "  [";
            for (int j = 0; j < cols_; ++j)
                std::cout << std::setw(width) << std::fixed << std::setprecision(prec) << at(i, j);
            std::cout << " ]\n";
        }
    }

    double Matrix::frobeniusNorm() const
    {
        double s = 0.0;
        for (auto v : data_)
            s += v * v;
        return std::sqrt(s);
    }

    double Matrix::vecNorm(const std::vector<double> &v)
    {
        double s = 0.0;
        for (auto x : v)
            s += x * x;
        return std::sqrt(s);
    }

    double Matrix::dot(const std::vector<double> &a, const std::vector<double> &b)
    {
        if (a.size() != b.size())
            throw std::invalid_argument("dot — vetores de tamanhos diferentes.");
        double s = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            s += a[i] * b[i];
        return s;
    }

    std::vector<double> Matrix::normalize(const std::vector<double> &v)
    {
        double n = vecNorm(v);
        if (n < 1e-14)
            throw std::runtime_error("normalize — vetor nulo.");
        std::vector<double> u(v.size());
        for (size_t i = 0; i < v.size(); ++i)
            u[i] = v[i] / n;
        return u;
    }

}