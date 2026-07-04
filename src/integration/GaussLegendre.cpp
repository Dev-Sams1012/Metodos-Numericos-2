#include "integration/GaussLegendre.hpp"
#include <stdexcept>
#include <sstream>

namespace nm
{

    GaussLegendre::GaussLegendre(std::shared_ptr<FunctionParser> f, int n)
        : Integrator(std::move(f))
    {
        setPoints(n);
    }

    void GaussLegendre::setPoints(int n)
    {
        if (n < 2 || n > 5)
            throw std::invalid_argument("GaussLegendre: n deve estar entre 2 e 5.");
        n_ = n;
        loadTable(n);
    }

    void GaussLegendre::loadTable(int n)
    {
        nodes_.clear();
        weights_.clear();

        switch (n)
        {
        case 2:

            nodes_ = {-0.5773502691896257, 0.5773502691896257};
            weights_ = {1.0, 1.0};
            break;

        case 3:

            nodes_ = {-0.7745966692414834, 0.0, 0.7745966692414834};
            weights_ = {0.5555555555555556, 0.8888888888888888,
                        0.5555555555555556};
            break;

        case 4:

            nodes_ = {
                -0.8611363115940526, -0.3399810435848563,
                0.3399810435848563, 0.8611363115940526};
            weights_ = {
                0.3478548451374538, 0.6521451548625461,
                0.6521451548625461, 0.3478548451374538};
            break;

        case 5:

            nodes_ = {
                -0.9061798459386640, -0.5384693101056831,
                0.0,
                0.5384693101056831, 0.9061798459386640};
            weights_ = {
                0.2369268850561891, 0.4786286704993665,
                0.5688888888888889,
                0.4786286704993665, 0.2369268850561891};
            break;
        }
    }

    double GaussLegendre::compute(double a, double b) const
    {

        const double mid = (a + b) / 2.0;
        const double half = (b - a) / 2.0;
        double result = 0.0;
        for (int k = 0; k < n_; ++k)
        {
            double xk = mid + half * nodes_[k];
            result += weights_[k] * eval(xk);
        }
        return half * result;
    }

    std::string GaussLegendre::name() const
    {
        std::ostringstream oss;
        oss << "Gauss-Legendre (n=" << n_ << ")";
        return oss.str();
    }

}