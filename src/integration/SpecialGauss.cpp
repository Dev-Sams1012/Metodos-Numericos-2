#include "integration/SpecialGauss.hpp"
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace nm
{

    SpecialGauss::SpecialGauss(std::shared_ptr<FunctionParser> f,
                               SpecialGaussType type, int n)
        : func_(std::move(f)), type_(type), n_(n)
    {
        loadTable();
    }

    void SpecialGauss::setType(SpecialGaussType t)
    {
        type_ = t;
        loadTable();
    }
    void SpecialGauss::setPoints(int n)
    {
        n_ = n;
        loadTable();
    }

    void SpecialGauss::loadTable()
    {
        switch (type_)
        {
        case SpecialGaussType::HERMITE:
            loadHermite(n_);
            break;
        case SpecialGaussType::LAGUERRE:
            loadLaguerre(n_);
            break;
        case SpecialGaussType::CHEBYSHEV:
            loadChebyshev(n_);
            break;
        }
    }

    void SpecialGauss::loadHermite(int n)
    {
        switch (n)
        {
        case 2:
            nodes_ = {-0.7071067811865476, 0.7071067811865476};
            weights_ = {0.8862269254527580, 0.8862269254527580};
            break;

        case 3:
            nodes_ = {-1.2247448713915890, 0.0, 1.2247448713915890};
            weights_ = {0.2954089751509193, 1.1816359006036774, 0.2954089751509193};
            break;

        case 4:
            nodes_ = {
                -1.6506801238857845, -0.5246476232752903,
                0.5246476232752903, 1.6506801238857845};
            weights_ = {
                0.08131283544724518, 0.8049140900055128,
                0.8049140900055128, 0.08131283544724518};
            break;

        case 5:
            nodes_ = {
                -2.0201828704560856, -0.9585724646137554,
                0.0,
                0.9585724646137554, 2.0201828704560856};

            weights_ = {
                0.019953242059045923, 0.39361932315224975,
                0.9453087204829419,
                0.39361932315224975, 0.019953242059045923};
            break;

        default:
            throw std::invalid_argument("Gauss-Hermite: apenas n=2,3,4,5 implementado.");
        }
    }

    void SpecialGauss::loadLaguerre(int n)
    {
        switch (n)
        {
        case 2:
            nodes_ = {0.5857864376269050, 3.4142135623730951};
            weights_ = {0.8535533905932738, 0.1464466094067262};
            break;

        case 3:
            nodes_ = {0.4157745567834791, 2.2942803602790417, 6.2899450829374791};
            weights_ = {0.7110930099291730, 0.2785177335692408, 0.0103892565015861};
            break;

        case 4:
            nodes_ = {
                0.3225476896193923, 1.7457611011583466,
                4.5366202969211280, 9.3950709123011331};
            weights_ = {
                0.6031541043416336, 0.3574186924377996,
                0.03888790851500538, 0.0005392947055613274};
            break;

        case 5:
            nodes_ = {
                0.26356031971814109, 1.4134030591065167,
                3.5964257710407243, 7.0858100058588376,
                12.640800844275783};

            weights_ = {
                0.52175561058280873, 0.39866681108315199,
                0.075942449681707595, 0.0036117586799220484,
                0.00002336997238577622};
            break;

        default:
            throw std::invalid_argument("Gauss-Laguerre: apenas n=2,3,4,5 implementado.");
        }
    }

    void SpecialGauss::loadChebyshev(int n)
    {
        if (n < 1)
            throw std::invalid_argument("Gauss-Chebyshev: n deve ser >= 1.");
        nodes_.resize(n);
        weights_.resize(n);
        const double pi = std::acos(-1.0);
        for (int k = 1; k <= n; ++k)
        {
            nodes_[k - 1] = std::cos((2.0 * k - 1.0) * pi / (2.0 * n));
            weights_[k - 1] = pi / n;
        }
    }

    double SpecialGauss::compute() const
    {
        if (!func_ || !func_->isReady())
            throw std::runtime_error("SpecialGauss: função não definida.");
        double result = 0.0;
        for (int k = 0; k < (int)nodes_.size(); ++k)
            result += weights_[k] * func_->eval(nodes_[k]);
        return result;
    }

    std::string SpecialGauss::name() const
    {
        std::ostringstream oss;
        switch (type_)
        {
        case SpecialGaussType::HERMITE:
            oss << "Gauss-Hermite   (n=" << n_ << ", peso e^{-x²})";
            break;
        case SpecialGaussType::LAGUERRE:
            oss << "Gauss-Laguerre  (n=" << n_ << ", peso e^{-x})";
            break;
        case SpecialGaussType::CHEBYSHEV:
            oss << "Gauss-Chebyshev (n=" << n_ << ", peso 1/√(1-x²))";
            break;
        }
        return oss.str();
    }

}