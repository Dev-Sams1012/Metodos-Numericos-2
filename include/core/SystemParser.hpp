#pragma once

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>
#include <muParser.h>

namespace nm
{

    class SystemParser
    {
    public:
        SystemParser() = default;
        ~SystemParser() = default;

        SystemParser(const SystemParser &) = delete;
        SystemParser &operator=(const SystemParser &) = delete;

        void setup(int n, const std::vector<std::string> &exprs);

        std::vector<double> eval(const std::vector<double> &s, double t) const;

        int size() const { return n_; }
        bool isReady() const { return ready_; }

        const std::string &expression(int i) const { return exprs_[i]; }

    private:
        int n_{0};
        bool ready_{false};
        std::vector<std::string> exprs_;

        mutable std::vector<double> sVars_;
        mutable double tVar_{0.0};

        std::vector<std::unique_ptr<mu::Parser>> parsers_;
    };

}