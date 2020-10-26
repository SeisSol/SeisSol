#include "Functions.h"

namespace seissol {
namespace functions {

uint64_t rangeProduct(uint64_t from, uint64_t to) {
    uint64_t product = 1;
    for (; from <= to; ++from) {
        product *= from;
    }
    return product;
}

double JacobiP(unsigned n, unsigned a, unsigned b, double x) {
    if (n == 0) {
        return 1.0;
    }
    double Pm_2;
    double Pm_1 = 1.0;
    double Pm = 0.5 * a - 0.5 * b + (1.0 + 0.5 * (a + b)) * x;
    double a2_b2 = static_cast<double>(a * a) - static_cast<double>(b * b);
    for (unsigned m = 2; m <= n; ++m) {
        Pm_2 = Pm_1;
        Pm_1 = Pm;
        Pm = ((2.0 * m + a + b - 1.0) * (a2_b2 + (2.0 * m + a + b) * (2.0 * m + a + b - 2.0) * x) *
                  Pm_1 -
              2.0 * (m + a - 1.0) * (m + b - 1.0) * (2.0 * m + a + b) * Pm_2) /
             (2.0 * m * (m + a + b) * (2.0 * m + a + b - 2.0));
    }
    return Pm;
}

double JacobiPDerivative(unsigned n, unsigned a, unsigned b, double x) {
    return (n == 0) ? 0.0 : 0.5 * (n + a + b + 1.0) * JacobiP(n - 1, a + 1, b + 1, x);
}

std::array<double, 5> SingularityFreeJacobiPFactors(unsigned m, unsigned a, unsigned b) {
    double c_0 = 2.0 * m + a + b;
    double c_1 = c_0 - 1.0;
    double c_2 = static_cast<double>(a * a) - static_cast<double>(b * b);
    double c_3 = c_0 * (c_0 - 2.0);
    double c_4 = 2.0 * (m + a - 1.0) * (m + b - 1.0) * c_0;
    double c_5 = 2.0 * m * (m + a + b) * (c_0 - 2.0);
    return {c_1, c_2, c_3, c_4, c_5};
}

double SingularityFreeJacobiPRecursion(double x, double y, std::array<double, 5> const& cm,
                                       double Pm_1, double Pm_2) {
    return (cm[0] * (cm[1] * y + cm[2] * x) * Pm_1 - cm[3] * y * y * Pm_2) / cm[4];
}

double SingularityFreeJacobiP(unsigned n, unsigned a, unsigned b, double x, double y) {
    if (n == 0) {
        return 1.0;
    }
    double Pm_2;
    double Pm_1 = 1.0;
    double Pm = (0.5 * a - 0.5 * b) * y + (1.0 + 0.5 * (a + b)) * x;
    for (unsigned m = 2; m <= n; ++m) {
        Pm_2 = Pm_1;
        Pm_1 = Pm;
        auto c = SingularityFreeJacobiPFactors(m, a, b);
        Pm = SingularityFreeJacobiPRecursion(x, y, c, Pm_1, Pm_2);
    }
    return Pm;
}

std::array<double, 3> SingularityFreeJacobiPAndDerivatives(unsigned n, unsigned a, unsigned b,
                                                           double x, double y) {
    if (n == 0) {
        return {1.0, 0.0, 0.0};
    }
    double Pm_2, ddxPm_2, ddyPm_2;
    double Pm_1 = 1.0, ddxPm_1 = 0.0, ddyPm_1 = 0.0;
    double Pm = SingularityFreeJacobiP(1, a, b, x, y);
    double ddxPm = 1.0 + 0.5 * (a + b);
    double ddyPm = 0.5 * (static_cast<double>(a) - static_cast<double>(b));
    for (unsigned m = 2; m <= n; ++m) {
        Pm_2 = Pm_1;
        Pm_1 = Pm;
        ddxPm_2 = ddxPm_1;
        ddxPm_1 = ddxPm;
        ddyPm_2 = ddyPm_1;
        ddyPm_1 = ddyPm;
        auto c = SingularityFreeJacobiPFactors(m, a, b);
        Pm = SingularityFreeJacobiPRecursion(x, y, c, Pm_1, Pm_2);
        ddxPm = (c[0] * (c[2] * Pm_1 + (c[1] * y + c[2] * x) * ddxPm_1) - c[3] * y * y * ddxPm_2) /
                c[4];
        ddyPm = (c[0] * (c[1] * Pm_1 + (c[1] * y + c[2] * x) * ddyPm_1) -
                 c[3] * (2.0 * y * Pm_2 + y * y * ddyPm_2)) /
                c[4];
    }
    return {Pm, ddxPm, ddyPm};
}

double TriDubinerP(std::array<unsigned, 2> const& i, std::array<double, 2> const& xi) {
    double r_num = 2.0 * xi[0] - 1.0 + xi[1];
    double s = 2.0 * xi[1] - 1.0;
    double theta = 1.0 - xi[1];

    double ti = SingularityFreeJacobiP(i[0], 0, 0, r_num, theta);
    double tij = SingularityFreeJacobiP(i[1], 2 * i[0] + 1, 0, s, 1.0);

    return ti * tij;
}

std::array<double, 2> gradTriDubinerP(std::array<unsigned, 2> const& i,
                                      std::array<double, 2> const& xi) {
    double r_num = 2.0 * xi[0] - 1.0 + xi[1];
    double s = 2.0 * xi[1] - 1.0;
    double theta = 1.0 - xi[1];

    auto ti = SingularityFreeJacobiPAndDerivatives(i[0], 0, 0, r_num, theta);
    auto tij = SingularityFreeJacobiPAndDerivatives(i[1], 2 * i[0] + 1, 0, s, 1.0);

    auto ddalpha = [&](double dr_num, double dtheta, double dt) {
        return (ti[1] * dr_num + ti[2] * dtheta) * tij[0] + ti[0] * tij[1] * dt;
    };

    return {ddalpha(2.0, 0.0, 0.0), ddalpha(1.0, -1.0, 2.0)};
}

double TetraDubinerP(std::array<unsigned, 3> const& i, std::array<double, 3> const& xi) {
    double r_num = 2.0 * xi[0] - 1.0 + xi[1] + xi[2];
    double s_num = 2.0 * xi[1] - 1.0 + xi[2];
    double t = 2.0 * xi[2] - 1.0;
    double sigmatheta = 1.0 - xi[1] - xi[2];
    double theta = 1.0 - xi[2];

    double ti = SingularityFreeJacobiP(i[0], 0, 0, r_num, sigmatheta);
    double tij = SingularityFreeJacobiP(i[1], 2 * i[0] + 1, 0, s_num, theta);
    double tijk = SingularityFreeJacobiP(i[2], 2 * i[0] + 2 * i[1] + 2, 0, t, 1.0);

    return ti * tij * tijk;
}

std::array<double, 3> gradTetraDubinerP(std::array<unsigned, 3> const& i,
                                        std::array<double, 3> const& xi) {
    double r_num = 2.0 * xi[0] - 1.0 + xi[1] + xi[2];
    double s_num = 2.0 * xi[1] - 1.0 + xi[2];
    double t = 2.0 * xi[2] - 1.0;
    double sigmatheta = 1.0 - xi[1] - xi[2];
    double theta = 1.0 - xi[2];

    auto ti = SingularityFreeJacobiPAndDerivatives(i[0], 0, 0, r_num, sigmatheta);
    auto tij = SingularityFreeJacobiPAndDerivatives(i[1], 2 * i[0] + 1, 0, s_num, theta);
    auto tijk = SingularityFreeJacobiPAndDerivatives(i[2], 2 * i[0] + 2 * i[1] + 2, 0, t, 1.0);

    auto ddalpha = [&](double dr_num, double dsigmatheta, double ds_num, double dtheta, double dt) {
        return (ti[1] * dr_num + ti[2] * dsigmatheta) * tij[0] * tijk[0] +
               ti[0] * (tij[1] * ds_num + tij[2] * dtheta) * tijk[0] +
               ti[0] * tij[0] * (tijk[1] * dt);
    };

    return {ddalpha(2.0, 0.0, 0.0, 0.0, 0.0), ddalpha(1.0, -1.0, 2.0, 0.0, 0.0),
            ddalpha(1.0, -1.0, 1.0, -1.0, 2.0)};
}

template <>
double DubinerP<1u>(std::array<unsigned, 1u> const& i, std::array<double, 1u> const& xi) {
    return JacobiP(i[0], 0, 0, 2.0 * xi[0] - 1.0);
}
template <>
double DubinerP<2u>(std::array<unsigned, 2u> const& i, std::array<double, 2u> const& xi) {
    return TriDubinerP(i, xi);
}
template <>
double DubinerP<3u>(std::array<unsigned, 3u> const& i, std::array<double, 3u> const& xi) {
    return TetraDubinerP(i, xi);
}

template <>
std::array<double, 1u> gradDubinerP<1u>(std::array<unsigned, 1u> const& i,
                                        std::array<double, 1u> const& xi) {
    return {JacobiPDerivative(i[0], 0, 0, 2.0 * xi[0] - 1.0)};
}
template <>
std::array<double, 2u> gradDubinerP<2u>(std::array<unsigned, 2u> const& i,
                                        std::array<double, 2u> const& xi) {
    return gradTriDubinerP(i, xi);
}
template <>
std::array<double, 3u> gradDubinerP<3u>(std::array<unsigned, 3u> const& i,
                                        std::array<double, 3u> const& xi) {
    return gradTetraDubinerP(i, xi);
}

} // namespace functions
} // namespace seissol
