#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>
#include "funkcja.h"

class Jet {
    double f, dx, dy, dxx, dxy, dyy;

public:
    Jet() : f(0.0), dx(0.0), dy(0.0), dxx(0.0), dxy(0.0), dyy(0.0) {}

    Jet(const double value,
        const double dx,
        const double dy,
        const double dxx,
        const double dxy,
        const double dyy) : f(value), dx(dx), dy(dy), dxx(dxx), dxy(dxy), dyy(dyy) {};

    Jet(const Jet& jet) : f(jet.f), dx(jet.dx), dy(jet.dy), dxx(jet.dxx), dxy(jet.dxy), dyy(jet.dyy) {};

    friend std::ostream &operator<<(std::ostream &ostream, const Jet& jet) {
        ostream << (jet.f == -0 ? 0 : jet.f) << " "
                << (jet.dx == -0 ? 0 : jet.dx) << " "
                << (jet.dy == -0 ? 0 : jet.dy) << " "
                << (jet.dxx == -0 ? 0 : jet.dxx) << " "
                << (jet.dxy == -0 ? 0 : jet.dxy) << " "
                << (jet.dyy == -0 ? 0 : jet.dyy);
        return ostream;
    }

    // operator przypisania
    inline Jet &operator=(const Jet& jet) = default;

    // negacja
    inline Jet operator-() const {
        return Jet(-f, -dx, -dy, -dxx, -dxy, -dyy);
    }

    // c - jet
    inline friend Jet operator-(const double& c, const Jet& jet) {
        return Jet(c, 0., 0., 0., 0., 0.) - jet;
    }

    // jet - c
    inline friend Jet operator-(const Jet& jet, const double& c) {
        return Jet(jet.f - c, jet.dx, jet.dy, jet.dxx, jet.dxy, jet.dyy);
    }

    // jet - jet
    inline friend Jet operator-(const Jet& u, const Jet& v) {
        return Jet(u.f - v.f, u.dx - v.dx, u.dy - v.dy,
                   u.dxx - v.dxx, u.dxy - v.dxy, u.dyy - v.dyy);
    }

    // c + jet
    inline friend Jet operator+(const double& c, const Jet& jet) {
        return Jet(c, 0., 0., 0., 0., 0.) + jet;
    }

    // jet + c
    inline friend Jet operator+(const Jet& jet, const double& c) {
        return Jet(jet.f + c, jet.dx, jet.dy, jet.dxx, jet.dxy, jet.dyy);
    }

    // jet + jet
    inline friend Jet operator+(const Jet& u, const Jet& v) {
        return Jet(u.f + v.f, u.dx + v.dx, u.dy + v.dy,
                   u.dxx + v.dxx, u.dxy + v.dxy, u.dyy + v.dyy);
    }

    // jet * jet
    inline friend Jet operator*(const Jet& u, const Jet& v) {
        return Jet( u.f * v.f,
                    u.f * v.dx + u.dx * v.f,
                    u.f * v.dy + u.dy * v.f,
                    u.f * v.dxx + 2 * u.dx * v.dx + u.dxx * v.f,
                    u.f * v.dxy + u.dx * v.dy + u.dy * v.dx + u.dxy * v.f,
                    u.f * v.dyy + 2 * u.dy * v.dy + u.dyy * v.f);
    }

    // jet * c
    inline friend Jet operator*(const Jet& jet, const double& c) {
        return c * jet;
    }

    // c * jet
    inline friend Jet operator*(const double& c, const Jet& jet) {
        return jet * Jet(c, 0., 0., 0., 0., 0.);
    }

    // jet / jet
    inline friend Jet operator/(const Jet& u, const Jet& v) {
        double quotient = u.f / v.f;
        double res_dx = (u.dx - (quotient * v.dx)) / v.f;
        double res_dy = (u.dy - (quotient * v.dy)) / v.f;
        double neg_v = -v.f;

        return Jet(quotient,
                   res_dx,
                   res_dy,
                   (2 * v.dx * res_dx + quotient * v.dxx -u.dxx) / neg_v,
                   (v.dx * res_dy + v.dy * res_dx + quotient * v.dxy - u.dxy) / neg_v,
                   (2 * v.dy * res_dy + quotient * v.dyy - u.dyy) / neg_v);
    }

    // jet / c
    inline friend Jet operator/(const Jet& jet, const double& c) {
        return (jet / Jet(c, 0., 0., 0., 0., 0.));
    }

    // c / jet
    inline friend Jet operator/(const double& c, const Jet& jet) {
        return Jet(c / jet.f, 0., 0., 0., 0., 0.);
    }

    inline friend Jet sin(const Jet& jet) {
        double sine = std::sin(jet.f);
        double cosine = std::cos(jet.f);

        return Jet(sine,
                   jet.dx * cosine,
                   jet.dy * cosine,
                   jet.dxx * cosine - jet.dx * jet.dx * sine,
                   jet.dxy * cosine - jet.dx * jet.dy * sine,
                   jet.dyy * cosine - jet.dy * jet.dy * sine);
    }

    inline friend Jet cos(const Jet& jet) {
        double sine = std::sin(jet.f);
        double cosine = std::cos(jet.f);

        return Jet(cosine,
                   -jet.dx * sine,
                   -jet.dy * sine,
                   -(jet.dxx * sine + jet.dx * jet.dx * cosine),
                   -(jet.dxy * sine + jet.dx * jet.dy * cosine),
                   -(jet.dyy * sine + jet.dy * jet.dy * cosine));
    }

    inline friend Jet exp(const Jet& jet) {
        double e = std::exp(jet.f);
        return Jet(e,
                   jet.dx * e,
                   jet.dy * e,
                   jet.dxx * e + jet.dx * jet.dx * e,
                   jet.dxy * e + jet.dx * jet.dy * e,
                   jet.dyy * e + jet.dy * jet.dy * e);
    }
};

int main() {
    unsigned int M;
    std::cin >> M;
    assert (M < 1000000 and M > 0);
    double x0, y0;
    std::cout << std::fixed << std::setprecision(15);

    for (int i = 0; i < M; ++i) {
        std::cin >> x0 >> y0;

        Jet jet_x(x0, 1., 0., 0., 0., 0.);
        Jet jet_y(y0, 0., 1., 0., 0., 0.);

        Jet res = funkcja(jet_x, jet_y);
        std::cout << res;
        std::cout << std::endl;
    }
    return 0;
}
