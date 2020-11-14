#ifndef FVMSCALAR1D_RECONSTRUCTION_HPP
#define FVMSCALAR1D_RECONSTRUCTION_HPP

#include <Eigen/Dense>
#include <cmath>
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

inline double sign(double a) { return copysign(1.0, a); }



//----------------SlopeLimiterABegin----------------
inline double minmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::min(std::abs(a), std::abs(b));
}

inline double maxmod(double a, double b){
    return 0.5 * (sign(a) + sign(b)) * std::max(std::abs(a), std::abs(b));
}

inline double minabs(double a, double b){
    if (abs(a) < abs(b)){
        return a;
    } else if (abs(a) == abs(b)){
        return (a+b)/2;
    } else if (abs(a) > abs(b)){
        return b;
    }
}

inline double tri_minmod(double a, double b, double c){
    return 1/3 * (sign(a) + sign(b) + sign(c)) * std::min(std::min(std::abs(a), std::abs(b)), std::abs(c));

}

//----------------SlopeLimiterAEnd----------------

//----------------SlopeLimiterBBegin----------------
struct MinMod {
    inline double operator()(double sL, double sR) const {
        return minmod(sL, sR);
    }
};

struct SuperBee{
    inline double operator()(double sL, double sR) const {
        auto SL = minmod(2*sL, sR);
        auto SR = minmod( sL, 2*sR);
        return maxmod(SL, SR);
    };
};

struct MinAbs{
    inline double operator()(double sL, double sR) const {
        return minabs(sL, sR);
    }
};

struct MonotonizedCentral{
    inline double operator()(double sL, double sR) const {
        return tri_minmod(2*sR, (sR + sL)/2, 2*sL);
    }
};

struct VanLeer{
    inline double operator()(double sL, double sR) const {
        return (std::abs(sR) * sL + std::abs(sL) * sR)/(std::abs(sR) + std::abs(sL));
    }

};
//----------------SlopeLimiterBEnd----------------


class PWConstantReconstruction {
  public:
    /// Compute the left and right trace at the interface i + 1/2.
    /** Note: This API is agnostic to the number of cell-averages required
     *        by the method. Therefore, reconstructions with different stencil
     *        sizes can implement this API; and this call can be used in parts
     *        of the code that do not need to know about the details of the
     *        reconstruction.
     */
    inline std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                                int i) const {
        return (*this)(u[i], u[i + 1]);
    }

    /// Compute the left and right trace at the interface.
    /** Piecewise constant reconstruction of the left and right trace only
     *  requires the cell-average to the left and right of the interface.
     *
     *  Note: Compared to the other overload this reduces the assumption on
     *        how the cell-averages are stored. This is useful when testing and
     *        generally makes the function useful in more situations.
     */
    inline std::pair<double, double> operator()(double ua, double ub) const {
        return {ua, ub};
    }
};


//----------------LinearRCBegin----------------
template <class SlopeLimiter>
class PWLinearReconstruction {
  public:
    explicit PWLinearReconstruction(const Grid &grid,const SlopeLimiter &slope_limiter);

    std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                         int i) const {
        return (*this)(u[i - 1], u[i], u[i + 1], u[i + 2]);
    }

    std::pair<double, double>
    operator()(double ua, double ub, double uc, double ud) const {
        double dx = grid.dx;

        double uL = 0.0;
        double uR = 0.0;
        auto sL1 = (ub-ua)/dx;
        auto sR1 = (uc-ub)/dx;
        auto sR2 = (ud-uc)/dx;
        uL = ub + dx/2 * slope_limiter(sL1, sR1);
        uR = uc - dx/2 * slope_limiter(sR1, sR2);

        return {uL, uR};
    }

  private:
    Grid grid;
    SlopeLimiter slope_limiter;
};
template <class SlopeLimiter>
PWLinearReconstruction<SlopeLimiter>::PWLinearReconstruction(
    const Grid &grid, const SlopeLimiter &slope_limiter)
    : grid(grid),
      slope_limiter(slope_limiter) {}
//----------------LinearRCEnd----------------


#endif // FVMSCALAR1D_RATE_OF_CHANGE_HPP
