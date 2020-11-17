#ifndef FVMSCALAR1D_NUMERICAL_FLUX_HPP
#define FVMSCALAR1D_NUMERICAL_FLUX_HPP

#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model. It is also unconditionally a
 * bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    Model model;
};

//Rusanov's

class Rusanov {
  public:

    explicit Rusanov(const Model &model) : model(model){}

    double operator()(double uL, double uR) const {

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        auto fL_d = model.max_eigenvalue(uL);
        auto fR_d = model.max_eigenvalue(uR);

        auto c = std::max(fL_d, fR_d);


        return 0.5 * (fL + fR) - c/2*(uR - uL);

    }

  private:
    Model model;
};
//

//Enquist-Osher
class Enquist_osher {
  public:

    explicit Enquist_osher(const Model &model) : model(model){}

    double operator()(double uL, double uR) const {

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);


        auto D = 0.0;

        if (uR > 0 && uL > 0){
            D = fR - fL;
        } else if (uR < 0 && uL < 0){
            D = fL - fR;
        } else if (uR > 0 && uL < 0){
            D = fL + fR;
        } else if (uR < 0 && uL > 0){
            D = fL + fR;
        }

        return 0.5 * (fL + fR) - 0.5 * D;

    }

  private:
    Model model;
};

//Godnunov

class Godunov {
  public:

    explicit Godunov(const Model &model) : model(model){}

    double operator()(double uL, double uR) const {

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        auto F = 0.0;
        if (uR >= uL && uL >0){
            F = fL;
        } else if (uR >= 0 && uL <= 0){
            F = 0.0;
        } else if (uR < 0 && uR >= uL){
            F = fR;
        } else if (uL > uR && uR >= 0){
            F = fL;
        } else if (uR < 0 && uL >= 0){
            F = std::max(fL, fR);
        } else if (uL > uR && uL < 0){
            F = fR;
        }

        return F;

    }

  private:
    Model model;
};

//Roe

class Roe{
  public:

    explicit Roe(const Model &model) : model(model){}

    double operator()(double uL, double uR) const {

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        auto A = 0.0;
        if (uR == uL){
            A = uL;
        } else {
            A = (fR-fL)/(uR-uL);
        }

        auto F = 0.0;
        if (A >= 0){
            F = fL;
        } else if (A < 0){
            F = fR;
        }

        return F;

    }

  private:
    Model model;
};



//----------------FluxLFBegin----------------
/// Lax-Friedrichs numerical flux.
/** This flux works for any model. */
class LaxFriedrichs {
  public:
    // Note: This version is a bit tricky. A numerical flux should be
    //       a function of the two trace values at the interface, i.e. what we
    //       call `uL`, `uR`. However, it requires 'dt' and 'dx'. Therefore,
    //       these need to be made available to the flux. This is one of the
    //       reasons why `SimulationTime`.
    LaxFriedrichs(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {
        double dx = grid.dx;
        double dt = simulation_time->dt;

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * (fL + fR) - 0.5 * dx/dt*(uR - uL);

    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};
//----------------FluxLFEnd----------------


#endif // FVMSCALAR1D_NUMERICAL_FLUX_HPP
