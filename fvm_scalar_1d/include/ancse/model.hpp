#ifndef FVMSCALAR1D_MODEL_HPP
#define FVMSCALAR1D_MODEL_HPP

#include <cmath>

// This is intentionally not generic, to avoid more boilerplate code.
// This example shows Burgers' equation. Change the example to something else
// you find interesting.
class Model {
  public:
    inline double flux(double u) const { return 0.5 * u * u; }

    inline double max_eigenvalue(double u) const { return std::abs(u); }

    inline double min_point() const{return 0.0;}
};

/// This is useful is you just need 'some' valid model.
Model make_dummy_model();

#endif // FVMSCALAR1D_MODEL_HPP
