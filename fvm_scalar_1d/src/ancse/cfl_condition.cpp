#include <ancse/cfl_condition.hpp>


//----------------StandardCFLConditionDefnBegin----------------
StandardCFLCondition::StandardCFLCondition(const Grid &grid,
                                           const Model &model,
                                           double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double StandardCFLCondition::operator()(const Eigen::VectorXd &u) const {

    auto dx = grid.dx;

    auto n_cells = grid.n_cells;
    auto n_ghost = grid.n_ghost;

    auto u_max = u[n_ghost];

    for (int i = n_ghost; i < n_cells-n_ghost; i++) {
        if (abs(u[i]) > u_max){
            u_max = abs(u[i]);
        }
    }


    auto dt = cfl_number*dx/u_max;

    return dt;
}
//----------------StandardCFLConditionDefnEnd----------------


std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const Model &model, double cfl_number) {
    // implement this 'factory' for your CFL condition.
    return std::make_shared<StandardCFLCondition>(grid, model, cfl_number);
    return nullptr;
}
