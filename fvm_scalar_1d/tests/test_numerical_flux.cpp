#include <gtest/gtest.h>
#include <ancse/numerical_flux.hpp>

template <class NumericalFlux>
void check_consistency(const NumericalFlux &nf) {
    auto model = make_dummy_model();

    auto to_check = std::vector<double>{1.0, 2.0, 3.3251, -1.332};
    for (auto u : to_check) {
        ASSERT_DOUBLE_EQ(model.flux(u), nf(u, u));
    }
}

TEST(TestCentralFlux, consistency) {
    auto model = make_dummy_model();
    auto central_flux = CentralFlux(model);

    check_consistency(central_flux);
}

TEST(TestRusanovFlux, consistency) {
    auto model = make_dummy_model();
    auto Rusanov_flux = Rusanov(model);

    check_consistency(Rusanov_flux);
}

TEST(TestRoeFlux, consistency){
    auto model = make_dummy_model();
    auto Roe_flux = Roe(model);

    check_consistency(Roe_flux);
}

TEST(TestEnquistOsherFlux, consistency){
    auto model = make_dummy_model();
    auto Enquist_Osher_flux = Enquist_osher(model);

    check_consistency(Enquist_Osher_flux);
}

TEST(TestGodunovFlux, consistency){
    auto model = make_dummy_model();
    auto Godunov_flux = Godunov(model);

    check_consistency(Godunov_flux);
}


TEST(TestLaxFriedrichsFlux, consistency) {
    auto model = make_dummy_model();
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid({0,10},100.0,2.0);
    auto LaxFriedrichs_flux = LaxFriedrichs(grid, model,simulation_time);

    check_consistency(LaxFriedrichs_flux);

}




