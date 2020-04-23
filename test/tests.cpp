

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "transfermatrix_test.hpp"

TEST_CASE("testing Trasfer matrix method on 1D heterostructures (GaAs/AlGaAs)") {
    CHECK(test_transfermatrix(20.0) == doctest::Approx(1.70565).epsilon(0.01));
}
