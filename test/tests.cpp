

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "transfermatrix_test.hpp"
#include "transfermatrix_test.hpp"

TEST_CASE("testing Trasfer matrix method on 1D heterostructures (GaAs/AlGaAs)") {
    CHECK(test_transfermatrix(20.0) == doctest::Approx(1.70565).epsilon(0.01));
}

TEST_CASE("testing split-operator method on 1D heterostructures (GaAs/AlGaAs)") {
    CHECK(test_splitoperator(20.0) == doctest::Approx(1.70523).epsilon(0.01));
}

TEST_CASE("testing if split-operator and transfer matrix method are equivalent ( 1% error )") {
    CHECK(test_transfermatrix(20.0) == doctest::Approx(test_splitoperator(20.0)).epsilon(0.01));
}
