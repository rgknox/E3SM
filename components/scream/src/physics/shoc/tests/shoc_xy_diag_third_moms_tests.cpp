#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestXYdiagThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   x_y_terms_diag_third_shoc_moment
  
    //  Test to be sure that very high values of w3
    //    are reduced but still the same sign   
  

    // Initialize data structure for bridging to F90
    SHOCXYdiagthirdmomsData SDS;

    
    // Call the fortran implementation    
    x_y_terms_diag_third_shoc_moment(SDS);    
    
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_xy_diag_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestXYdiagThirdMoms;
  
  TestStruct::run_property();
}

TEST_CASE("shoc_xy_diag_third_moms_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestXYdiagThirdMoms;

}

} // namespace
