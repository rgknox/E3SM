#include "catch2/catch.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"
#include "diagnostics/potential_temperature.hpp"
#include "diagnostics/register_diagnostics.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/field/field_utils.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

//-----------------------------------------------------------------------------------------------//
template<typename DeviceT>
void run(std::mt19937_64& engine, int int_ptype)
{
  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<Real>;
  using Pack       = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<Pack>;
  using rview_1d   = typename KT::template view_1d<Real>;

  const     int packsize = SCREAM_PACK_SIZE;
  constexpr int num_levs = packsize*2 + 1; // Number of levels to use for tests, make sure the last pack can also have some empty slots (packsize>1).
  const     int num_mid_packs = ekat::npack<Pack>(num_levs);

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // Create a grids manager - single column for these tests
  const int ncols = 1;
  auto gm = create_gm(comm,ncols,num_levs);

  // Kokkos Policy
  auto policy = ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(ncols, num_mid_packs);

  // Input (randomized) views
  view_1d temperature("temperature",num_mid_packs),
          pressure("pressure",num_mid_packs),
          condensate("condensate",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<Real*>(v.data()),v.size()*packsize);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Construct the Diagnostic
  ekat::ParameterList params;
  register_diagnostics();
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  std::string ptype = int_ptype == 0 ? "Tot" : "Liq";
  params.set("Temperature Kind", ptype);
  auto diag = diag_factory.create("PotentialTemperature",comm,params);
  diag->set_grids(gm);


  // Set the required fields for the diagnostic.
  std::map<std::string,Field> input_fields;
  for (const auto& req : diag->get_required_field_requests()) {
    Field f(req.fid);
    auto & f_ap = f.get_header().get_alloc_properties();
    f_ap.request_allocation(packsize);
    f.allocate_view();
    const auto name = f.name();
    f.get_header().get_tracking().update_time_stamp(t0);
    diag->set_required_field(f.get_const());
    REQUIRE_THROWS(diag->set_computed_field(f));
    input_fields.emplace(name,f);
  }

  // Initialize the diagnostic
  diag->initialize(t0,RunType::Initial);

  // Run tests
  {
    // Construct random data to use for test
    // Get views of input data and set to random values
    const auto& T_mid_f = input_fields["T_mid"];
    const auto& T_mid_v = T_mid_f.get_view<Pack**>();
    const auto& p_mid_f = input_fields["p_mid"];
    const auto& p_mid_v = p_mid_f.get_view<Pack**>();
    const auto& q_mid_f = input_fields["qc"];
    const auto& q_mid_v = q_mid_f.get_view<Pack**>();
    for (int icol=0;icol<ncols;icol++) {
      const auto& T_sub = ekat::subview(T_mid_v,icol);
      const auto& p_sub = ekat::subview(p_mid_v,icol);
      const auto& q_sub = ekat::subview(q_mid_v,icol);
      ekat::genRandArray(dview_as_real(temperature), engine, pdf_temp);
      ekat::genRandArray(dview_as_real(pressure),    engine, pdf_pres);
      ekat::genRandArray(dview_as_real(condensate),  engine, pdf_pres);
      Kokkos::deep_copy(T_sub,temperature);
      Kokkos::deep_copy(p_sub,pressure);
      Kokkos::deep_copy(q_sub,condensate);
    }

    // Run diagnostic and compare with manual calculation
    diag->compute_diagnostic();
    const auto& diag_out = diag->get_diagnostic();
    Field theta_f = diag_out.clone();
    theta_f.deep_copy<double,Host>(0.0);
    theta_f.sync_to_dev();
    const auto& theta_v = theta_f.get_view<Pack**>();
    Kokkos::parallel_for("", policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int icol = team.league_rank();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team,num_mid_packs), [&] (const Int& jpack) {
        auto theta = PF::calculate_theta_from_T(T_mid_v(icol,jpack),p_mid_v(icol,jpack));
        if (int_ptype==1) {
          theta_v(icol,jpack) = theta - (theta/T_mid_v(icol,jpack)) * (PC::LatVap/PC::Cpair) * q_mid_v(icol,jpack);
        } else { theta_v(icol,jpack) = theta; }
      });
      team.team_barrier();
    });
    Kokkos::fence();
    REQUIRE(views_are_equal(diag_out,theta_f));
  }
 
  // Finalize the diagnostic
  diag->finalize(); 

} // run()

TEST_CASE("potential_temp_test", "potential_temp_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    for (const int int_ptype : {0,1}) {
    run<Device>(engine, int_ptype);
  }
  }
  printf("ok!\n");

  printf("\n");

} // TEST_CASE

} // namespace
