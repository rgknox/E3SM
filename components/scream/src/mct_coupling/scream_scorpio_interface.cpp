#include "scream_scorpio_interface.hpp"
#include "scream_config.h"

#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_types.hpp"

#include <string>

using scream::Real;
using scream::Int;
extern "C" {

// Fortran routines to be called from C++
  void grid_write_data_array_c_real_1d(const std::string &filename, const std::string &varname, const Int dim1_length, const Real*& hbuf);
  void grid_write_data_array_c_real_2d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int dim2_length, const Real*& hbuf);
  void grid_write_data_array_c_real_3d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Real*& hbuf);
  void grid_write_data_array_c_real_4d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int dim4_length, const Real*& hbuf);
  void grid_write_data_array_c_int_1d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int*& hbuf);
  void grid_write_data_array_c_int_2d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int dim2_length, const Int*& hbuf);
  void grid_write_data_array_c_int_3d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int*& hbuf);
  void grid_write_data_array_c_int_4d(const std::string &filename, const std::string &varname, const Int dim1_length, const Int dim2_length, const Int dim3_length, const Int dim4_length, const Int*& hbuf);
  void eam_init_pio_subsystem_c(const int mpicom, const int compid, const bool local);
  void eam_pio_finalize_c();
  void register_outfile_c(const std::string (&filename));
  void sync_outfile_c(const std::string (&filename));
  void pio_update_time_c(const std::string (&filename),const Real time);
  void register_dimension_c(const std::string &filename, const std::string &shortname, const std::string &longname, const int length);
  void register_variable_c(const std::string& filename,const std::string& shortname, const std::string& longname, const int numdims, const std::string*& var_dimensions, const int dtype, const std::string& pio_decomp_tag);
  void eam_pio_enddef_c();

} // extern C

namespace scream {
namespace scorpio {
/* ----------------------------------------------------------------- */
void eam_init_pio_subsystem(const int mpicom, const int compid, const bool local) {
  eam_init_pio_subsystem_c(mpicom,compid,local);
}
/* ----------------------------------------------------------------- */
void eam_pio_finalize() {
  eam_pio_finalize_c();
}
/* ----------------------------------------------------------------- */
/* ----------------------------------------------------------------- */
void register_outfile(const std::string& filename) {

  register_outfile_c(filename);
}
/* ----------------------------------------------------------------- */
void sync_outfile(const std::string& filename) {

  sync_outfile_c(filename);
}
/* ----------------------------------------------------------------- */
void pio_update_time(const std::string& filename, const Real time) {

  pio_update_time_c(filename,time);
}
/* ----------------------------------------------------------------- */
void register_dimension(const std::string &filename, const std::string& shortname, const std::string& longname, const int length) {

  register_dimension_c(filename, shortname, longname, length);
}
/* ----------------------------------------------------------------- */
void register_variable(const std::string &filename, const std::string& shortname, const std::string& longname, const int numdims, const std::string* var_dimensions, const int dtype, const std::string& pio_decomp_tag) {

  register_variable_c(filename, shortname, longname, numdims, var_dimensions, dtype, pio_decomp_tag);
}
/* ----------------------------------------------------------------- */
void eam_pio_enddef() {
  eam_pio_enddef_c(); 
}
/* ----------------------------------------------------------------- */
//void grid_write_data_array(const std::string &filename, const std::string &varname, const Int buf_length, const Int &hbuf) {
//
//  grid_write_data_array_c_int(filename,varname,buf_length,hbuf);
//
//};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const Int dim1_length, const Real* hbuf) {

  grid_write_data_array_c_real_1d(filename,varname,dim1_length,hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_2d(filename,varname,dim_length[0],dim_length[1],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_3d(filename,varname,dim_length[0],dim_length[1],dim_length[2],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,4>& dim_length, const Real* hbuf) {

  grid_write_data_array_c_real_4d(filename,varname,dim_length[0],dim_length[1],dim_length[2],dim_length[3],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const Int dim1_length, const Int* hbuf) {

  grid_write_data_array_c_int_1d(filename,varname,dim1_length,hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,2>& dim_length, const Int* hbuf) {

  grid_write_data_array_c_int_2d(filename,varname,dim_length[0],dim_length[1],hbuf);

};
/* ----------------------------------------------------------------- */
void grid_write_data_array(const std::string &filename, const std::string &varname, const std::array<Int,3>& dim_length, const Int* hbuf) {

  grid_write_data_array_c_int_3d(filename,varname,dim_length[0],dim_length[1],dim_length[2],hbuf);

};
/* ----------------------------------------------------------------- */

} // namespace scorpio
} // namespace scream
