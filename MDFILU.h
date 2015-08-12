
#ifndef __MDFILU__H__
#define __MDFILU__H__

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/base/std_cxx11/array.h>

#include <fstream>

#define VERBOSE_OUTPUT

using namespace dealii;

#define USE_TRILINOS_LA
namespace LA
{
#ifdef USE_PETSC_LA
  using namespace dealii::LinearAlgebraPETSc;
#else
  using namespace dealii::LinearAlgebraTrilinos;
#endif
}

#ifndef __NSVector__DEFINED__
typedef LA::MPI::Vector NSVector;
#define __NSVector__DEFINED__
#endif

typedef SparseMatrixEZ<double> DynamicMatrix;
typedef unsigned short local_index_type;
typedef unsigned long global_index_type;
typedef double data_type;
typedef bool flag_type;

#define N_INDICATOR 3
class Indicator: public std_cxx11::array<data_type,N_INDICATOR>
{
public:
  void init();
  int operator- (const Indicator &op) const;
};

template<typename Matrix>
void get_indices_of_non_zeros (
  const Matrix &matrix,
  const global_index_type row_to_factor,
  const std::vector<flag_type> &row_factored,
  std::vector<global_index_type> &incides_need_update,
  const bool except_pivot);

void compute_discarded_value (
  const unsigned int row_to_factor,
  const DynamicMatrix &LU,
  const DynamicMatrix &fill_in_level,
  const std::vector<flag_type> &row_factored,
  const unsigned int fill_in_threshold,
  Indicator &return_value);

global_index_type find_min_discarded_value (
  const std::vector<Indicator> &indicators,
  const std::vector<flag_type> &row_factored);

void MDF_reordering_and_ILU_factoring (
  const LA::MPI::SparseMatrix &system_matrix,
  DynamicMatrix &LU,
  std::vector<global_index_type> &permutation);

#endif
