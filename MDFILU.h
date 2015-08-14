
#ifndef __MDFILU__H__
#define __MDFILU__H__

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/std_cxx11/array.h>
#include <Epetra_Operator.h>
#include <Epetra_Map.h>
#include <Sacado.hpp>


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
typedef LA::MPI::SparseMatrix SourceMatrix;

// #define UseTrilinosMatrix
#ifdef UseTrilinosMatrix
typedef LA::MPI::SparseMatrix DynamicMatrix;
#else
typedef SparseMatrixEZ<double> DynamicMatrix;
#endif

typedef unsigned short local_index_type;
typedef unsigned long global_index_type;
typedef double data_type;
typedef Vector<data_type> MDFVector;
typedef bool flag_type;

class MDFILU : public Epetra_Operator
{
public:
  MDFILU (const SourceMatrix &matrix,
          const global_index_type estimated_row_length_in,
          const global_index_type fill_in_threshold_in);

  int apply (const data_type *const in, data_type *const out) const;
  int apply_inverse (const data_type *const in, data_type *const out) const;

  int apply (const MDFVector &in, MDFVector &out) const;
  int apply_inverse (const MDFVector &in, MDFVector &out) const;

  // Interface for debug code
  const std::vector<global_index_type> &get_permutation() const;
  const DynamicMatrix &get_LU() const;
  ~MDFILU();

  // Virtual functions from Epetra_Operator

  virtual int Apply (const Epetra_MultiVector &, Epetra_MultiVector &) const;
  virtual int ApplyInverse (const Epetra_MultiVector &, Epetra_MultiVector &) const;
  virtual double NormInf() const;
  virtual bool UseTranspose() const;
  virtual bool HasNormInf() const;

  virtual const char *Label() const;

  virtual const Epetra_Comm &Comm() const;
  virtual const Epetra_Map &OperatorDomainMap() const;
  virtual const Epetra_Map &OperatorRangeMap() const;

  virtual int SetUseTranspose (const bool);

private:
  const global_index_type invalid_index;
  const data_type very_large_number;

#define N_INDICATOR 3
  class Indicator: public std_cxx11::array<data_type,N_INDICATOR>
  {
  public:
    void init();
    int operator- (const Indicator &op) const;
  };

  void get_indices_of_non_zeros (
    const global_index_type row_to_factor,
    std::vector<global_index_type> &incides_need_update,
    const bool except_pivot) const;

  void compute_discarded_value (const unsigned int row_to_factor);

  global_index_type find_min_discarded_value() const;
  void MDF_reordering_and_ILU_factoring();


  const global_index_type degree;
  const global_index_type estimated_row_length_in;
  const global_index_type fill_in_threshold;
  DynamicMatrix LU;
  // Record fill-in level for all non-zero entries, we need this to compute
  // level for new fill-ins.
  // SparseMatrixEZ<double> fill_in_level (degree,degree,degree);
  DynamicMatrix fill_in_level;

  // Where is the k-th row and column in LU
  std::vector<global_index_type> permute_logical_to_storage;
  // For k-th row and column in LU, where is it in the permuted matrix
  std::vector<global_index_type> permuta_storage_to_logical;

  std::vector<Indicator> indicators;

  // During factoring procedure, we need to go through all un-factored entries that connected
  // with this row, i.e., for all k that a(i_row, k) \ne 0 and a(k, i_row) \ne 0.
  // That's why we need the flag array row_factored.
  std::vector<flag_type> row_factored;


  // Data for Epetra_Operator interface
  bool use_transpose;
  const bool has_norm_infty;
  const Epetra_Comm *epetra_comm;
  const Epetra_Map operator_domain_map;
  const Epetra_Map operator_range_map;

  const static char label[];
};

#endif
//     of #ifndef __MDFILU__H__
