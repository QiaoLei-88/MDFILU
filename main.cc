#include "MDFILU.h"
#include <fstream>



int main (int argc, char *argv[])
{
// #ifdef VERBOSE_OUTPUT
//   debugStream.open("debug.out");
// #endif

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, ::numbers::invalid_unsigned_int);
  const unsigned int degree (10);
  LA::MPI::SparseMatrix system_matrix (degree, degree, /*max_entries_per_row*/degree);
  DynamicMatrix LU (degree, degree, /*max_entries_per_row*/degree);
  std::vector<global_index_type> permutation (degree);

  // Set value for system_matrix
  std::ifstream fin ("matrix.dat");
  for (unsigned i=0; i<degree; ++i)
    for (unsigned j=0; j<degree; ++j)
      {
        data_type value;
        fin >> value;
        system_matrix.set (i,j,value);
      }
  fin.close();
  system_matrix.compress (VectorOperation::insert);
  {
    // Out put system_matrix
    std::ofstream fout ("matrix.out");
    system_matrix.print (fout);
    fout.close();
  }
  MDF_reordering_and_ILU_factoring (system_matrix, LU, permutation);

  // Out put LU
  {
    std::ofstream fout ("LU.out");
    LU.print (fout);
    fout.close();
  }

  // Test the program by multiply LU back
  // Ignore sparsity pattern fist.
  {
    const unsigned int estimated_row_length (10);
    DynamicMatrix A (degree, degree, estimated_row_length);
    // Compute LD*U
    for (unsigned i_row=0; i_row<degree; ++i_row)
      for (unsigned j_col=0; j_col<degree; ++j_col)
        {
          data_type value = 0;
          unsigned vmult_max_index = 0;
          // Diagonal values of L is always 1 thus not been stored.
          // Recover its effect manually.
          if (j_col>=i_row)
            {
              value = LU.el (i_row,j_col);
              vmult_max_index = i_row;
            }
          else
            {
              vmult_max_index = j_col+1;
            }
          for (unsigned k=0; k<vmult_max_index; ++k)
            {
              value += LU.el (i_row,k) * LU.el (k,j_col);
            }
          A.set (i_row, j_col, value);
        }
    std::ofstream fout ("A.out");
    A.print (fout);
    fout.close();
  }
  // Make preconditioner from LU
// #ifdef VERBOSE_OUTPUT
//   debugStream.close();
// #endif
  return (0);
}
