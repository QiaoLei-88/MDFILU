#include "MDFILU.h"
#include <fstream>



int main (int argc, char *argv[])
{
// #ifdef VERBOSE_OUTPUT
//   debugStream.open("debug.out");
// #endif

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, ::numbers::invalid_unsigned_int);
  const global_index_type degree (10);
  const unsigned int estimated_row_length (10);
  LA::MPI::SparseMatrix system_matrix (degree, degree, /*max_entries_per_row*/estimated_row_length);
  DynamicMatrix LU (degree, degree, /*max_entries_per_row*/estimated_row_length);
  std::vector<global_index_type> permutation (degree);

  // Set value for system_matrix
  std::ifstream fin ("matrix.dat");
  for (global_index_type i=0; i<degree; ++i)
    for (global_index_type j=0; j<degree; ++j)
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
    const data_type tolerance (1e-12);
    DynamicMatrix A (degree, degree, estimated_row_length);
    // Compute LD*U
    for (global_index_type i=0; i<degree; ++i)
      {
        const global_index_type i_row = permutation[i];
        for (global_index_type j=0; j<degree; ++j)
          {
            const global_index_type j_col = permutation[j];
            data_type value = 0;
            global_index_type vmult_max_index = 0;
            // Diagonal values of L is always 1 thus not been stored.
            // Recover its effect manually.
            if (j>=i)
              {
                value = LU.el (i_row,j_col);
                vmult_max_index = i;
              }
            else
              {
                vmult_max_index = j+1;
              }
            for (global_index_type k=0; k<vmult_max_index; ++k)
              {
                const global_index_type k_permuted = permutation[k];
                value += LU.el (i_row,k_permuted) * LU.el (k_permuted,j_col);
              }
            if (std::abs (value) > tolerance)
              {
                A.set (i_row, j_col, value);
              }
          }
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
