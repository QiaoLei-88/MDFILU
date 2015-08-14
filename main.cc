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
  SourceMatrix system_matrix (degree, degree, /*max_entries_per_row*/estimated_row_length);

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
  MDFILU mdfilu (system_matrix, estimated_row_length, 20);

  // Out put LU
  {
    std::ofstream fout ("LU.out");
    mdfilu.get_LU().print (fout);
    fout.close();
  }

  // Test the MDFILU class by multiplying LU back
  // For now the sparsity pattern is ignored.
  {
    const data_type tolerance (1e-12);
    DynamicMatrix A (degree, degree, estimated_row_length);
    // Compute LD*U
    for (global_index_type i=0; i<degree; ++i)
      {
        const global_index_type i_row = mdfilu.get_permutation()[i];
        for (global_index_type j=0; j<degree; ++j)
          {
            const global_index_type j_col = mdfilu.get_permutation()[j];
            data_type value = 0;
            global_index_type vmult_max_index = 0;
            // Diagonal values of L is always 1 thus not been stored.
            // Recover its effect manually.
            if (j>=i)
              {
                value = mdfilu.get_LU().el (i_row,j_col);
                vmult_max_index = i;
              }
            else
              {
                vmult_max_index = j+1;
              }
            for (global_index_type k=0; k<vmult_max_index; ++k)
              {
                const global_index_type k_permuted = mdfilu.get_permutation()[k];
                value += mdfilu.get_LU().el (i_row,k_permuted) * mdfilu.get_LU().el (k_permuted,j_col);
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

// #ifdef VERBOSE_OUTPUT
//   debugStream.close();
// #endif
  return (0);
}
