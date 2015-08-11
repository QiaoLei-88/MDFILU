#include "MDFILU.h"
#include <fstream>

int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, ::numbers::invalid_unsigned_int);
  const unsigned int degree (10);
  LA::MPI::SparseMatrix system_matrix (degree, degree, /*max_entries_per_row*/degree);
  DynamicMatrix LU (degree, degree, /*max_entries_per_row*/degree);

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
  // MDF_reordering_and_ILU_factoring (system_matrix, LU);
  LU.copy_from (system_matrix);

  // Out put LU
  {
    std::ofstream fout ("LU.out");
    LU.print (fout);
    fout.close();
  }

  // Make preconditioner from LU

  return (0);
}
