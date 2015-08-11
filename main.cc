#include "MDFILU.h"


int main()
{
  const unsigned int degree (10);
  LA::MPI::SparseMatrix system_matrix (degree, degree);
  DynamicMatrix LU (degree, degree);

  // Set value for system_matrix

  // Out put system_matrix

  MDF_reordering_and_ILU_factoring (system_matrix, LU);

  // Out put LU

  // Make preconditioner from LU

  return (0);
}
