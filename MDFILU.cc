
#include "MDFILU.h"

void Indicator::init()
{
  for (unsigned int i=0; i<N_INDICATOR; ++i)
    {
      (*this)[i] = 0.0;
    }
  return;
}

int Indicator::operator- (const Indicator &op) const
{
  for (unsigned int i=0; i<N_INDICATOR; ++i)
    {
      if ((*this)[i] < op[i])
        {
          return (-1);
        }
      if ((*this)[i] > op[i])
        {
          return (1);
        }
    }
  return (0);
}

template<typename Matrix>
void get_indices_of_non_zeros (
  const Matrix &matrix,
  const global_index_type row_to_factor,
  const std::vector<flag_type> &row_factored,
  std::vector<global_index_type> &incides_need_update,
  const bool except_pivot)
{
  global_index_type i_col = 0;
  for (typename Matrix::const_iterator iter_col = matrix.begin (row_to_factor);
       iter_col != matrix.end (row_to_factor);
       ++iter_col, ++i_col)
    {
      const global_index_type i_col = iter_col->column();
      if (i_col == row_to_factor && except_pivot)
        {
          // If we do not want to count on pivot, jump over
          continue;
        }

      if (!row_factored[i_col])
        {
          incides_need_update.push_back (i_col);
        }
    }
  return;
}


void compute_discarded_value (
  const unsigned int row_to_factor,
  const DynamicMatrix &LU,
  const DynamicMatrix &fill_in_level,
  const std::vector<flag_type> &row_factored,
  const unsigned int fill_in_threshold,
  Indicator &return_value)
{
#ifdef VERBOSE_OUTPUT
  if (row_to_factor == 0)
    {
      std::ofstream f_level ("fill_level_cdv.out");
      fill_in_level.print (f_level);
      f_level.close();
    }
#endif
  const int prior_discarded_value (0);
  const int prior_n_discarded (1);
  const int prior_n_fill (2);
  // const int prior_index(3);

  return_value.init();
  const data_type pivot = LU.el (row_to_factor,row_to_factor);

  if (pivot==0.0)
    {
      return_value[prior_n_fill] = 1<<30;
      return_value[prior_n_discarded] = 1<<30;
      return_value[prior_discarded_value] = 1e+200;
      return;
    }

  // compute discarded value for i_row := row_to_factor.
  // During this procedure, we need to go through all un-factored entries that connected
  // with this row, i.e., for all k that a(i_row, k) \ne 0 and a(k, i_row) \ne 0.
  // That's why we need the flag array row_factored.

  // Find number of rows need to go through. The value is all non-zero
  // entries except the pivot and factored rows.
  std::vector<global_index_type> incides_need_update;
  const bool except_pivot (true);
  get_indices_of_non_zeros<DynamicMatrix>
  (LU, row_to_factor, row_factored, incides_need_update, except_pivot);

  const data_type pivot_neg_inv = -1.0/pivot;
  const global_index_type n_row_need_update = incides_need_update.size();

  for (global_index_type i=0; i<n_row_need_update; ++i)
    {
      const global_index_type i_row = incides_need_update[i];
      for (global_index_type j=0; j<n_row_need_update; ++j)
        {
          const global_index_type j_col = incides_need_update[j];
          // Check fill-in level
          const data_type new_fill_in_level = fill_in_level.el (i_row, j_col);
#ifdef VERBOSE_OUTPUT
          std::ofstream f_level ("fill_level_verbos.out", std::fstream::app);
          f_level << row_to_factor << "\t"
                  << i_row << "\t" << j_col << "\t"
                  << fill_in_level.el (i_row, j_col) << std::endl;
          f_level.close();
#endif
          if (new_fill_in_level == 0.0 /* fill in level for new entry*/)
            {
              ++return_value[prior_n_fill];
            }

          // Make sure that the provided fill_in_threshold consists with
          // the internal definition, i.e., has a offset one. See documentation
          // above for details
          if (new_fill_in_level > fill_in_threshold)
            {
              // Element will be discarded
              const data_type update = pivot_neg_inv *
                                       LU.el (row_to_factor,j_col) * LU.el (i_row,row_to_factor);
              return_value[prior_discarded_value] += update*update;
              ++return_value[prior_n_discarded];
            }
        } // For each column need update
    } // For each row need update

  return;
}

// Determine the next row to be factored by finding out the one with minimum
// indicator form rows that have not been factored.
global_index_type find_min_discarded_value (
  const std::vector<Indicator> &indicators,
  const std::vector<flag_type> &row_factored)
{
  global_index_type candidate (0);
  bool need_init_candidate (true);
  for (global_index_type i=0; i<indicators.size(); ++i)
    {
      if (row_factored[i])
        {
          continue;
        }
      // Set first un-factored row as candidate if it is not initialized
      if (need_init_candidate)
        {
          candidate = i;
          need_init_candidate = false;
        }

      // This means will not update selected row when indicators is equal,
      // which means small index is with higher priority in tie situation.
      // One can switch to large index first by using " <=0 ".
      if ((indicators[i] - indicators[candidate]) < 0)
        {
          candidate = i;
        }
    }
  return (candidate);
}


void MDF_reordering_and_ILU_factoring (
  const LA::MPI::SparseMatrix &system_matrix,
  DynamicMatrix &LU,
  std::vector<global_index_type> &permutation)
{
#ifdef VERBOSE_OUTPUT
  std::ofstream debugStream ("debug.out");
#endif
  const global_index_type fill_in_threshold (20);
  const global_index_type degree = system_matrix.m();
  // Record fill-in level for all non-zero entries, we need this to compute
  // level for new fill-ins.
  // SparseMatrixEZ<double> fill_in_level (degree,degree,degree);
  DynamicMatrix fill_in_level (degree,degree,degree);

  //initialize the LU matrix
  LU.copy_from (system_matrix);

  std::vector<Indicator> indicators (degree);
  std::vector<flag_type> row_factored (degree, false);

  // Initialize::BEGIN
  // Compute Initial fill in level
  for (global_index_type i_row=0; i_row<degree; ++i_row)
    {
      // Get indices of non-zero's of the target row

      std::vector<global_index_type> incides_of_non_zeros;
      get_indices_of_non_zeros<DynamicMatrix>
      (LU, i_row, row_factored, incides_of_non_zeros, /*except_pivot=*/ false);

      const global_index_type n_non_zero_in_row = incides_of_non_zeros.size();
      // global_index_type i_col = 0;
      // for (LA::MPI::SparseMatrix::const_iterator iter_col (system_matrix.begin (i_row));
      //      iter_col != system_matrix.end (i_row);
      //      ++iter_col, ++i_col)
      //   {
      //     incides_of_non_zeros[i_col] = iter_col->column();
      //   }

      // Set initial fill in level.
      // Because the matrix access function SparseMatrixEZ.el(i,j) will return a zero when element
      // (i,j) is out of the sparsity pattern, here I set the fill-in level with offset 1.
      // That is to say when fill-in level equals
      //    0  :  level infinite, new fill-in
      //    1  :  level 0 in article, original entry
      //    2  :  level 1 fill in
      //    ... so on the same.

      // for (global_index_type j_col=0; j_col<system_matrix.size(); ++j_col)
      for (global_index_type i_nz=0; i_nz<n_non_zero_in_row; ++i_nz)
        {
          fill_in_level.set (i_row, incides_of_non_zeros[i_nz], 1);
          // // a(i,j) exists?
          // if (system_matrix. (i,j) == 0.0)
          //   {
          //     fill_in_level.add (i, j, numbers::invalid_unsigned_int);
          //   }
          // else
          //   {
          //     fill_in_level.add (i, j, 0);
          //   }
        }
    }

  // Compute initial discarded value, must be done after all fill-in level have
  // been set.
  for (global_index_type i_row=0; i_row<degree; ++i_row)
    {
      compute_discarded_value (i_row,
                               LU,
                               fill_in_level,
                               row_factored,
                               fill_in_threshold,
                               indicators[i_row]);
    }

  // Initialize::END
#ifdef VERBOSE_OUTPUT
  {
    std::ofstream f_level ("fill_level.out");
    fill_in_level.print (f_level);
    f_level.close();
  }
#endif

  // Factoring the matrix
  for (global_index_type n_row_factored=0; n_row_factored<degree; ++n_row_factored)
    {
      // Find the row with minimal discarded value
      const global_index_type row_to_factor
        = find_min_discarded_value (indicators,row_factored);

#ifdef VERBOSE_OUTPUT
      for (global_index_type i=0; i<degree; ++i)
        {
          debugStream << indicators.at (i).at (0) << "\t ";
        }
      debugStream << std::endl;
      for (global_index_type i=0; i<degree; ++i)
        {
          debugStream << indicators.at (i).at (1) << "\t ";
        }
      debugStream << std::endl;
      for (global_index_type i=0; i<degree; ++i)
        {
          debugStream << indicators.at (i).at (2) << "\t ";
        }
      debugStream << std::endl;
#endif

      row_factored[row_to_factor] = true;
      permutation[n_row_factored] = row_to_factor;

#ifdef VERBOSE_OUTPUT
      debugStream << "row_to_factor: " << row_to_factor << std::endl;
#endif

      // Find number of rows need to go through. The value is all non-zero
      // entries except the pivot and factored rows.
      std::vector<global_index_type> incides_need_update;
      const data_type pivot = LU.el (row_to_factor,row_to_factor);
      Assert (pivot != 0.0, ExcMessage ("Zero pivot encountered!"));
      const bool except_pivot (true);
      get_indices_of_non_zeros<DynamicMatrix>
      (LU, row_to_factor, row_factored, incides_need_update, except_pivot);

      const data_type pivot_inv = 1.0/pivot;
      const global_index_type n_row_need_update = incides_need_update.size();

      for (global_index_type i=0; i<n_row_need_update; ++i)
        {
          const global_index_type i_row = incides_need_update[i];

          // Update current column, i.e., lower triangle part
          const data_type value = LU.el (i_row, row_to_factor);
          LU.set (i_row,row_to_factor, value*pivot_inv);

          // Update the remaining matrix
          for (global_index_type j=0; j<n_row_need_update; ++j)
            {
              const global_index_type j_col = incides_need_update[j];
              // Check fill-in level
              unsigned int new_fill_in_level = fill_in_level.el (i_row, j_col);
              if (new_fill_in_level == 0 /* fill in level for new entry*/)
                {
                  // Because we set original entry to level 1, what we need to set is
                  // (lv(i)-1) + (lv(j)-1) + 1 + 1 = lv(i) + lv(j).
                  new_fill_in_level = fill_in_level.el (row_to_factor,j_col) +
                                      fill_in_level.el (i_row,row_to_factor);
                }

              const unsigned int fill_in_threshold (3);
              // Make sure that the provided fill_in_threshold consists with
              // the internal definition, i.e., has a offset one. See documentation
              // above for details
              if (new_fill_in_level <= fill_in_threshold)
                {
                  // Element accepted
                  const data_type update = -LU.el (row_to_factor,j_col) * LU.el (i_row,row_to_factor);
                  LU.add (i_row, j_col, update);
                  // Update fill-level if this is a new entry
                  fill_in_level.set (i_row, j_col, new_fill_in_level);
                }
            } // For each column need update
          compute_discarded_value (i_row,
                                   LU,
                                   fill_in_level,
                                   row_factored,
                                   fill_in_threshold,
                                   indicators[i_row]);
        } // For each row need update
    } // For each row in matrix
  return;
}

// template get_indices_of_non_zeros<LA::MPI::SparseMatrix>;
template
void get_indices_of_non_zeros<DynamicMatrix> (
  const DynamicMatrix &matrix,
  const global_index_type row_to_factor,
  const std::vector<flag_type> &row_factored,
  std::vector<global_index_type> &incides_need_update,
  const bool except_pivot);
