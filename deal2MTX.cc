#include <fstream>

int main (int argc, char *argv[])
{

  unsigned n_row (0);
  unsigned n_column (0);
  int n_none_zero (-1);

  {
    std::ifstream fin ("matrix.out");
    while (!fin.eof())
      {
        ++n_none_zero;
        unsigned i,j;
        double ignore;
        fin >> i;
        fin >> j;
        fin >> ignore;

        n_row = std::max (i, n_row);
        n_column = std::max (j, n_column);
      }
    fin.close();
  }


  std::ofstream fout ("matrix.MTX");
  fout <<"%%MatrixMarket matrix coordinate real general\n";
  fout << n_row+1 << "\t" << n_column +1 << "\t"
       << n_none_zero << std::endl;
  std::ifstream fin ("matrix.out");
  for (int line=0; line < n_none_zero; ++line)
    {
      unsigned i,j;
      double value;
      fin >> i;
      fin >> j;
      // i,j transpose
      fout << j + 1 << "\t";
      fout << i + 1 << "\t";
      fin >> value;
      fout << value << std::endl;
    }

  fin.close();
  fout.close();

  return (0);
}
