#include <iostream>
#include "eval.hpp"
#include "png.hpp"
using namespace std;
using namespace middlebury;

int main(int argc, char *argv[]) {

  std::string path_disp, path_gt, path_mask;
  if(argc == 4){
    path_disp = argv[1];
    path_gt = argv[2];
    path_mask = argv[3];
  }
  else{
    cout << "usage: ./eval disparity disparity_gt mask" << endl;
    return 1;
  }

  Image mask = readPNG(path_mask);
  int rows_disp, rows_gt, cols_disp, cols_gt;

  // Load pfm images
  float* disp = ReadFilePFM(path_disp.c_str(), &rows_disp, &cols_disp);
  float* gt = ReadFilePFM(path_gt.c_str(), &rows_gt, &cols_gt);

  int maxdisp = 64;
  int rounddisp = 1;

  StereoScore score = evaldisp(disp, gt, mask, rows_disp, cols_disp, maxdisp, rounddisp);
  cout << setw(12) << "coverage"
       << setw(12) << "bad-0.5"
       << setw(12) << "bad-1"
       << setw(12) << "invalid"
       //<< setw(12) << "total-bad"
       << setw(12) << "avg-err" << endl;
  cout << setw(12) << score.coverage
       << setw(12) << score.bad05
       << setw(12) << score.bad1
       << setw(12) << score.invalid
       << setw(12) << score.avgErr << endl;

  delete[] disp;
  delete[] gt;

  return 0;
}
