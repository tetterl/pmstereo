#ifndef HELPERS_HPP
#define HELPERS_HPP
/// Common functions used by the algorithm, independet of implementation
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string.h>
#include <vector>

#include "parameters.hpp"
#include "png.hpp"

using namespace std;

// Global random number generator used for intitialization, fixed random seed
std::mt19937 gen(0);

/**
 * @brief The View struct
 *
 * Common view struct with maximum of information.
 * Used for common functions, used for intitialization
 * of specific implementations
 */
struct CommonView{
  // All data in this struct is stored in row order
  // Image in HWC format i.e. row order and channels interleaved (R1G1B1,R2G2B2, ...)
  uint8_t* i;
  // Gradient
  float* g;
  // Planes (format: ABCN1N2N3P1P2P3)
  // ABC: plane coeffs
  // N1N2N3 :  normal coeffs
  // P1P2P3 :  point coeffs
  float* p;
  // disparity
  float* d;
  // cost
  float* c;
  // disparity(used in post processing)
  uint8_t* v;
};

/**
 * @brief allocate_common_view
 * @param v
 * @param rows
 * @param cols
 */
void allocate_common_view(CommonView& v, int rows, int cols) {
  // Images
  v.i = (uint8_t*)malloc(rows*cols*3);
  // Planes
  v.p = (float*)malloc(rows*cols*9*sizeof(float));
  v.g = (float*)malloc(rows*cols*sizeof(float));
  // Cost
  v.c = (float*)malloc(rows*cols*sizeof(float));
  // Disparity
  v.d = (float*)malloc(rows*cols*sizeof(float));
  // Disparity Validity
  v.v = (uint8_t*)malloc(rows*cols);
}

/**
 * @brief deallocate_common_view
 * @param v
 */
void deallocate_common_view(CommonView& v) {
  delete[] v.i;
  delete[] v.p;
  delete[] v.g;
  delete[] v.d;
  delete[] v.v;
  delete[] v.c;
}


// Modified from Middlebury benchmark code
// Stores floating point image as .pfm file

void WriteFilePFM(CommonView& v, const char* filename, int rows, int cols, float scalefactor=1/255.0) {
    // Write a PFM file
    int nBands = 1;

    // Open the file
    FILE *stream = fopen(filename, "wb");
    if (stream == 0){
      cout << "ERR: could not open file for writing" << endl;
      return ;
    }

    // sign of scalefact indicates endianness, see pfms specs
    //if (littleendian())
    scalefactor = -scalefactor; //since we are using x86, yes it is.

    // write the header: 3 lines: Pf, dimensions, scale factor (negative val == little endian)
    fprintf(stream, "Pf\n%d %d\n%f\n", cols, rows, scalefactor);

    int n = cols;
    // write rows -- pfm stores rows in inverse order!
    for (int y =rows-1; y >= 0; y--) {
      float* ptr = (float *)&v.d[y*cols];
      if ((int)fwrite(ptr, sizeof(float), n, stream) != n){
        cout << "ERR: File is too short" << endl;
        return;
      }
    }

    // close file
    if (fclose(stream)){
        cout << "ERR: error closing file" << endl;
        return ;
    }
}
/**
 * @brief Converts an image from the png.hpp header into a simple array
 *
 * @param iimg  The input image
 * @param oimg  The output image
 */
void convert_img_to_array(const Image& iimg, uint8_t* oimg){
  for(int i=0;i<iimg.width*iimg.height*iimg.channels;i++){
     oimg[i] = (uint8_t)iimg.data[i];
  }
}

/**
 * @brief Converts an RGB image into black and white
 *
 * @param in    The input image
 * @param out   The output image
 * @param rows  Rows
 * @param cols  Cols
 */
void rgb_to_bw(uint8_t* in, uint8_t* out, int rows, int cols){
  for(int i=0;i<rows*cols;i++){

    out[i] = 0.3*in[3*i] + 0.59*in[3*i+1] + 0.11*in[3*i+2];
    //out[i] = ( in[3*i] + in[3*i+1] + in[3*i+2] ) / 3;
  }
}

/**
 * @brief Naive sobel filter. Implemented as shown
 *        in https://en.wikipedia.org/wiki/Sobel_operator,
 *
 * @param in        Input image
 * @param gradient  Gradient output
 * @param rows      rows
 * @param cols      columns
 */
void sobelNaive(uint8_t* in, float *gradient, int rows, int cols) {
  //assert( width % 16 == 0 && "width must be multiple of 16!" );
  uint8_t* ptr = in;

  uint8_t* p11 = ptr + 0*cols;
  uint8_t* p12 = ptr + 0*cols + 1;
  uint8_t* p13 = ptr + 0*cols + 2;

  uint8_t* p21 = ptr + 1 * cols;
  uint8_t* p22 = ptr + 1 * cols + 1;
  uint8_t* p23 = ptr + 1 * cols + 2;

  uint8_t* p31 = ptr + 2 * cols ;
  uint8_t* p32 = ptr + 2 * cols + 1;
  uint8_t* p33 = ptr + 2 * cols + 2;

  // output pointer
  float* optr = gradient + 1 * cols + 1;
  // Apply 3x3 sobel filter to image less pixel border of 1 (to avoid treating boundary) (unoptimized)
  for (int iy = 1; iy <rows - 1; iy++) {
    for (int ix = 0; ix < cols ; ix++) {
      int sx = ( *p11 + *p31 + 2 * *p21 - *p13 - 2 * *p23  - *p33) ;
      int sy = ( *p11 + *p13 + 2 * *p12 - *p31 - 2 * *p32  - *p33) ;
      float g = sqrt(sx * sx + sy * sy);
      *optr = g; //(int(g) < 255.f) ? g : 255.f;
      p11++; p12++; p13++; p21++; p22++; p23++; p31++; p32++; p33++; optr++;
    }
  }
}

void compute_gradients(CommonView& v, int rows, int cols) {
  uint8_t* ibw = (uint8_t*)malloc(rows*cols*sizeof(float));
  rgb_to_bw(v.i,ibw,rows,cols);
  sobelNaive(ibw,v.g,rows, cols);
  delete[] ibw;
}

/**
 * @brief Random plane initialization as described in the pmstereo paper
 *
 * @param planes  Pointer to a plane array
 * @param rows    Number of rows
 * @param cols    Number of cols
 */
void randomInit(float* planes, int rows, int cols){
  std::uniform_real_distribution<float> z_dis(0.f, float(max_disp));
  std::uniform_real_distribution<float> n_dis(-1.f, 1.f);
  for(int y=0;y<rows;y++){
    for(int x=0;x<cols;x++){
      //uniform z in [0,MAX_DISP)
      float z0 = z_dis(gen);
      //random normal
      float nx = n_dis(gen);
      float ny = n_dis(gen);
      float nz = n_dis(gen);
      //make unit length
      float n = sqrt(nx*nx + ny*ny + nz * nz);
      nx/=n;
      ny/=n;
      nz/=n;
      //get Plane params
      float a = -nx/nz;
      float b = -ny/nz;
      float c = (nx*x+ny*y+nz*z0)/nz;
      planes[(y*cols*9)+x*9]   = a;
      planes[(y*cols*9)+x*9+1] = b;
      planes[(y*cols*9)+x*9+2] = c;
      planes[(y*cols*9)+x*9+3] = nx;
      planes[(y*cols*9)+x*9+4] = ny;
      planes[(y*cols*9)+x*9+5] = nz;
      planes[(y*cols*9)+x*9+6] = x;
      planes[(y*cols*9)+x*9+7] = y;
      planes[(y*cols*9)+x*9+8] = z0;
    }
  }
}

/**
 * @brief fill_invalid_pixels
 * @param v
 * @param x
 * @param y
 * @param rows
 * @param cols
 */
void fill_invalid_pixels(CommonView& v, int x, int y, int rows, int cols){
  int x_lft = x - 1;
  int x_rgt = x + 1;

  // If left neighbor invalid
  //while(!validity(y, x_lft) && x_lft >= 0)
  while(v.v[y*cols+x_lft]==0 && x_lft >= 0)
    --x_lft;

  // If right neighbor invalid
  //while(!validity(y, x_rgt) && x_lft < cols)
  while(v.v[y*cols+x_rgt]==0 && x_rgt < cols)
    ++x_rgt;

  int best_plane_x = x;

  if(x_lft >= 0 && x_rgt < cols)
  {
    float disp_l = v.p[9*(y*cols+x_lft)]*x+v.p[9*(y*cols+x_lft)+1]*y+v.p[9*(y*cols+x_lft)+2];
    float disp_r = v.p[9*(y*cols+x_rgt)]*x+v.p[9*(y*cols+x_rgt)+1]*y+v.p[9*(y*cols+x_rgt)+2];

    best_plane_x = (disp_l < disp_r) ? x_lft : x_rgt;
  }
  else if(x_lft >= 0)
    best_plane_x = x_lft;
  else if(x_rgt < cols)
    best_plane_x = x_rgt;

  //planes(y, x) = planes(y, best_plane_x);
  float* plane_old = &v.p[(y*cols*9)+x*9];
  float* plane_new = &v.p[(y*cols*9)+best_plane_x*9];
  memcpy(plane_old,plane_new,9*sizeof(float));
}

/**
 * @brief Converts the planes contained in a view v into disparities
 *
 * @param v     a view struct
 * @param rows  total rows
 * @param cols  total columns
 */
void PlanesToDisparity(CommonView& v, int rows, int cols){
  int idx=0;
  for(int y=0; y < rows; ++y)
    for(int x=0; x < cols; ++x){
      float d = v.p[9*idx+0]*x+v.p[9*idx+1]*y+v.p[9*idx+2];
      // Limit to disparity range [0,MAX_DISP]
      v.d[idx] = std::max(std::min(float(max_disp), d),0.f);
      idx++;
    }
}

/**
 * @brief Computes L1 norm of two 3-vectors v1 and v2
 *
 * @param v1  3-vector
 * @param v2  3-vector
 *
 * @return    l1 norm (scalar)
 */
inline int l1norm_naive(uint8_t* v1, uint8_t* v2){
  return std::abs(v1[0] - v2[0]) + std::abs(v1[1] - v2[1]) + std::abs(v1[2] - v2[2]);
}

inline bool inside(int x, int y, int lbx, int lby, int ubx, int uby) {
    return lbx <= x && x < ubx && lby <= y && y < uby;
}

/**
 * @brief Apply weighted median filter to a pixel
 *        see Definition 2.2 in https://ieeexplore.ieee.org/document/486465
 *        YIN, Lin, et al. Weighted median filters: a tutorial. IEEE Transactions on Circuits and Systems II: Analog and Digital Signal Processing, 1996, 43. Jg., Nr. 3, S. 157-192.
 *
 * @param v
 * @param cx
 * @param cy
 * @param rows
 * @param cols
 * @param ws
 */
void WeightedMedianFilter(CommonView& v, int cx, int cy, int rows, int cols, int ws){
  int half = ws / 2;
  float w_tot = 0;

  std::vector<std::pair<float, float>> disps_weights_pairs;
  disps_weights_pairs.reserve(ws * ws);

  float gamma = 10.f;
  for(int y = cy-half; y <= cy + half; ++y)
    for(int x = cx-half; x <= cx + half; ++x)
      if(inside(x, y, 0, 0, cols, rows)){
        int inorm = l1norm_naive(&v.i[(y-cy+half)*cols*3+(x-cx+half)*3], &v.i[cy*cols*3+cx*3]);
        float wt = std::exp(-inorm * GAMMA_INV);
        w_tot+=wt;
        disps_weights_pairs.push_back(std::make_pair(v.d[y*cols+x], wt));
      }

  // Sort disp_w based on their disparity value
  std::sort(disps_weights_pairs.begin(), disps_weights_pairs.end());

  // This will be the cumulative sum at which the weights
  // satisfy sum w_i_before < 0.5 and sum w_i_after < 0.5
  float med_w = w_tot / 2.0f;
  float w = 0;
  for (auto dw = disps_weights_pairs.begin(); dw < disps_weights_pairs.end(); ++dw) {
    w += dw->second;
    if (w >= med_w) {
      v.d[cy*cols+cx] = dw->first;
      break;
    }
  }
}

/**
 * @brief Post processing of both views
 *
 * @param v1
 * @param v2
 * @param rows
 * @param cols
 */
void PostProcessing(CommonView& v1, CommonView& v2, int rows, int cols){

  // checking pixels-plane disparity validity
  for(int y=0; y < rows; ++y) {
    for(int x=0; x < cols; ++x) {
      int x_rgt_match = std::round(std::max(0.f, std::min((float)cols, float(x) - v1.d[y*cols+x])));
      // 0: invalid pixel
      v1.v[y*cols+x] = (std::abs(v1.d[y*cols+x] - v2.d[y*cols+x_rgt_match])  <= 1.) ? 255 : 0;

      int x_lft_match = std::round(std::max(0.f, std::min((float)rows, float(x) + v2.d[y*cols+x])));
      v2.v[y*cols+x] = (std::abs(v2.d[y*cols+x] - v1.d[y*cols+x_lft_match]) <= 1.) ? 255 : 0;
    }
  }

  // fill-in holes related to invalid pixels
  for(int y=0; y < rows; y++) {
    for (int x=0; x < cols; x++) {
      if (v1.v[y*cols+x] == 0)
        fill_invalid_pixels(v1,x,y,rows,cols);
      if (v2.v[y*cols+x] == 0)
        fill_invalid_pixels(v2,x,y,rows,cols);
    }
  }
  // Update disparities based on new plane estimates
  PlanesToDisparity(v1, rows, cols);
  PlanesToDisparity(v2, rows, cols);

  Image disp1 = Image(rows,cols,1,8, v1.d);
  Image disp2 = Image(rows,cols,1,8, v2.d);
  //writePNG(disp1, "disp1_before_median.png");
  //writePNG(disp2, "disp2_before_median.png");

  // apply weighted median filter on hole-filled disparities
  for(int x=0; x<cols; ++x) {
    for(int y=0; y<rows; ++y) {
      if (v1.v[y*cols+x] == 0)
        WeightedMedianFilter(v1, x,y, rows, cols, WINDOW_SIZE);
      if (v2.v[y*cols+x] == 0)
        WeightedMedianFilter(v2, x,y, rows, cols, WINDOW_SIZE);
    }
  }
}

/**
 * @brief PreProcessing
 * @param img1
 * @param img2
 * @param v1
 * @param v2
 * @param rows
 * @param cols
 * @param outputDirectory
 */
void PreProcessing(const Image& img1, const Image& img2,
                   CommonView& v1, CommonView& v2,
                   int rows, int cols,
                   std::string outputDirectory) {
  allocate_common_view(v1, rows, cols);
  allocate_common_view(v2, rows, cols);
  // Convert to RGB arrays
  convert_img_to_array(img1, v1.i);
  convert_img_to_array(img2, v2.i);
  // Random planes, normals, points (note that the p1 and p2 arrays contain plane params, normals AND points)
  cout << "PM: computing random planes" << endl;
  randomInit(v1.p,rows,cols);
  randomInit(v2.p,rows,cols);
  // Compute gradients on grayscale images
  compute_gradients(v1, rows, cols);
  compute_gradients(v2, rows, cols);
  // Convert array to Image(to save to file)
  Image grad1 = Image(rows,cols,1,8, v1.g);
  writePNG(grad1, outputDirectory + "/grad1.png");
  cout << "PM: computed gradients" << endl;
}

#endif // HELPERS_HPP
