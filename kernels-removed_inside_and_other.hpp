#include <algorithm>
#include <cstring> // for memcpy
#include <iostream>
#include <limits>
#include <random>

#include "base_kernel.hpp"
#include "helpers.hpp" // needed for init
#include "parameters.hpp"

/**
 * Based on compressed_planes
 */

namespace pm{
namespace removed_inside_and_other{

// Global random number generator, fixed random seed
std::mt19937 gen(42);

/**
 * @brief The View struct
 *
 * Implementations specif struct, same as common View struct in baseline implementation
 */
struct KernelView{
  // All data in this struct is stored in row order
  // Image in HWC format i.e. row order and channels interleaved (R1G1B1,R2G2B2, ...)
  uint8_t* i;
  // Gradient
  float* g;
  // Planes (format: ABC)
  // ABC: plane coeffs
  float* p;
  // cost
  float* c;
};

/**
 * @brief Checks if a pixel x,y lies within
 *        the bounding rectangle spanned by lbx, lby,ubx,uby.
 *
 * @param x
 * @param y
 * @param lbx
 * @param lby
 * @param ubx
 * @param uby
 *
 * @return true if inside
 */
inline bool inside(int x, int y, int lbx, int lby, int ubx, int uby) {
		return lbx <= x && x < ubx && lby <= y && y < uby;
}

/**
 * @brief Computes the l1 norm of wm*v2a+(1-wm)*v2b-v1
 *
 * @param v1    3-vector
 * @param v2a   3-vector component a
 * @param v2b   3-vector component b
 * @param wm    weighting factor to mix components a and b
 *
 * @return  ||wm*v2a+(1-wm)*v2b-v1||_1
 */
inline float weighted_l1norm(uint8_t* v1, uint8_t* v2a, uint8_t* v2b, float wm){
  float vp1 =  wm * v2a[0] + (1.f - wm) * v2b[0] - v1[0];
  float vp2 =  wm * v2a[1] + (1.f - wm) * v2b[1] - v1[1];
  float vp3 =  wm * v2a[2] + (1.f - wm) * v2b[2] - v1[2];
  return std::abs(vp1) + std::abs(vp2) + std::abs(vp3);
}

/**
 * @brief Computes the cost function m.
 *
 * @param wv      Working view
 * @param ov      Other view
 * @param fp    current plane
 * @param x     current pixel x coord
 * @param y     current pixel y coord
 * @param rows  rows in view
 * @param cols  cols in view
 * @param cpv   indicates which one is the work view. false: left, false: right
 *
 * @return      matching cost
 */
float mcost(KernelView& wv, KernelView& ov, float* fp, int x, int y, int rows, int cols, int sign,
            int qy_start, int qy_end, int qx_start, int qx_end){
#ifdef TIME_MCOST
  myInt64 start = start_tsc();
#endif
  float cost = 0.f;
  for(int qy = qy_start; qy <= qy_end; ++qy){
    float disp = fp[0] * qx_start + fp[1] * qy + fp[2];
    for(int qx = qx_start; qx <= qx_end; ++qx){
        if(disp < 0 || disp > max_disp){
          return std::numeric_limits<float>::infinity();
        }else{
          float match = qx + sign * disp;
          match = match > cols - 2 ? cols - 2 : match < 0 ? 0 : match;

          int qdx = (int)match;
          int qdy = qy;
          float wm = 1.f - (match - qdx);

          // Weight between p and q
          int inorm = l1norm_naive(&wv.i[y*cols*3+x*3], &wv.i[qy*cols*3+qx*3]);
          float w = std::exp(-inorm/GAMMA);

          // intensity difference between this and other view
          float iqnorm = weighted_l1norm(&wv.i[qy*cols*3+qx*3], &ov.i[qdy*cols*3+qdx*3], &ov.i[qdy*cols*3+(qdx+1)*3], wm);

          // gradient intensity difference between this and other view
          float ovg  = wm*ov.g[qdy*cols+qdx] + (1.f-wm)*ov.g[qdy*cols+qdx+1];
          float iqgnorm = std::abs(wv.g[qy*cols+qx] - ovg);

          float rho=(1.f-ALPHA)*std::min(iqnorm,TAUCOL) + ALPHA*std::min(iqgnorm,TAUGRAD);

          cost += w*rho;
        }
        disp += fp[0];
    }
  }
#ifdef TIME_MCOST
  myInt64 end = stop_tsc(start);
  pm::mcost_total_time_ += end;
  pm::mcost_calls_ += 1;
#endif
  
  return cost;
}

/**
 * @brief Initial  cost evaluation for all planes
 *
 * @param v1    View1 (left)
 * @param v2    View2 (right)
 * @param rows  number of rows in view
 * @param cols  number of cols in view
 * @param cpv   indicates which one is the work view. false: left, false: right
 */
void evaluatePlanesCost(KernelView& wv, KernelView& ov, int rows, int cols, int sign){
  for(int y=0; y<rows; ++y)
    for(int x=0; x<cols; ++x){
      int HALF_WIN = WINDOW_SIZE/2;
      int qy_start = y - HALF_WIN >= 0 ? y - HALF_WIN : 0;
      int qy_end = y + HALF_WIN < rows ? y + HALF_WIN : rows - 1;
      int qx_start = x - HALF_WIN >= 0 ? x - HALF_WIN : 0;
      int qx_end = x + HALF_WIN < cols ? x + HALF_WIN : cols - 1;
      float* cc = &(wv.c[(y*cols)+x]);
      float* fp = &(wv.p[(y*cols*3)+x*3]);
      *cc = mcost(wv,ov,fp,x,y,rows,cols,sign, qy_start, qy_end, qx_start, qx_end);
    }
}
/**
 * @brief Spatial propagation
 *
 * @param wv      Working view
 * @param ov      Other view
 * @param x       current x-coord
 * @param y       current y coord
 * @param rows    number of rows in view
 * @param cols    number of cols in view
 * @param sign    sign for adding/subtracting of disparity
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void SpatialPropagation(KernelView& wv, KernelView& ov, int x, int y, int rows, int cols, int sign,
                        bool isEven, int qy_start, int qy_end, int qx_start, int qx_end){
  int n1x,n1y,n2x,n2y;
  bool n1in,n2in;
  if(isEven){ //odd iteration: right and lower neighbor
    n1x = x-1;
    n1y = y;
    n2x = x;
    n2y = y-1;
  }else{
    n1x = x+1;
    n1y = y;
    n2x = x;
    n2y = y+1;
  }
  n1in = inside(n1x, n1y, 0, 0, cols, rows);
  n2in = inside(n2x, n2y, 0, 0, cols, rows);

  //old plane, old cost
  float* plane_old = &(wv.p[(y*cols*3)+x*3]);
  float* cost_old = &(wv.c[y*cols+x]);
  if(n1in){
    //neighbor planes
    float* plane_new = &(wv.p[(n1y*cols*3)+n1x*3]);
    float cost_new =  mcost(wv,ov,plane_new,x,y,rows,cols,sign, qy_start, qy_end, qx_start, qx_end);

    if(cost_new < *cost_old){
      memcpy(plane_old,plane_new,3*sizeof(float));
      *cost_old = cost_new;
    }
  }
  if(n2in){
    //neighbor planes
    float* plane_new = &(wv.p[(n2y*cols*3)+n2x*3]);
    float cost_new =  mcost(wv,ov,plane_new,x,y,rows,cols,sign, qy_start, qy_end, qx_start, qx_end);

    if(cost_new < *cost_old){
      memcpy(plane_old,plane_new,3*sizeof(float));
      *cost_old = cost_new;
    }
  }
}

/**
 * @brief View propagation step
 *
 * @param wv      Working view
 * @param ov      Other view
 * @param x       current x-coord
 * @param y       current y coord
 * @param rows    number of rows in view
 * @param cols    number of cols in view
 * @param sign    sign for adding/subtracting of disparity
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void ViewPropagation(KernelView& wv, KernelView& ov, int x, int y, int rows, int cols, int sign,
                     bool isEven, int qy_start, int qy_end, int qx_start, int qx_end){
  // current plane
  float* fp = &wv.p[(y * cols * 3) + x * 3];
  // TODO: precompute candidates instead of iterating over whole epipolar line
  // check epipolar line in other view
  for (int x_other = 0; x_other < cols; ++x_other) {
      float* fpother = &ov.p[(y * cols * 3) + x_other * 3];
      float z = fpother[0] * x_other + fpother[1] * y + fpother[2];
      // compute matching point in work view, note the minus
      int mx = roundf(x_other - sign * z);
      if (mx != x) {    // no match
          continue;
      }
      else {
        int my = y;
        // Copy over same normal. thus a,b will be the same, c will change:
        float c = fpother[0] * mx + fpother[1] * my + z;
        float new_plane[3] = {fp[0], fp[1], c};

        float* cost_old = &(wv.c[y * cols + x]);
        float  cost_new = mcost(wv, ov, new_plane, x, y, rows, cols, sign, qy_start, qy_end, qx_start, qx_end);
        if(cost_new < *cost_old){
          //Update the plane
          memcpy(fp, &new_plane, 3 * sizeof(float));
          *cost_old = cost_new;
        }
      }
  }
}

/**
 * @brief Plane refinement step
 *
 * @param wv      Working view
 * @param ov      Other view
 * @param x       current x-coord
 * @param y       current y coord
 * @param rows    number of rows in view
 * @param cols    number of cols in view
 * @param sign    sign for adding/subtracting of disparity
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void PlaneRefinement(KernelView& wv, KernelView& ov, int x, int y, int rows, int cols, int sign,
                     bool isEven, int qy_start, int qy_end, int qx_start, int qx_end){
  float max_dz  = max_disp / 2.f;
  float max_dn = 1.0f;
  float end_dz = 0.1f;

  //Current pixel's plane and matching cost
  float* plane_old = &wv.p[(y*cols*3)+x*3];
  float* cost_old = &(wv.c[y*cols+x]);

  float z_old = plane_old[0] * x + plane_old[1] * y + plane_old[2];

  // get normal: (-a, -b, 1).normalize()
  float norm_inv = 1.f / std::sqrt(plane_old[0]*plane_old[0] + plane_old[1]*plane_old[1] + 1.f);
  float nx_old = -plane_old[0] * norm_inv;
  float ny_old = -plane_old[1] * norm_inv;
  float nz_old = norm_inv;

  //Buffer for new plane proposal
  float plane[3];
  // Searching a random plane starting from the actual one
  while(max_dz >= end_dz)
  {
    std::uniform_real_distribution<float> dz_dis(-max_dz, +max_dz);
    std::uniform_real_distribution<float> dn_dis(-max_dn, +max_dn);

    // New point
    float z = z_old + dz_dis(gen); //delta_z

    // New normal
    float nx = nx_old + dn_dis(gen);
    float ny = ny_old + dn_dis(gen);
    float nz = nz_old + dn_dis(gen);
    nz = nz == 0.f ? 1e-18f : nz;

    //Normalize new normal
    float n = sqrt(nx * nx + ny * ny + nz * nz);
    nx = nx / n;
    ny = ny / n;
    nz = nz / n;

    // Plane params
    plane[0] = -nx / nz;
    plane[1] = -ny / nz;
    plane[2] = (nx * x + ny * y + nz * z) / nz;

    // test the new plane
    // old_cost can be moved out of loop, only need it first time
    float cost_new = mcost(wv,ov,plane,x,y,rows,cols,sign, qy_start, qy_end, qx_start, qx_end);

    if(cost_new < *cost_old){
      memcpy(plane_old,&plane,3*sizeof(float));
      *cost_old = cost_new;
      z_old = z;
      nx_old = nx;
      ny_old = ny;
      nz_old = nz;
    }

    max_dz /= 2.0f;
    max_dn /= 2.0f;
  }
}

/**
 * @brief Processes a single pixel
 *
 * @param wv      Working view
 * @param ov      Other view
 * @param x       x-coord
 * @param y       y-coord
 * @param rows    rows in view
 * @param cols    cols in view
 * @param sign    sign for adding/subtracting of disparity
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void processPixel(KernelView& wv, KernelView& ov, int x, int y, int rows, int cols, int sign, bool isEven){
  int HALF_WIN = WINDOW_SIZE/2;
  int qy_start = y - HALF_WIN >= 0 ? y - HALF_WIN : 0;
  int qy_end = y + HALF_WIN < rows ? y + HALF_WIN : rows - 1;
  int qx_start = x - HALF_WIN >= 0 ? x - HALF_WIN : 0;
  int qx_end = x + HALF_WIN < cols ? x + HALF_WIN : cols - 1;
  SpatialPropagation(wv,ov, x,y,rows,cols,sign,isEven, qy_start, qy_end, qx_start, qx_end);
  ViewPropagation(wv,ov, x,y,rows,cols,sign,isEven, qy_start, qy_end, qx_start, qx_end);
  PlaneRefinement(wv,ov, x,y,rows,cols,sign,isEven, qy_start, qy_end, qx_start, qx_end);
}

void process(KernelView& v1, KernelView& v2, int rows, int cols) {
  // Eval plane's cost
  // TODO: should this be moved to a helper?
  evaluatePlanesCost(v1, v2, rows, cols, -1);
  evaluatePlanesCost(v2, v1, rows, cols, 1);
  std::cout << "PM: evaluated plane cost" << std::endl;

  for(int it = 0; it < 3; it++){
    std::cout << "Iteration " << it << std::endl;
    bool isOdd = it&1;
    bool isEven = !isOdd;

    for(int work_view=0; work_view < 2; ++work_view){
      int sign = (work_view == false) ? -1 : 1;	// -1 processing left, +1 processing right
      // Work view
      KernelView& wv = (work_view == false) ? v1 : v2;
      // The "other view"
      KernelView& ov = (work_view == false) ? v2 : v1;
      if(isEven){
        // Top down
        for(int y=0;y<rows;y++){
          if(( y % 50 ) == 0) std::cout << "y:" << y << " / " <<  rows << std::endl;
          for(int x=0;x<cols;x++){
            processPixel(wv,ov, x,y,rows,cols,sign,isEven);
          }
        }
      }else{
        // Bottom up
        for(int y=rows-1; y>=0;--y){
          if(( y % 50 ) == 0) std::cout << "y:" << y << " / " <<  rows << std::endl;
          for(int x=cols-1;x>=0;--x){
            processPixel(wv,ov, x,y,rows,cols,sign,isEven);
          }
        }
      }
    }
  }
}

class Kernel : public BaseKernel{
public:
  // init
  Kernel(const CommonView& v1, const CommonView& v2, int rows, int cols) : BaseKernel(rows, cols) {
    // use some data of CommonView
    // Images
    v1_.i = v1.i;
    v2_.i = v2.i;
    // Gradients
    v1_.g = v1.g;
    v2_.g = v2.g;
    // Planes
    v1_.p = (float*)malloc(rows*cols*3*sizeof(float));
    v2_.p = (float*)malloc(rows*cols*3*sizeof(float));
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        v1_.p[y * cols * 3 + x * 3 + 0] = v1.p[y * cols * 9 + x * 9 + 0];
        v1_.p[y * cols * 3 + x * 3 + 1] = v1.p[y * cols * 9 + x * 9 + 1];
        v1_.p[y * cols * 3 + x * 3 + 2] = v1.p[y * cols * 9 + x * 9 + 2];
      }
    }
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        v2_.p[y * cols * 3 + x * 3 + 0] = v2.p[y * cols * 9 + x * 9 + 0];
        v2_.p[y * cols * 3 + x * 3 + 1] = v2.p[y * cols * 9 + x * 9 + 1];
        v2_.p[y * cols * 3 + x * 3 + 2] = v2.p[y * cols * 9 + x * 9 + 2];
      }
    }
    // Costs
    v1_.c = v1.c;
    v2_.c = v2.c;
  }

  void run_patch_match() {
    process(v1_, v2_, rows_, cols_);
  }

  void update_common_view(CommonView& v1, CommonView& v2) const {
    for (int y = 0; y < rows_; ++y) {
      for (int x = 0; x < cols_; ++x) {
        v1.p[y * cols_ * 9 + x * 9 + 0] = v1_.p[y * cols_ * 3 + x * 3 + 0];
        v1.p[y * cols_ * 9 + x * 9 + 1] = v1_.p[y * cols_ * 3 + x * 3 + 1];
        v1.p[y * cols_ * 9 + x * 9 + 2] = v1_.p[y * cols_ * 3 + x * 3 + 2];

      }
    }
    for (int y = 0; y < rows_; ++y) {
      for (int x = 0; x < cols_; ++x) {
        v2.p[y * cols_ * 9 + x * 9 + 0] = v2_.p[y * cols_ * 3 + x * 3 + 0];
        v2.p[y * cols_ * 9 + x * 9 + 1] = v2_.p[y * cols_ * 3 + x * 3 + 1];
        v2.p[y * cols_ * 9 + x * 9 + 2] = v2_.p[y * cols_ * 3 + x * 3 + 2];
      }
    }
  }

  std::pair<float, float> test_mcost(int x, int y) {
    int HALF_WIN = WINDOW_SIZE/2;
    int qy_start = y - HALF_WIN >= 0 ? y - HALF_WIN : 0;
    int qy_end = y + HALF_WIN < rows_ ? y + HALF_WIN : rows_ - 1;
    int qx_start = x - HALF_WIN >= 0 ? x - HALF_WIN : 0;
    int qx_end = x + HALF_WIN < cols_ ? x + HALF_WIN : cols_ - 1;
    float* fp = &(v1_.p[(y*cols_*3)+x*3]);
    float cost_left = mcost(v1_, v2_, fp, x, y, rows_, cols_, -1, qy_start, qy_end, qx_start, qx_end);
    fp = &(v2_.p[(y*cols_*3)+x*3]);
    float cost_right = mcost(v2_, v1_, fp, x, y, rows_, cols_, +1, qy_start, qy_end, qx_start, qx_end);
    return {cost_left, cost_right};
  }
  float get_W_mcost(void){
    return 1.0 * 48. * WINDOW_SIZE * WINDOW_SIZE;
  }
  float get_Q_mcost(void){
    return 1.0 * 29. * WINDOW_SIZE * WINDOW_SIZE;
  }
  float peakperf_mcost(void){
    // The following two measures ensure that all pixels 
    // of the cost window are evaluated (upper bound)
    // 
    // Set center coordinate of cost window s.t.
    // complete cost window is within image
    int x = 50;
    int y = 50;
    // Yields disparity = 0 in any case -> match within image
    float fp[] = {0,0,0,0}; 
    int HALF_WIN = WINDOW_SIZE/2;
    int qy_start = y - HALF_WIN >= 0 ? y - HALF_WIN : 0;
    int qy_end = y + HALF_WIN < rows_ ? y + HALF_WIN : rows_ - 1;
    int qx_start = x - HALF_WIN >= 0 ? x - HALF_WIN : 0;
    int qx_end = x + HALF_WIN < cols_ ? x + HALF_WIN : cols_ - 1;

    return mcost(v2_, v1_, fp, x, y, rows_, cols_, +1, qy_start, qy_end, qx_start, qx_end);
  }

  // destructor
  ~Kernel() {
    // nothing since we work directly on the Common View data
    delete[] v1_.p;
    delete[] v2_.p;
  }

private:
  KernelView v1_;
  KernelView v2_;
};

}//namespace removed_inside_and_other
}//namespace pm
