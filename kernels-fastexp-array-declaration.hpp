#include <algorithm>
#include <cstring> // for memcpy
#include <iostream>
#include <limits>
#include <random>

#include "base_kernel.hpp"
#include "helpers.hpp" // needed for init
#include "parameters.hpp"
#include "fast_exp_array_declaration.cpp"

namespace pm{
namespace fastexp_array_declaration{

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
  // Planes (format: ABCN1N2N3P1P2P3)
  // ABC: plane coeffs
  // N1N2N3 :  normal coeffs
  // P1P2P3 :  point coeffs
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
 * @param v1    View1 (left)
 * @param v2    View2 (right)
 * @param fp    current plane
 * @param x     current pixel x coord
 * @param y     current pixel y coord
 * @param rows  rows in view
 * @param cols  cols in view
 * @param cpv   indicates which one is the work view. false: left, false: right
 *
 * @return      matching cost
 */
float mcost(KernelView& v1, KernelView& v2, float* fp, int x, int y, int rows, int cols, bool cpv){
#ifdef TIME_MCOST
  myInt64 start = start_tsc();
#endif
  int sign = (cpv == false) ? -1 : 1;	// -1 processing left, +1 processing right
  // Work view
  KernelView& wv = (cpv == false) ? v1 : v2;
  // The "other view"
  KernelView& ov = (cpv == false) ? v2 : v1;

  int HALF_WIN = WINDOW_SIZE/2;
  float cost = 0.f;
  for(int qy = y - HALF_WIN; qy <= y + HALF_WIN; ++qy){
    for(int qx = x - HALF_WIN; qx <= x + HALF_WIN; ++qx)
      if(inside(qx, qy, 0, 0, cols, rows)){

        //Coords of q'(in the other frame)
        //y coordinate will be the same, only x is affected by disparity.
        float disp =  (fp[0]*float(qx) + fp[1]*float(qy) + fp[2]);
        if(disp < 0 || disp > max_disp){
          return std::numeric_limits<float>::infinity();
        }else{
          float match = qx + sign * disp;

          if (match > cols - 2) {
            match = cols - 2;
          }
          if (match < 0) {
            match = 0;
          }

          int qdx = (int)match;
          int qdy = qy;
          float wm = 1.f - (match - qdx);

          // Weight between p and q
          int inorm = l1norm_naive(&wv.i[y*cols*3+x*3], &wv.i[qy*cols*3+qx*3]);
          float w = fast_exp(-inorm);

          // intensity difference between this and other view
          float iqnorm = weighted_l1norm(&wv.i[qy*cols*3+qx*3], &ov.i[qdy*cols*3+qdx*3], &ov.i[qdy*cols*3+(qdx+1)*3], wm);

          // gradient intensity difference between this and other view
          float ovg  = wm*ov.g[qdy*cols+qdx] + (1.f-wm)*ov.g[qdy*cols+qdx+1];
          float iqgnorm = std::abs(wv.g[qy*cols+qx] - ovg);

          float rho=(1.f-ALPHA)*std::min(iqnorm,TAUCOL) + ALPHA*std::min(iqgnorm,TAUGRAD);

          cost += w*rho;
        }
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
void evaluatePlanesCost(KernelView& v1, KernelView& v2, int rows, int cols, bool cpv){
  KernelView& wv = (cpv == false) ? v1 : v2;
  KernelView& ov = (cpv == false) ? v2 : v1;
  for(int y=0; y<rows; ++y)
    for(int x=0; x<cols; ++x){
        float* cc = &(wv.c[(y*cols)+x]);
        float* fp = &(wv.p[(y*cols*9)+x*9]);
				*cc = mcost(v1,v2,fp,x,y,rows,cols,cpv);
    }
}
/**
 * @brief Spatial propagation
 *
 * @param v1      Left view
 * @param v2      Right view
 * @param x       current x-coord
 * @param y       current y coord
 * @param rows    number of rows in view
 * @param cols    number of cols in view
 * @param cpv     indicates which one is the work view. false: left, false: right
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void SpatialPropagation(KernelView& v1, KernelView& v2, int x, int y, int rows, int cols, bool cpv, bool isEven){

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
  KernelView& wv = (cpv == false) ? v1 : v2;
  KernelView& ov = (cpv == false) ? v2 : v1;

  //old plane, old cost
  float* plane_old = &(wv.p[(y*cols*9)+x*9]);
  float* cost_old = &(wv.c[y*cols+x]);
  if(n1in){
    //neighbor planes
    float* plane_new = &(wv.p[(n1y*cols*9)+n1x*9]);
    float cost_new =  mcost(v1,v2,plane_new,x,y,rows,cols,cpv);

    if(cost_new < *cost_old){
      memcpy(plane_old,plane_new,9*sizeof(float));
      *cost_old = cost_new;
    }
  }
  if(n2in){
    //neighbor planes
    float* plane_new = &(wv.p[(n2y*cols*9)+n2x*9]);
    float cost_new =  mcost(v1,v2,plane_new,x,y,rows,cols,cpv);

    if(cost_new < *cost_old){
      memcpy(plane_old,plane_new,9*sizeof(float));
      *cost_old = cost_new;
    }
  }
}

/**
 * @brief View propagation step
 *
 * @param v1      Left view
 * @param v2      Right view
 * @param x       current x-coord
 * @param y       current y coord
 * @param rows    number of rows in view
 * @param cols    number of cols in view
 * @param cpv     indicates which one is the work view. false: left, false: right
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void ViewPropagation(KernelView& v1, KernelView& v2, int x, int y, int rows, int cols, bool cpv, bool isEven){
  int sign = (cpv == false) ? -1 : 1;	// -1 processing left, +1 processing right
  KernelView& wv = (cpv == false) ? v1 : v2;
  KernelView& ov = (cpv == false) ? v2 : v1;
  // current plane
  float* fp = &wv.p[(y * cols * 9) + x * 9];
  // TODO: precompute candidates instead of iterating over whole epipolar line
  // check epipolar line in other view
  for (int x_other = 0; x_other < cols; ++x_other) {
      float* fpother = &ov.p[(y * cols * 9) + x_other * 9];
      float z = fpother[0] * x_other + fpother[1] * y + fpother[2];
      // compute matching point in work view, note the minus
      int mx = roundf(x_other - sign * z);
      if (mx != x) {    // no match
          continue;
      }
      else {
        int my = y;
        // Copy over same normal. thus a,b will be the same, c will change:
        float c = (fpother[3] * mx + fpother[4] * my + fpother[5] * z) / fpother[5];
        float new_plane[9] = {fp[0], fp[1], c, fp[3], fp[4], fp[5], (float)mx, (float)my, z};

        float* cost_old = &(wv.c[y * cols + x]);
        float  cost_new = mcost(v1, v2, new_plane, x, y, rows, cols, cpv);
        if(cost_new < *cost_old){
          //Update the plane
          memcpy(fp, &new_plane, 9 * sizeof(float));
          *cost_old = cost_new;
        }
      }
  }
}

/**
 * @brief Plane refinement step
 *
 * @param v1      Left view
 * @param v2      Right view
 * @param x       current x-coord
 * @param y       current y coord
 * @param rows    number of rows in view
 * @param cols    number of cols in view
 * @param cpv     indicates which one is the work view. false: left, false: right
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void PlaneRefinement(KernelView& v1, KernelView& v2, int x, int y, int rows, int cols, bool cpv, bool isEven){
  int sign = (cpv==0) ? -1 : 1;	// -1 processing left, +1 processing right
  KernelView& wv = (cpv == false) ? v1 : v2;
  KernelView& ov = (cpv == false) ? v2 : v1;

  float max_dz  = max_disp / 2.f;
  float max_dn = 1.0f;
  float end_dz = 0.1f;

  //Current pixel's plane and matching cost
  float* plane_old = &wv.p[(y*cols*9)+x*9];
  float* cost_old = &(wv.c[y*cols+x]);

  //Buffer for new plane proposal
  float plane[9];
  // Searching a random plane starting from the actual one
  while(max_dz >= end_dz)
  {
    std::uniform_real_distribution<float> dz_dis(-max_dz, +max_dz);
    std::uniform_real_distribution<float> dn_dis(-max_dn, +max_dn);

    float zold = plane_old[0] * x + plane_old[1] * y + plane_old[2];

    // New point
    plane[6] = x;
    plane[7] = y;
    plane[8] = zold + dz_dis(gen); //delta_z

    // New normal
    float nx = plane_old[3] + dn_dis(gen);
    float ny = plane_old[4] + dn_dis(gen);
    float nz = plane_old[5] + dn_dis(gen);

    nz = nz == 0.f ? 1e-18f : nz;

    //Normalize new normal
    float n = sqrt(nx*nx + ny*ny + nz * nz);
    plane[3] = nx/n;
    plane[4] = ny/n;
    plane[5] = nz/n;

    // Plane params
    plane[0] = -plane[3]/plane[5];
    plane[1] = -plane[4]/plane[5];
    plane[2] = (plane[3]*plane[6]+plane[4]*plane[7]+plane[5]*plane[8])/plane[5];

    // test the new plane
    // old_cost can be moved out of loop, only need it first time
    float cost_new = mcost(v1,v2,plane,x,y,rows,cols,cpv);

    if(cost_new < *cost_old){
      memcpy(plane_old,&plane,9*sizeof(float));
      *cost_old = cost_new;
    }

    max_dz /= 2.0f;
    max_dn /= 2.0f;
  }
}

/**
 * @brief Processes a single pixel
 *
 * @param v1      Left view
 * @param v2      Right view
 * @param x       x-coord
 * @param y       y-coord
 * @param rows    rows in view
 * @param cols    cols in view
 * @param cpv     indicates which one is the work view. false: left, false: right
 * @param isEven  indicates if this is an even iteration (decides which neighbors we look at.
 */
void processPixel(KernelView& v1, KernelView& v2, int x, int y, int rows, int cols, bool cpv, bool isEven){
  SpatialPropagation(v1,v2, x,y,rows,cols,cpv,isEven);
  ViewPropagation(v1,v2, x,y,rows,cols,cpv,isEven);
  PlaneRefinement(v1,v2, x,y,rows,cols,cpv,isEven);
}

void process(KernelView& v1, KernelView& v2, int rows, int cols) {
  // Eval plane's cost
  // TODO: should this be moved to a helper?
  evaluatePlanesCost(v1, v2, rows, cols, 0);
  evaluatePlanesCost(v1, v2, rows, cols, 1);
  std::cout << "PM: evaluated plane cost" << std::endl;

  for(int it = 0; it < 3; it++){
    std::cout << "Iteration " << it << std::endl;
    bool isOdd = it&1;
    bool isEven = !isOdd;

    for(int work_view=0; work_view < 2; ++work_view){
      if(isEven){
        // Top down
        for(int y=0;y<rows;y++){
          if(( y % 50 ) == 0) std::cout << "y:" << y << " / " <<  rows << std::endl;
          for(int x=0;x<cols;x++){
            processPixel(v1,v2, x,y,rows,cols,work_view,isEven);
          }
        }
      }else{
        // Bottom up
        for(int y=rows-1; y>=0;--y){
          if(( y % 50 ) == 0) std::cout << "y:" << y << " / " <<  rows << std::endl;
          for(int x=cols-1;x>=0;--x){
            processPixel(v1,v2, x,y,rows,cols,work_view,isEven);
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
    // use data of CommonView (no copying)
    // Images
    v1_.i = v1.i;
    v2_.i = v2.i;
    // Gradients
    v1_.g = v1.g;
    v2_.g = v2.g;
    // Planes
    v1_.p = v1.p;
    v2_.p = v2.p;
    // Costs
    v1_.c = v1.c;
    v2_.c = v2.c;
  }

  void run_patch_match() {
    process(v1_, v2_, rows_, cols_);
  }

  void update_common_view(CommonView& v1, CommonView& v2) const {
    // nothing since we work directly on the Common View data
  }

  std::pair<float, float> test_mcost(int x, int y) {
    float* fp = &(v1_.p[(y*cols_*9)+x*9]);
    float cost_left = mcost(v1_, v2_, fp, x, y, rows_, cols_, false);
    fp = &(v2_.p[(y*cols_*9)+x*9]);
    float cost_right = mcost(v1_, v2_, fp, x, y, rows_, cols_, true);
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

    return mcost(v1_, v2_, fp, x, y, rows_, cols_, true);
  }
  // destructor
  ~Kernel() {
    // nothing since we work directly on the Common View data
  }

private:
  KernelView v1_;
  KernelView v2_;
};

}//namespace fastexp_array
}//namespace pm
