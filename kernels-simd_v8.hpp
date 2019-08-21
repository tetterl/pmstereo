#include <algorithm>
#include <cstring> // for memcpy
#include <iostream>
#include <limits>
#include <random>

#include <immintrin.h>

#include "base_kernel.hpp"
#include "helpers.hpp" // needed for init
#include "parameters.hpp"

/**
 * Based on simd_v7
 * merged gather -> slower
 */

namespace pm{
namespace simd_v8{

float precomputed_vals[771];
void setup_exp() {
  precomputed_vals[0] = std::exp(0);
  for (int i = 1; i < 771; ++i) {
    precomputed_vals[i] = std::exp(-i * GAMMA_INV);
  }
}
inline float fast_exp(int val) {
  // Assume value is in [0,n]
  return precomputed_vals[val];
}


// Global random number generator, fixed random seed
std::mt19937 gen(42);

// global array for precomputed weights
// note that the weights are stored sequentially according to the access order
// in boundary regions the weight matrix is not fully filled with valid weights
// the valid weights don't necessarily form a rectangle
float weights[WINDOW_SIZE * WINDOW_SIZE + 8] __attribute__ ((aligned (16))) ;

/**
 * @brief The View struct
 *
 * Implementations specif struct, same as common View struct in baseline implementation
 */
struct KernelView{
  // All data in this struct is stored in row order
  // row order and channels interleaved (R1G1B1A1,R2G2B2A1, ...)
  // the alpha channel is always equal to zero
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

  // check if disparities out of range
  float disp_11 = fp[0] * qx_start + fp[1] * qy_start + fp[2];
  float disp_21 = fp[0] * qx_end + fp[1] * qy_start + fp[2];
  float disp_12 = fp[0] * qx_start + fp[1] * qy_end + fp[2];
  float disp_22 = fp[0] * qx_end + fp[1] * qy_end + fp[2];
  if(
     disp_11 < 0 || disp_11 > max_disp ||
     disp_21 < 0 || disp_21 > max_disp ||
     disp_12 < 0 || disp_12 > max_disp ||
     disp_22 < 0 || disp_22 > max_disp
     ) {
    return std::numeric_limits<float>::infinity();
  }

  // TODO: store globally or similar?!
  __m256 ones = _mm256_set1_ps(1);
  __m256i onesi = _mm256_set1_epi32(1);
  // from pixel x to pixel x+1 we have a disparity increase of +/-fp[0] and for the shift +1
  __m256 match_increase = _mm256_set1_ps(8*(1+sign*fp[0]));
  __m256i eightsi = _mm256_set1_epi32(8);
  __m256i trues = _mm256_cmpeq_epi32(eightsi, eightsi);
  __m256 fp0 = _mm256_set1_ps(fp[0]);
  __m256 fp1 = _mm256_set1_ps(fp[1]);
  __m256 fp2 = _mm256_set1_ps(fp[2]);
  __m256 zerosf = _mm256_setzero_ps();
  __m256i mcols = _mm256_set1_epi32(cols);
  __m256 colsm2 = _mm256_set1_ps(cols - 2);
  __m256 minuszeros = _mm256_set1_ps(-0.f);
  __m256 signs = _mm256_set1_ps(sign);
  __m256 taucol = _mm256_set1_ps(TAUCOL);
  __m256 taugrad = _mm256_set1_ps(TAUGRAD);
  __m256i qx_ends = _mm256_set1_epi32(qx_end);
  __m256i offset = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
  __m256i perm0 = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
  __m256i perm1 = _mm256_set_epi32(3, 3, 3, 3, 2, 2, 2, 2);
  __m256i perm2 = _mm256_set_epi32(5, 5, 5, 5, 4, 4, 4, 4);
  __m256i perm3 = _mm256_set_epi32(7, 7, 7, 7, 6, 6, 6, 6);
  __m256i perm_sum = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);

  // accumulators for cost
  // float cost_img_diff = 0.f;
  // float cost_grad_diff = 0.f;
  __m256 cost_img_diff = _mm256_setzero_ps();
  __m256 cost_grad_diff = _mm256_setzero_ps();
  float disp_tmp = fp[1] * qy_start + fp[2];
  int idx = 0;
  __m256i qdycols = _mm256_set1_epi32(qy_start*cols);
  __m256i qx_starts = _mm256_set1_epi32(qx_start);
  __m256i qxis_init = _mm256_add_epi32(qx_starts, offset);
  __m256 qxs_init = _mm256_cvtepi32_ps(qxis_init);
  for(int qy = qy_start; qy <= qy_end; ++qy){
    __m256i qxis = qxis_init;
    // float disp = fp[0] * qx_start + disp_tmp;
    __m256 qys = _mm256_set1_ps(qy);
    __m256 disp = _mm256_fmadd_ps(fp1, qys, _mm256_fmadd_ps(fp0, qxs_init, fp2));  // reordered to hide latency
    __m256 match_unclamped = _mm256_fmadd_ps(signs, disp, qxs_init);
    idx = (qy - qy_start) * (qx_end - qx_start + 1);
    for(int qx = qx_start; qx <= qx_end; qx += 8){
      // construct mask since we might access elements not in the window, thus clean-up code can be prevented
      __m256 mask = reinterpret_cast<const __m256>(_mm256_xor_si256(trues, _mm256_cmpgt_epi32(qxis, qx_ends))); // mask = qxs <= qx_ends

      __m256 match = _mm256_max_ps(zerosf, _mm256_min_ps(colsm2, match_unclamped));  // match = match > cols - 2 ? cols - 2 : match < 0 ? 0 : match;

      __m256i qdx = _mm256_cvtps_epi32(_mm256_floor_ps(match)); // int qdx = (int)match;
      __m256 inv_fac = _mm256_sub_ps(match, _mm256_cvtepi32_ps(qdx)); // float inv_fac = 1.f - wm = match - qdx;
      __m256 wm = _mm256_sub_ps(ones, inv_fac); // float wm = 1.f - (match - qdx) = 1.f - inv_fac;

      __m256 w = _mm256_loadu_ps(weights + idx);  // float w = weights[idx];
      w = _mm256_and_ps(w, mask);

      /// Start:
      /// l1 norm computation (maximum 3/4 usage of simd lanes)
      ///
      // pixel 0-3
      __m256i idx_left = _mm256_add_epi32(qdycols, qdx);  // omit multiplication by four since we have scale
      __m256i idx_right = _mm256_add_epi32(idx_left, onesi);  // omit multiplication by four since we have scale
      // [R0,G0,B0,A0, R1,G1,B1,A1, R2,G2,B2,A2, R3,G3,B3,A3]
      // left: left matching point for weighting
      // scale = 4 since RGBA uses 4 bytes
      __m256i RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_left = _mm256_i32gather_epi32(reinterpret_cast<const int*>(ov.i), idx_left, 4);
      __m128i RGBA0_RGBA1_RGBA2_RGB3_left = _mm256_castsi256_si128(RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_left);
      __m256i RGBA0_RGBA1_left = _mm256_cvtepu8_epi32(RGBA0_RGBA1_RGBA2_RGB3_left); // convert uint8 to int32 (lower 128 lane)
      __m128i RGBA2_RGB3_RGBA0_RGBA1_left = _mm_shuffle_epi32(RGBA0_RGBA1_RGBA2_RGB3_left, _MM_SHUFFLE(1, 0, 3, 2));
      __m256i RGBA2_RGBA3_left = _mm256_cvtepu8_epi32(RGBA2_RGB3_RGBA0_RGBA1_left); // convert uint8 to int32 (lower 128 lane)
      __m256 RGBA0_RGBA1_left_f = _mm256_cvtepi32_ps(RGBA0_RGBA1_left); // convert int32 to float
      __m256 RGBA2_RGBA3_left_f = _mm256_cvtepi32_ps(RGBA2_RGBA3_left); // convert int32 to float

      // right: right matching point for weighting
      // scale = 4 since RGBA uses 4 bytes
      __m256i RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_right = _mm256_i32gather_epi32(reinterpret_cast<const int*>(ov.i), idx_right, 4);
      __m128i RGBA0_RGBA1_RGBA2_RGB3_right = _mm256_castsi256_si128(RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_right);
      __m256i RGBA0_RGBA1_right = _mm256_cvtepu8_epi32(RGBA0_RGBA1_RGBA2_RGB3_right); // convert uint8 to int32 (lower 128 lane)
      __m128i RGBA2_RGB3_RGBA0_RGBA1_right = _mm_shuffle_epi32(RGBA0_RGBA1_RGBA2_RGB3_right, _MM_SHUFFLE(1, 0, 3, 2));
      __m256i RGBA2_RGBA3_right = _mm256_cvtepu8_epi32(RGBA2_RGB3_RGBA0_RGBA1_right); // convert uint8 to int32 (lower 128 lane)
      __m256 RGBA0_RGBA1_right_f = _mm256_cvtepi32_ps(RGBA0_RGBA1_right); // convert int32 to float
      __m256 RGBA2_RGBA3_right_f = _mm256_cvtepi32_ps(RGBA2_RGBA3_right); // convert int32 to float

      // work view
      __m256i RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_wv = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(wv.i + qy*cols*4 + qx*4)); // unaligned load needed
      __m128i RGBA0_RGBA1_RGBA2_RGB3_wv = _mm256_castsi256_si128(RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_wv);
      __m256i RGBA0_RGBA1_wv = _mm256_cvtepu8_epi32(RGBA0_RGBA1_RGBA2_RGB3_wv); // convert uint8 to int32
      __m128i RGBA2_RGB3_RGBA0_RGBA1_wv = _mm_shuffle_epi32(RGBA0_RGBA1_RGBA2_RGB3_wv, _MM_SHUFFLE(1, 0, 3, 2));
      __m256i RGBA2_RGBA3_wv = _mm256_cvtepu8_epi32(RGBA2_RGB3_RGBA0_RGBA1_wv); // convert uint8 to int32
      __m256 RGBA0_RGBA1_wv_f = _mm256_cvtepi32_ps(RGBA0_RGBA1_wv); // convert int32 to float
      __m256 RGBA2_RGBA3_wv_f = _mm256_cvtepi32_ps(RGBA2_RGBA3_wv); // convert int32 to float

      // wm = [wm0, wm1, wm2, ..., wm7]
      // Goal: [wm0, wm0, wm0, wm0, wm1, wm1, wm1, wm1], [wm2, ..., wm3, ...]
      __m256 inv_wm00001111 = _mm256_permutevar8x32_ps(inv_fac, perm0); // AVX2
      __m256 wm00001111 = _mm256_sub_ps(ones, inv_wm00001111);
      __m256 inv_wm22223333 = _mm256_permutevar8x32_ps(inv_fac, perm1); // AVX2
      __m256 wm22223333 = _mm256_sub_ps(ones, inv_wm22223333);


      // float vp1 =  wm * ov.i[qdy*cols*4+qdx*4 + 0] + inv_fac * ov.i[qdy*cols*4+qdx*4 + 4] - wv.i[qy*cols*4+qx*4 + 0];
      // float vp2 =  wm * ov.i[qdy*cols*4+qdx*4 + 1] + inv_fac * ov.i[qdy*cols*4+qdx*4 + 5] - wv.i[qy*cols*4+qx*4 + 1];
      // float vp3 =  wm * ov.i[qdy*cols*4+qdx*4 + 2] + inv_fac * ov.i[qdy*cols*4+qdx*4 + 6] - wv.i[qy*cols*4+qx*4 + 2];
      __m256 vp_p0_1 = _mm256_fmadd_ps(wm00001111, RGBA0_RGBA1_left_f, _mm256_fmsub_ps(inv_wm00001111, RGBA0_RGBA1_right_f, RGBA0_RGBA1_wv_f));
      __m256 vp_p2_3 = _mm256_fmadd_ps(wm22223333, RGBA2_RGBA3_left_f, _mm256_fmsub_ps(inv_wm22223333, RGBA2_RGBA3_right_f, RGBA2_RGBA3_wv_f));

      // float iqnorm = std::abs(vp1) + std::abs(vp2) + std::abs(vp3);
      // weighted intensity difference of rgba, absolute value
      // [wr0,wg0,wb0,wa0, wr1,wg1,wb1,wa1]
      __m256 abs_vp_p0_1 = _mm256_andnot_ps(minuszeros, vp_p0_1);
      // [wr2,wg2,wb2,wa2, wr3,wg3,wb3,wa3]
      __m256 abs_vp_p2_3 = _mm256_andnot_ps(minuszeros, vp_p2_3);
      // sum up
      // [wrg0,wba0, wrg2,wba2, wrg1,wba1, wrg3,wba3]
      __m256 abs_vp_p0213 = _mm256_hadd_ps(abs_vp_p0_1, abs_vp_p2_3);

      /// repeated structure from above (only index changes)
      // pixel 4-7
      // [R4,G4,B4,A4, R5,G5,B5,A5, R6,G6,B6,A6, R7,G7,B7,A7]
      // left: left matching point for weighting
      // scale = 4 since RGBA uses 4 bytes
      __m128i RGBA4_RGBA5_RGBA6_RGB7_left = _mm256_extracti128_si256(RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_left, 0x1); // AVX2, extract upper 128 lane
      __m256i RGBA4_RGBA5_left = _mm256_cvtepu8_epi32(RGBA4_RGBA5_RGBA6_RGB7_left); // convert uint8 to int32 (lower 128 lane)
      __m128i RGBA6_RGB7_RGBA4_RGBA5_left = _mm_shuffle_epi32(RGBA4_RGBA5_RGBA6_RGB7_left, _MM_SHUFFLE(1, 0, 3, 2));
      __m256i RGBA6_RGBA7_left = _mm256_cvtepu8_epi32(RGBA6_RGB7_RGBA4_RGBA5_left); // convert uint8 to int32 (lower 128 lane)
      __m256 RGBA4_RGBA5_left_f = _mm256_cvtepi32_ps(RGBA4_RGBA5_left); // convert int32 to float
      __m256 RGBA6_RGBA7_left_f = _mm256_cvtepi32_ps(RGBA6_RGBA7_left); // convert int32 to float

      // right: right matching point for weighting
      // scale = 4 since RGBA uses 4 bytes
      __m128i RGBA4_RGBA5_RGBA6_RGB7_right = _mm256_extracti128_si256(RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_right, 0x1); // AVX2, extract upper 128 lane
      __m256i RGBA4_RGBA5_right = _mm256_cvtepu8_epi32(RGBA4_RGBA5_RGBA6_RGB7_right); // convert uint8 to int32 (lower 128 lane)
      __m128i RGBA6_RGB7_RGBA4_RGBA5_right = _mm_shuffle_epi32(RGBA4_RGBA5_RGBA6_RGB7_right, _MM_SHUFFLE(1, 0, 3, 2));
      __m256i RGBA6_RGBA7_right = _mm256_cvtepu8_epi32(RGBA6_RGB7_RGBA4_RGBA5_right); // convert uint8 to int32 (lower 128 lane)
      __m256 RGBA4_RGBA5_right_f = _mm256_cvtepi32_ps(RGBA4_RGBA5_right); // convert int32 to float
      __m256 RGBA6_RGBA7_right_f = _mm256_cvtepi32_ps(RGBA6_RGBA7_right); // convert int32 to float

      // work view
      __m128i RGBA4_RGBA5_RGBA6_RGB7_wv = _mm256_extracti128_si256(RGBA0_RGBA1_RGBA2_RGB3_RGBA4_RGBA5_RGBA6_RGB7_wv, 0x1); // AVX2, extract upper 128 lane
      __m256i RGBA4_RGBA5_wv = _mm256_cvtepu8_epi32(RGBA4_RGBA5_RGBA6_RGB7_wv); // convert uint8 to int32
      __m128i RGBA6_RGB7_RGBA4_RGBA5_wv = _mm_shuffle_epi32(RGBA4_RGBA5_RGBA6_RGB7_wv, _MM_SHUFFLE(1, 0, 3, 2));
      __m256i RGBA6_RGBA7_wv = _mm256_cvtepu8_epi32(RGBA6_RGB7_RGBA4_RGBA5_wv); // convert uint8 to int32
      __m256 RGBA4_RGBA5_wv_f = _mm256_cvtepi32_ps(RGBA4_RGBA5_wv); // convert int32 to float
      __m256 RGBA6_RGBA7_wv_f = _mm256_cvtepi32_ps(RGBA6_RGBA7_wv); // convert int32 to float

      // wm = [wm0, wm1, wm2, ..., wm7]
      // Goal: [wm4, wm4, wm4, wm4, wm5, wm5, wm5, wm5], [wm6, ..., wm7, ...]
      __m256 inv_wm44445555 = _mm256_permutevar8x32_ps(inv_fac, perm2); // AVX2
      __m256 wm44445555 = _mm256_sub_ps(ones, inv_wm44445555);
      __m256 inv_wm66667777 = _mm256_permutevar8x32_ps(inv_fac, perm3); // AVX2
      __m256 wm66667777 = _mm256_sub_ps(ones, inv_wm66667777);

      // float vp1 =  wm * ov.i[qdy*cols*4+qdx*4 + 0] + inv_fac * ov.i[qdy*cols*4+qdx*4 + 4] - wv.i[qy*cols*4+qx*4 + 0];
      // float vp2 =  wm * ov.i[qdy*cols*4+qdx*4 + 1] + inv_fac * ov.i[qdy*cols*4+qdx*4 + 5] - wv.i[qy*cols*4+qx*4 + 1];
      // float vp3 =  wm * ov.i[qdy*cols*4+qdx*4 + 2] + inv_fac * ov.i[qdy*cols*4+qdx*4 + 6] - wv.i[qy*cols*4+qx*4 + 2];
      __m256 vp_p4_5 = _mm256_fmadd_ps(wm44445555, RGBA4_RGBA5_left_f, _mm256_fmsub_ps(inv_wm44445555, RGBA4_RGBA5_right_f, RGBA4_RGBA5_wv_f));
      __m256 vp_p6_7 = _mm256_fmadd_ps(wm66667777, RGBA6_RGBA7_left_f, _mm256_fmsub_ps(inv_wm66667777, RGBA6_RGBA7_right_f, RGBA6_RGBA7_wv_f));

      // float iqnorm = std::abs(vp1) + std::abs(vp2) + std::abs(vp3);
      // weighted intensity difference of rgba, absolute value
      // [wr4,wg4,wb4,wa4, wr5,wg5,wb5,wa5]
      __m256 abs_vp_p4_5 = _mm256_andnot_ps(minuszeros, vp_p4_5);
      // [wr6,wg6,wb6,wa6, wr7,wg7,wb7,wa7]
      __m256 abs_vp_p6_7 = _mm256_andnot_ps(minuszeros, vp_p6_7);
      // sum up
      // [wrg4,wba4, wrg6,wba6, wrg5,wba5, wrg4,wba4]
      __m256 abs_vp_p4657 = _mm256_hadd_ps(abs_vp_p4_5, abs_vp_p6_7);

      // final sum up
      // [wrgba0, wrgba2, wrgba4, wrgba6, wrgba1, wrgba3, wrgba5, wrgba7]
      __m256 abs_vp_p02461357 = _mm256_hadd_ps(abs_vp_p0213, abs_vp_p4657);
      // shuffle in 01234567 order
      __m256 iqnorm = _mm256_permutevar8x32_ps(abs_vp_p02461357, perm_sum);
      /// End:
      /// l1 norm computation (maximum 3/4 usage of simd lanes)
      ///

      // gradient intensity difference between this and other view
      // float ovg  = wm*ov.g[qdy*cols+qdx] + inv_fac*ov.g[qdy*cols+qdx+1];
      __m256i qdy_cols_qdx_left = _mm256_add_epi32(qdycols, qdx);
      __m256i qdy_cols_qdx_right = _mm256_add_epi32(qdy_cols_qdx_left, onesi);
      __m256 ovg_left = _mm256_mask_i32gather_ps(zerosf, ov.g, qdy_cols_qdx_left, mask, 4);
      __m256 ovg_right = _mm256_mask_i32gather_ps(zerosf, ov.g, qdy_cols_qdx_right, mask, 4);
      __m256 ovg = _mm256_fmadd_ps(wm, ovg_left, _mm256_mul_ps(inv_fac, ovg_right));

      // float iqgnorm = std::abs(wv.g[qy*cols+qx] - ovg);
      __m256 wvg = _mm256_loadu_ps(wv.g + qy*cols + qx);  // potentially loads too much data -> use mask
      wvg = _mm256_and_ps(wvg, mask);  // masking
      __m256 diff_g = _mm256_sub_ps(wvg, ovg);
      __m256 iqgnorm = _mm256_andnot_ps(minuszeros, diff_g);  // absolute value

      // cost_img_diff += w * std::min(iqnorm,TAUCOL);
      // cost_grad_diff += w * std::min(iqgnorm,TAUGRAD);
      __m256 clamp_iqnorm = _mm256_min_ps(iqnorm, taucol);
      // TODO: maybe mask taucol, less dependencies
      clamp_iqnorm = _mm256_and_ps(clamp_iqnorm, mask);
      __m256 clamp_iqgnrom = _mm256_min_ps(iqgnorm, taugrad);
      // TODO: maybe mask taugrad, less dependencies
      clamp_iqgnrom = _mm256_and_ps(clamp_iqgnrom, mask);
      cost_img_diff = _mm256_fmadd_ps(w, clamp_iqnorm, cost_img_diff);
      cost_grad_diff = _mm256_fmadd_ps(w, clamp_iqgnrom, cost_grad_diff);

      match_unclamped = _mm256_add_ps(match_unclamped, match_increase);
      qxis = _mm256_add_epi32(qxis, eightsi);
      idx += 8; // ++idx
    }
    disp_tmp += fp[1];
    qdycols = _mm256_add_epi32(qdycols, mcols);
  }
  // accumulators
  // reduce to mm256 array
  __m256 iiggiigg = _mm256_hadd_ps(cost_img_diff, cost_grad_diff);
  __m256 weighting = _mm256_set_ps(ALPHA, ALPHA, ONEMINUSALPHA, ONEMINUSALPHA, ALPHA, ALPHA, ONEMINUSALPHA, ONEMINUSALPHA);
  iiggiigg = _mm256_mul_ps(iiggiigg, weighting);
  // use: https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86
  __m128 vlow  = _mm256_castps256_ps128(iiggiigg);
  __m128 vhigh = _mm256_extractf128_ps(iiggiigg, 1); // high 128
         vlow  = _mm_add_ps(vlow, vhigh);     // add the low 128
  // hsum_ps_sse3
  __m128 shuf = _mm_movehdup_ps(vlow);        // broadcast elements 3,1 to 2,0
  __m128 sums = _mm_add_ps(vlow, shuf);
  shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
  sums        = _mm_add_ss(sums, shuf);
#ifdef TIME_MCOST
  myInt64 end = stop_tsc(start);
  pm::mcost_total_time_ += end;
  pm::mcost_calls_ += 1;
#endif
  
  return _mm_cvtss_f32(sums);
}

void precompute_weights(KernelView& wv, int x, int y, int rows, int cols,
                        int qy_start, int qy_end, int qx_start, int qx_end) {
  int idx = 0;
  for(int qy = qy_start; qy <= qy_end; ++qy){
    for(int qx = qx_start; qx <= qx_end; ++qx){
      // Weight between p and q
      int inorm = l1norm_naive(&wv.i[y*cols*4+x*4], &wv.i[qy*cols*4+qx*4]);
      weights[idx] = fast_exp(inorm);
      ++idx;
    }
  }
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
      precompute_weights(wv, x, y, rows, cols, qy_start, qy_end, qx_start, qx_end);
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
  precompute_weights(wv, x, y, rows, cols, qy_start, qy_end, qx_start, qx_end);
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
    // in the CommonView RGB format is used here we use the RGBA format for easier SIMD processing
    v1_.i = static_cast<uint8_t *>(aligned_alloc(16,(rows*cols*4 + 16)*sizeof(uint8_t)));
    v2_.i = static_cast<uint8_t *>(aligned_alloc(16,(rows*cols*4 + 16)*sizeof(uint8_t)));
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        v1_.i[y * cols * 4 + x * 4 + 0] = v1.i[y * cols * 3 + x * 3 + 0];
        v1_.i[y * cols * 4 + x * 4 + 1] = v1.i[y * cols * 3 + x * 3 + 1];
        v1_.i[y * cols * 4 + x * 4 + 2] = v1.i[y * cols * 3 + x * 3 + 2];
        v1_.i[y * cols * 4 + x * 4 + 3] = 0;
      }
    }
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        v2_.i[y * cols * 4 + x * 4 + 0] = v2.i[y * cols * 3 + x * 3 + 0];
        v2_.i[y * cols * 4 + x * 4 + 1] = v2.i[y * cols * 3 + x * 3 + 1];
        v2_.i[y * cols * 4 + x * 4 + 2] = v2.i[y * cols * 3 + x * 3 + 2];
        v2_.i[y * cols * 4 + x * 4 + 3] = 0;
      }
    }
    // Gradients
    v1_.g = static_cast<float *>(aligned_alloc(16, (rows*cols + 8)*sizeof(float)));
    v2_.g = static_cast<float *>(aligned_alloc(16, (rows*cols + 8)*sizeof(float)));
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        v1_.g[y * cols + x] = v1.g[y * cols + x];
        v2_.g[y * cols + x] = v2.g[y * cols + x];
      }
    }
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
    // precompute exp values
    setup_exp();
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
    precompute_weights(v1_, x, y, rows_, cols_, qy_start, qy_end, qx_start, qx_end);
    float cost_left = mcost(v1_, v2_, fp, x, y, rows_, cols_, -1, qy_start, qy_end, qx_start, qx_end);
    fp = &(v2_.p[(y*cols_*3)+x*3]);
    precompute_weights(v2_, x, y, rows_, cols_, qy_start, qy_end, qx_start, qx_end);
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
    delete[] v1_.i;
    delete[] v2_.i;
    delete[] v1_.g;
    delete[] v2_.g;
    delete[] v1_.p;
    delete[] v2_.p;
  }

private:
  KernelView v1_;
  KernelView v2_;
};

}//namespace simd_v8
}//namespace pm
