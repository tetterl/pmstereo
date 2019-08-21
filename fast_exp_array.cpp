#include <cmath>
#include "parameters.hpp"
namespace pm {
namespace fastexp_array {

float precomputed_vals[771];
void setup_exp() {
  precomputed_vals[0] = std::exp(0);
  for (int i=1; i<771; i++) {
    precomputed_vals[i] = std::exp(-i/GAMMA);
  }
}
float fast_exp(int val) {
  // Assume value is in [-n,0]
  if (val<=0 || val>=-770) {
      return precomputed_vals[-val];
  }
  return std::exp(val/GAMMA);
}
} // namespace fastexp_array
} // namespace pm
