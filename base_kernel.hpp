#ifndef BASE_KERNEL_HPP
#define BASE_KERNEL_HPP

#include <utility>

#include "helpers.hpp" // needed for init

#include "tsc_x86.h"

namespace pm {
std::size_t mcost_calls_ = 0;
myInt64 mcost_total_time_ = 0;

class BaseKernel {
public:
  // force children to implement destructor
  virtual ~BaseKernel() = 0;
  /**
   * @brief run_patch_match
   */
  virtual void run_patch_match() = 0;
  /**
   * @brief update_common_view
   *
   * only needs to update values a,b,c in planes
   * @param v1
   * @param v2
   */
  virtual void update_common_view(CommonView& v1, CommonView& v2) const = 0;

  /**
   * @brief test_mcost
   *
   * calls mcost on a given pixel for testing
   * @param v1
   * @param v2
   */
  virtual std::pair<float, float> test_mcost(int x, int y) = 0; // TODO: make const

  /**
   * @brief Single call to mcost s.t. it uses it's maximum runtime
   *
   * @return the cost
   */
  virtual float peakperf_mcost(void) = 0;

  /**
   * @brief Return the theoretical bound on the work, 
   *        deterimned by counting the ops.
   *
   * @return  the work for one call to mcost. 
   */
  virtual float get_W_mcost(void)=0;

  /**
   * @brief Return an upper bound of data transferred. 
   *
   * @return the number of bytes transferred
   */
  virtual float get_Q_mcost(void)=0;

  void reset_mcost_counters() {
    mcost_calls_ = 0;
    mcost_total_time_ = 0;
  }
  /**
   * @brief Get the number of times mcost was called
   * @return number of times mcost was called
   */ 
  std::size_t get_mcost_calls() {
    return mcost_calls_;
  }
 
  /**
   * @breif Get the cumulative amount of time to run mcost when mcost does not exit early
   */ 
  myInt64 get_mcost_total_time() {
    return mcost_total_time_;
  }
private:
  // disable constructor
  BaseKernel();
protected:
  // disable constructor
  BaseKernel(int rows, int cols) : rows_(rows), cols_(cols) {}
  int rows_;
  int cols_;


};
BaseKernel::~BaseKernel() {}
}
#endif // BASE_KERNEL_HPP
