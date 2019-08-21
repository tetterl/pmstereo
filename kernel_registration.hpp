#ifndef KERNEL_REGISTRATION_HPP
#define KERNEL_REGISTRATION_HPP

#include <functional>
#include <map>
#include <string>

#include "base_kernel.hpp"
// include all kernel files
#include "kernels-baseline.hpp"
#include "kernels-compressed_planes.hpp"
#include "kernels-removed_inside_and_other.hpp"
#include "kernels-weight_precomputation.hpp"
#include "kernels-fastexp-instr.hpp"
#include "kernels-fastexp-array.hpp"
#include "kernels-fastexp-array-declaration.hpp"
#include "kernels-fast_exp.hpp"
#include "kernels-weight_precomputation_fast_exp.hpp"
#include "kernels-split-up-cost.hpp"
#include "kernels-remove_if.hpp"
#include "kernels-inline_l1_norm.hpp"
#include "kernels-rgba.hpp"
#include "kernels-simd_v0.hpp"
#include "kernels-simd_v1.hpp"
#include "kernels-simd_v2.hpp"
#include "kernels-simd_v3.hpp"
#include "kernels-simd_v4.hpp"
#include "kernels-simd_v6.hpp"
#include "kernels-simd_v7.hpp"
#include "kernels-simd_v8.hpp"
#include "kernels-simd_v9.hpp"
#include "kernels-simd_v10.hpp"
#include "kernels-simd_v11.hpp"
#include "kernels-simd_v12.hpp"
#include "kernels-simd_v13.hpp"
#include "kernels-simd_v14.hpp"
#include "kernels-simd_v15.hpp"
#include "kernels-simd_v16.hpp"
#include "kernels-simd_v17.hpp"
#include "kernels-simd_v18.hpp"
#include "kernels-simd_v21.hpp"
#include "kernels-simd_v22.hpp"
#include "kernels-simd_v23.hpp"

// global kernel factory
#define SIGNATURE_KERNEL_CONSTRUCTOR const CommonView& v1, const CommonView& v2, int rows, int cols
std::map<std::string, std::function<pm::BaseKernel*(SIGNATURE_KERNEL_CONSTRUCTOR)>> kernel_map{
    // include all kernel implementations
    {"baseline", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::baseline::Kernel(v1, v2, rows, cols);}},
    {"compressed_planes", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::compressed_planes::Kernel(v1, v2, rows, cols);}},
    {"removed_inside_and_other", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::removed_inside_and_other::Kernel(v1, v2, rows, cols);}},
    {"split_up_cost", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::split_up_cost::Kernel(v1, v2, rows, cols);}},
    {"weight_precomputation", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::weight_precomputation::Kernel(v1, v2, rows, cols);}},
    {"fast_exp", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::fast_exp::Kernel(v1, v2, rows, cols);}},
    {"weight_precomputation_fast_exp", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::weight_precomputation_fast_exp::Kernel(v1, v2, rows, cols);}},
    {"remove_if", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::remove_if::Kernel(v1, v2, rows, cols);}},
    {"inline_l1_norm", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::inline_l1_norm::Kernel(v1, v2, rows, cols);}},
    {"rgba", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::rgba::Kernel(v1, v2, rows, cols);}},
    {"simd_v0", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v0::Kernel(v1, v2, rows, cols);}},
    {"simd_v1", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v1::Kernel(v1, v2, rows, cols);}},
    {"simd_v2", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v2::Kernel(v1, v2, rows, cols);}},
    {"simd_v3", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v3::Kernel(v1, v2, rows, cols);}},
    {"simd_v4", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v4::Kernel(v1, v2, rows, cols);}},
    {"simd_v6", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v6::Kernel(v1, v2, rows, cols);}},
    {"simd_v7", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v7::Kernel(v1, v2, rows, cols);}},
    {"simd_v8", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v8::Kernel(v1, v2, rows, cols);}},
    {"simd_v9", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v9::Kernel(v1, v2, rows, cols);}},
    {"simd_v10", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v10::Kernel(v1, v2, rows, cols);}},
    {"simd_v11", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v11::Kernel(v1, v2, rows, cols);}},
    {"simd_v12", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v12::Kernel(v1, v2, rows, cols);}},
    {"simd_v13", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v13::Kernel(v1, v2, rows, cols);}},
    {"simd_v14", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v14::Kernel(v1, v2, rows, cols);}},
    {"simd_v15", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v15::Kernel(v1, v2, rows, cols);}},
    {"simd_v16", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v16::Kernel(v1, v2, rows, cols);}},
    {"simd_v17", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v17::Kernel(v1, v2, rows, cols);}},
    {"simd_v18", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v18::Kernel(v1, v2, rows, cols);}},
    {"simd_v21", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v21::Kernel(v1, v2, rows, cols);}},
    {"simd_v22", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v22::Kernel(v1, v2, rows, cols);}},
    {"simd_v23", [](SIGNATURE_KERNEL_CONSTRUCTOR) {return new pm::simd_v23::Kernel(v1, v2, rows, cols);}},
  };

// Iterator for the kernel map(used to list all available methods)
typedef std::map<std::string, std::function<pm::BaseKernel*(SIGNATURE_KERNEL_CONSTRUCTOR)>>::iterator KernelIterator;
#endif // KERNEL_REGISTRATION_HPP
