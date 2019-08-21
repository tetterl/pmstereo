import sys
import os

if len(sys.argv) != 3:
    print("Usage: " + sys.argv[0] + " <lower_bound> <upper_bound>")
    print("Generates a fast exponential c file for values between the lower and upper bounds inclusive")

with open("fast_exp_instr.cpp", "w") as fast_exp_file:
    fast_exp_file.write("#include <cmath>"+os.linesep)
    fast_exp_file.write("#include \"parameters.hpp\""+os.linesep)
    fast_exp_file.write("namespace pm {"+os.linesep)
    fast_exp_file.write("namespace fastexp_instr {"+os.linesep)
    fast_exp_file.write("float fast_exp(int val) {"+os.linesep)

    fast_exp_file.write("  switch(val) {"+os.linesep)
    for val in range(int(sys.argv[1]), int(sys.argv[2])+1):
        fast_exp_file.write("    case "+str(val)+":"+os.linesep)
        fast_exp_file.write("      return std::exp("+str(val)+"/GAMMA);"+os.linesep)

    fast_exp_file.write("    default:"+os.linesep)
    fast_exp_file.write("      return std::exp(val/GAMMA);"+os.linesep)
    fast_exp_file.write("  }"+os.linesep)
    fast_exp_file.write("}"+os.linesep)
    fast_exp_file.write("} // namespace fastexp_instr"+os.linesep)
    fast_exp_file.write("} // namespace pm"+os.linesep)
