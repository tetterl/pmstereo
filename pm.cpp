#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <fstream>
#include <list>

#include "base_kernel.hpp"
#include "parameters.hpp"
#include "png.hpp"
#include "kernel_registration.hpp"
#include "eval.hpp"
#include "tsc_x86.h"

using namespace std;

int main(int argc, char** argv) {
  std::string kernel("baseline");
  std::string path1, path2, outputDirectory=".";
  std::string listkernels = "";
  std::string gt_path = "";
  std::string mask_path = "";

  // extract the binaries filename (we'll need it later for the timing logs)
  std::string path = string(argv[0]);
  std::string base_filename = path.substr(path.find_last_of("/\\") + 1);

  // Single precision peak performance (flops per cycle)
  float sp_pp = 32.;


  // Get args
  std::vector<std::string> args = { argv, argv + argc };
  bool unit_test = true;
  int arridx=0;
  int numRequiredArgs=2; // (left and right image paths are required args)
  // Exits if argument list doesn't contain arg after named arg
  auto argcheck = [&] () {
    if(args.size() < arridx + 2){
      //cout << args.size() << "  " << arridx << endl;
      cout << "ERR: malformed argument input, aborting" << endl;
      cout  << "./pm --left img1 --right img2 --output . --kernel baseline --maxdisp 64 --unit-test 1 --gt groundtruth --mask maskfile" << endl;
      cout << "Order of arguments can be arbitrary, --left and --right are required args" << endl;
      exit(EXIT_FAILURE);
    }
  };
  for (auto && par : args) {
    if (par == "--left"){
      argcheck();
      path1 = args[arridx+1];
      numRequiredArgs--;
    }
    // This is a special mode which serves to list all available kernels
    // for automatic timing measurements by an external script
    else if(par == "--listkernels"){
      argcheck();
      listkernels = args[arridx+1];
    }
    else if (par == "--right"){
      argcheck();
      path2 = args[arridx+1];
      numRequiredArgs--;
    }
    else if (par == "--output"){
      argcheck();
      outputDirectory = args[arridx+1];
    }
    else if (par == "--kernel"){
      argcheck();
      kernel = args[arridx+1];
    }
    else if (par == "--maxdisp"){
      argcheck();
      max_disp =stoi(args[arridx+1]);
    }
    else if (par == "--unit-test"){
      argcheck();
      unit_test = stoi(args[arridx+1]);
    }
    else if (par == "--sppp"){
      argcheck();
      sp_pp = float(stoi(args[arridx+1]));
    }
    // supplying groundtruth does evaluation and writes 
    // the resulting scores into the same file as the timings.
    else if (par == "--gt"){
      argcheck();
      gt_path = args[arridx+1];
    }
    else if(par == "--mask"){
      argcheck();
      mask_path = args[arridx+1];
    }
    arridx++;
  }
  // Ensure that both mask and ground truth are provided 
  // for evaluation
  if(gt_path != "" && mask_path == "" 
      || gt_path == "" && mask_path != ""){
    cout << "ERR: For evaluation within ./pm, you need to specificy both --mask AND --gt. aborting." << endl;
    exit(EXIT_FAILURE);
  }
  if(listkernels != ""){
    KernelIterator ki;

    for(ki = kernel_map.begin(); ki != kernel_map.end(); ki++){
      cout << ki->first << ",";
    }
    cout << endl;
    exit(EXIT_SUCCESS);
  }
  if(numRequiredArgs > 0){
    cout << "ERR: --left and --right image paths are required, aborting!" << endl;
    exit(EXIT_FAILURE);
  }


  Image img1 = readPNG(path1);
  Image img2 = readPNG(path2);

  if (img1.width == -1 || img1.height == -1 || img2.width == -1 || img2.height == -1) {
    // Error message is already printed for png reader
    return 1;
  }

  if(!( img1.width == img2.width && img1.height == img2.height && img1.channels == img2.channels)){
    cout << "ERR: image dimensions mismatch!" << endl;
    return 1;
  }

  int rows = img1.height;
  int cols = img1.width;
  cout << "PM: imgs loaded" << endl;
  cout << "PM: dimensions: " << cols << " x " << rows << endl;


  // Structs with both views (each have image, gradient, plane and disparity)
  CommonView v1, v2;
  PreProcessing(img1, img2, v1, v2, rows, cols, outputDirectory);

  /// unit-test
  if (unit_test) {
      float rel_tol = 2e-3;
      cout << "PM: ###### unit-testing of all kernels ######" << endl;
    // pixels to test
    std::vector<std::pair<int, int>> test_pixels;
    test_pixels.reserve(cols*rows);
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++ x) {
        test_pixels.push_back({y,x});
      }
    }

    // reference solution
    pm::BaseKernel* baseline_kernel = kernel_map.at("baseline")(v1, v2, rows, cols);
    std::vector<std::pair<float, float>> baseline_cost; // pair for left and right image as working view
    for (auto& test_pixel : test_pixels) {
      baseline_cost.push_back(baseline_kernel->test_mcost(test_pixel.second, test_pixel.first));
    }

    // comparison with other kernels
    for (auto& kernel_creator : kernel_map) {
      cout << "PM: test kernel: " << kernel_creator.first<< endl;
      pm::BaseKernel* test_kernel = kernel_creator.second(v1, v2, rows, cols);
      std::vector<std::pair<float, float>> test_kernel_cost; // pair for left and right image as working view
      for (auto& test_pixel : test_pixels) {
        test_kernel_cost.push_back(test_kernel->test_mcost(test_pixel.second, test_pixel.first));
      }
      bool failed = false;
      for (int i = 0; i < test_pixels.size(); ++i) {
        // left working view
        float diff_left = std::abs(baseline_cost[i].first - test_kernel_cost[i].first);
        bool fail_left = diff_left / baseline_cost[i].first > rel_tol;
        if (fail_left) {
          cout << "PM: \033[1;31mFAILED\033[0m for left view of kernel: " << kernel_creator.first
               << ", with pixel: {y,x}={" << test_pixels[i].first << "," << test_pixels[i].second << "}" << endl;
          cout << "baseline_cost= " << baseline_cost[i].first << ", this_cost=" << test_kernel_cost[i].first << endl;
        }
        float diff_right = std::abs(baseline_cost[i].second - test_kernel_cost[i].second);
        bool fail_right = diff_right / baseline_cost[i].second > rel_tol;
        if (fail_right) {
          cout << "PM: \033[1;31mFAILED\033[0m for right view of kernel: " << kernel_creator.first
               << ", with pixel: {y,x}={" << test_pixels[i].first << "," << test_pixels[i].second << "}" << endl;
          cout << "baseline_cost= " << baseline_cost[i].second << ", this_cost=" << test_kernel_cost[i].second << endl;
        }
        failed += fail_left + fail_right;
      }
      if (!failed) {
        cout << "PM: \033[1;32mPASSED\033[0m for kernel: " << kernel_creator.first << endl;
      }
    }
    cout << "PM: ###### unit-testing of all kernels finished, " << test_pixels.size() << " points tested ######" << endl;
  }
  
  // perf-test (fast and thus run in any case)
  float perf_w, perf_q, flop_per_cycle;
  {
    // Get a kernel with given view
    pm::BaseKernel* peakperf_kernel = kernel_map.at(kernel)(v1, v2, rows, cols);
    
    // Setup common views with deterministic plane values
    int num_runs = 10;
    double cycles = 0.;
    int REP=100;
    myInt64 start, end;
    start = start_tsc();
    for (size_t i = 0; i < num_runs; i++) {
        peakperf_kernel->peakperf_mcost(); //initially hardcoded.
    }
    end = stop_tsc(start);

    cycles = (double)end;
    
    // Actual measurement
    list<double> cyclesList;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
          peakperf_kernel->peakperf_mcost(); //initially hardcoded.
        }
        end = stop_tsc(start);
        cycles = ((double)end) / num_runs;
        cyclesList.push_back(cycles);
    }

    cyclesList.sort();
    cycles = cyclesList.front();    
    perf_w = peakperf_kernel->get_W_mcost();
    perf_q = peakperf_kernel->get_Q_mcost();
    flop_per_cycle =  perf_w / cycles;

    cout << "cycles: " << cycles << " W: " << perf_w << endl;
    cout << "floppercycle:" << flop_per_cycle << endl;
    cout << "I:" << perf_w/perf_q << endl;
  }

  // The kernel we are using
  cout << "PM: kernel " << kernel << endl;

  /// ------------ Core start ------------
  pm::BaseKernel* active_kernel;
  try {
      active_kernel = kernel_map.at(kernel)(v1, v2, rows, cols);
  } catch (const std::out_of_range& e) {
      cout << "ERR: specified kernel not found. Have you registered it?" << endl;
      return 1;
  }
  active_kernel->reset_mcost_counters();

  auto start = std::chrono::high_resolution_clock::now();
  active_kernel->run_patch_match();
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end-start;
  double duration = diff.count();
  cout << "PM: kernel " << kernel << endl;
  cout << "PM: time for process = " << duration << "s" << endl;
#ifdef TIME_MCOST
  cout << "PM: mcost num calls: " << active_kernel->get_mcost_calls() << endl;
  cout << "PM: mcost total time: " << active_kernel->get_mcost_total_time() << endl;
  double avg_mcost_time = double(active_kernel->get_mcost_total_time())/double(active_kernel->get_mcost_calls());
  cout << "PM: avg mcost time: " << avg_mcost_time  << endl;
#else
  std::string avg_mcost_time = "";
#endif

  active_kernel->update_common_view(v1, v2);
  delete active_kernel;
  /// ------------ Core end ------------


  cout << "PM: Computing disparities" << endl;
  // We need to do this once before PostProcessing
  PlanesToDisparity(v1, rows, cols);
  PlanesToDisparity(v2, rows, cols);

  Image disp1 = Image(rows,cols,1,8, v1.d);
  Image disp2 = Image(rows,cols,1,8, v2.d);
  writePNG(disp1, outputDirectory + "/disp1.png");
  writePNG(disp2, outputDirectory + "/disp2.png");

  // Do post processing.
  // This operates directly on the disparities and does not update the planes.
  // --> Do not call PlanesToDisparity again after this!
  cout << "PM: Post processing" << endl;
  PostProcessing(v1, v2, rows, cols);

  Image disp1_clean = Image(rows,cols,1,8, v1.d);
  Image disp2_clean = Image(rows,cols,1,8, v2.d);

  // Discrete disparities for visual inspection
  writePNG(disp1_clean, outputDirectory + "/disp1_clean.png");
  writePNG(disp2_clean, outputDirectory + "/disp2_clean.png");

  // Floating point disparities for evaluation
  // By default, store it in the build folder.
  // else save to given path.
  // We only store the left disparity because
  // the ground truth does not offer the right disparity
  WriteFilePFM(v1, (outputDirectory + "/disp1.pfm").c_str(), rows,cols);
  WriteFilePFM(v2, (outputDirectory + "/disp2.pfm").c_str(), rows,cols);
  
  // Write timings and and (optionally) evaluation to file (append)
  ofstream evalFile;
  cout << "PM: writing evaluation file" << endl;
  evalFile.open (outputDirectory + "/evaluation.txt", std::ofstream::out | std::ofstream::app);
  // Write timings + evaluation
  if(gt_path != "" && mask_path != ""){
    int rows_gt, cols_gt;
    float* gt = middlebury::ReadFilePFM(gt_path.c_str(), &rows_gt, &cols_gt);

    Image mask = readPNG(mask_path);
    middlebury::StereoScore score = middlebury::evaldisp(v1.d, gt, mask, rows, cols, 64, 1);
    evalFile << base_filename << "," << kernel << "," <<  rows << "," << cols << "," << duration  << ","
      <<  perf_w << "," << perf_q << "," << flop_per_cycle << "," << avg_mcost_time << "," << sp_pp << ","
      << score.coverage << "," << score.bad05 << "," << score.bad1 << "," 
      << score.invalid << "," << score.avgErr <<  endl;
    delete[] gt;
  }else{ // Write timings only
    evalFile << base_filename << "," << kernel << "," <<  rows << "," << cols << "," << duration  << ","
      <<  perf_w << "," << perf_q << "," << flop_per_cycle << "," << avg_mcost_time << "," << sp_pp
      << endl;
  }
  evalFile.close();

  cout << "PM: cleaning up, deleting image arrays" << endl;
  deallocate_common_view(v1);
  deallocate_common_view(v2);
  
  

  return 0;
}
