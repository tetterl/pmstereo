/*    Adapted from MiddEval (middlebury stereo evaluation toolkit)
 *    Original Authro Daniel Scharstein
 *    GNU Licence
 */
#ifndef MIDDLEBURY_EVAL
#define MIDDLEBURY_EVAL
#include <iomanip>
#include <cmath>
#include "png.hpp"
using namespace std;
namespace middlebury{
class StereoScore{
  public:
    float coverage;
    float bad05;
    float bad1;
    float invalid;
    float avgErr;
    StereoScore(float coverage, float bad05, float bad1, float invalid, float avgErr) :
      coverage(coverage), bad05(bad05), bad1(bad1), invalid(invalid), avgErr(avgErr){}
}; 
void skip_space(FILE *fp) {
    char c;
    do {
        c = getc(fp);
    } while (c == '\n' || c == ' ' || c == '\t' || c == '\r');
    ungetc(c, fp);
}
void skip_comment(FILE *fp) {
    char c;
    while ((c=getc(fp)) == '#')
        while (getc(fp) != '\n') ;
    ungetc(c, fp);
}
int littleendian() {
    int intval = 1;
    uint8_t *uval = (uint8_t *)&intval;
    return uval[0] == 1;
}
StereoScore evaldisp(float* disp,  float* gt_dat, Image& mask, int rows, int cols, int maxdisp, int rounddisp){
  int n = 0;
  int bad1 = 0;
  int bad05 = 0;
  int invalid = 0;
  float serr = 0;
  int width  = cols;
  int height = rows;
  for (int y = 0; y < height; y++) {
  	for (int x = 0; x < width; x++) {
      float gt = gt_dat[width*y+x];
			//float gt = gtdisp.Pixel(x, y, 0);
	    if (gt == INFINITY) // unknown
    		continue;
      float d = disp[y*width+x];
	    int valid = (d != INFINITY);
	    if (valid) {
	     	float maxd = maxdisp; // max disp range
		    d = std::max(0.f, std::min(maxd, d)); // clip disps to max disp range
	    }
	    if (valid && rounddisp)
		    d = std::round(d);

	    float err = std::abs(d - gt);
	    if (mask.data[y*width+x] != 255) { // don't evaluate pixel
      } else {
        n++;
        if (valid) {
          serr += err;
          if (err > 1.) {
            bad1++;
          }
          if (err > 0.5) {
            bad05++;
          }
        } else {// invalid (i.e. hole in sparse disp map)
          invalid++;
        }
	    }//
    }//inner for
  }//outer for
  float bad1percent =  100.0*bad1/n;
  float bad05percent =  100.0*bad05/n;
  float invalidpercent =  100.0*invalid/n;
  //float totalbadpercent =  100.0*(bad+invalid)/n;
  float avgErr = serr / (n - invalid); // CHANGED 10/14/2014 -- was: serr / n
  //printf("mask  bad%.1f  invalid  totbad   avgErr\n", badthresh);

  StereoScore score(100.0*n/(width*height), bad05percent, 
      bad1percent, invalidpercent, avgErr);
  return score;
}
void read_header(FILE *fp, const char *imtype, char c1, char c2, int *width, int *height, int *nbands, int thirdArg) {
    // read the header of a pnmfile and initialize width and height
    char c;

	if (getc(fp) != c1 || getc(fp) != c2)
    cout << "wrong magic code" << endl;
	skip_space(fp);
	skip_comment(fp);
	skip_space(fp);
	fscanf(fp, "%d", width);
	skip_space(fp);
	fscanf(fp, "%d", height);
	if (thirdArg) {
		skip_space(fp);
		fscanf(fp, "%d", nbands);
	}
    // skip SINGLE newline character after reading image height (or third arg)
	c = getc(fp);
    if (c == '\r')      // <cr> in some files before newline
        c = getc(fp);
    if (c != '\n') {
        if (c == ' ' || c == '\t' || c == '\r')
          cout << "newline expected in file" << endl;
        else
          cout << "white space expected" << endl;
  }
}
float*  ReadFilePFM(const char* filename, int* rows, int* cols) {
  // Open the file and read the header
  FILE *fp = fopen(filename, "rb");
  if (fp == 0)
    cout << "couldn't open file" << endl;

  int width, height, nBands;
  read_header(fp, "PFM", 'P', 'f', &width, &height, &nBands, 0);

  skip_space(fp);

  float scalef;
  fscanf(fp, "%f", &scalef);  // scale factor (if negative, little endian)

  // skip SINGLE newline character after reading third arg
  char c = getc(fp);
  if (c == '\r')      // <cr> in some files before newline
      c = getc(fp);
  if (c != '\n') {
      if (c == ' ' || c == '\t' || c == '\r')
        cout << "ERR: newline expected in file after scale factor"<<endl;
      else
        cout << "whitespace expected in file after scale factor"<<endl;
  }


  float* data = (float*)malloc(width*height*sizeof(float));
  // Set the image shape
  *rows = height;
  *cols = width;


  int littleEndianFile = (scalef < 0);
  int littleEndianMachine = littleendian();
  int needSwap = (littleEndianFile != littleEndianMachine);

  for (int y = height-1; y >= 0; y--) { // PFM stores rows top-to-bottom!!!!
	  int n = width;
		//float* ptr = (float *) img.PixelAddress(0, y, 0);
	  float* ptr = &data[y*width];

	  if ((int)fread(ptr, sizeof(float), n, fp) != n)
      cout << "File is too short" << endl;

	  if (needSwap) { // if endianness doesn't agree, swap bytes
	    uint8_t* ptr = (uint8_t *) &data[y*width];
	    int x = 0;
	    uint8_t tmp = 0;
	    while (x < n) {
        tmp = ptr[0]; ptr[0] = ptr[3]; ptr[3] = tmp;
        tmp = ptr[1]; ptr[1] = ptr[2]; ptr[2] = tmp;
        ptr += 4;
        x++;
	    }
	  }
  }
  if (fclose(fp))
    cout << "ReadFile: error closing file" << endl;
  return data;
}
}//middlebury namespace
#endif
