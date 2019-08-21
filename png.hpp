#ifndef __PNGWRAPPER
#define __PNGWRAPPER
#include <png.h>
#include <vector>
#include <iostream>
using namespace std;

struct Image{
  // Image dimensions
  int height=-1;
  int width=-1;
  // 8 or 16
  int bits=-1;
  // 1 or 3(RGB)
  int channels=-1;
  // Row major pixel data. If RGB, R1G1B1R2G2B2...
  // 8bit data is stored in 16bit type as well.
  // Endianness-handling is up to the user.
  std::vector<uint16_t> data;

  // Ensures images has non-zero dimensions and that
  // the vector length matches those dimensions
  bool isValid(void){
    return (channels > 0 && height > 0 && width > 0 && (bits == 8 || bits == 16) &&
        ( (channels == 1 && data.size() == width * height)
          || (channels == 3 && data.size() == width * height * 3) ));
  }
  Image(int rows, int cols, int channels, int bits) : height(rows), width(cols), channels(channels), bits(bits){
    data.resize(rows*cols*channels);
  }
  Image(int rows, int cols, int channels, int bits, uint8_t* img)
    : height(rows), width(cols), channels(channels), bits(bits){
    data.resize(rows*cols*channels);
    // Copy data
    for(int i=0;i<rows*cols*channels;i++){
      data[i] = img[i];
    }
  }
  Image(int rows, int cols, int channels, int bits, float* img)
    : height(rows), width(cols), channels(channels), bits(bits){
    data.resize(rows*cols*channels);
    // Copy data
    for(int i=0;i<rows*cols*channels;i++){
      data[i] = uint8_t(img[i]);
    }
  }
  Image(){};
};

/**
 * @brief Reads 8 or 16 bit 1-channel and 8bit 3 channel (RGB) PNG.
 *
 * @param filename
 *
 * @return the image
 */
Image readPNG(std::string filename) {
  unsigned char header[8];
  png_byte colorType;
  png_byte bitDepth;

  png_structp pngPtr;
  png_infop infoPtr;
  png_bytep * rowPointers;

  Image img;

  // open file and test for it being a png
  FILE *fp = fopen(filename.c_str(), "rb");
  if (!fp) {
    cout << "ERR: File" << filename << " could not be opened for reading" << endl;
    return img;
  }
  size_t res = fread(header, 1, 8, fp);
  if (png_sig_cmp(header, 0, 8)) {
    cout << "ERR: File" << filename <<  " is not recognized as a PNG file" << endl;
    return img;
  }
  // initialize stuff
  pngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!pngPtr) {
    cout <<  "ERR: png_create_read_struct failed" << endl;
    return img;
  }

  infoPtr = png_create_info_struct(pngPtr);
  if (!infoPtr) {
    cout << "ERR: png_create_info_struct failed" << endl;
    return img;
  }

  if (setjmp(png_jmpbuf(pngPtr))) {
    cout << "ERR: Error during init_io" << endl;
    return img;
  }

  png_init_io(pngPtr, fp);
  png_set_sig_bytes(pngPtr, 8);

  png_read_info(pngPtr, infoPtr);

  img.width = png_get_image_width(pngPtr, infoPtr);
  img.height = png_get_image_height(pngPtr, infoPtr);
  colorType = png_get_color_type(pngPtr, infoPtr);
  img.bits = png_get_bit_depth(pngPtr, infoPtr);

  int numberOfPasses = png_set_interlace_handling(pngPtr);
  png_read_update_info(pngPtr, infoPtr);

  // read file
  if (setjmp(png_jmpbuf(pngPtr))) {
    cout << "ERR: Error during read_image" << endl;
    return img;
  }

  rowPointers = (png_bytep*) malloc(sizeof(png_bytep) * img.height);
  for (int y = 0; y < img.height; y++)
    rowPointers[y] = (png_byte*) malloc(png_get_rowbytes(pngPtr, infoPtr));

  png_read_image(pngPtr, rowPointers);

  fclose(fp);

  switch (png_get_color_type(pngPtr, infoPtr)) {
    case PNG_COLOR_TYPE_GRAY:
       img.channels = 1;
       img.data.resize(img.width*img.height);
       break;
    case PNG_COLOR_TYPE_RGB:
       img.channels = 3;
       img.data.resize(img.width*img.height*3);
       break;
    case PNG_COLOR_TYPE_RGBA:
       img.channels = -1;
       return img;
       break;
    default:
       img.channels=-1;
       return img;
       break;
  }

  if (img.bits == 16){
    if(img.channels != 1) {
      img.channels = -1;
      cout << "ERR: 16bits channels are only supported for single channel images" << endl;
    }else{
      for (int y = 0; y < img.height; y++) {
        png_byte* row = rowPointers[y];
        for (int x = 0; x < img.width; x++) {
          img.data[img.width*y+x]  = ((uint16_t)row[x * 2] << 8) + row[x * 2 + 1];
        }
      }
    }
  } else {
    for (int y = 0; y < img.height; y++) {
      png_byte* row = rowPointers[y];
      for (int x = 0; x < img.width; x++) {
        int offset = y * img.width + x;
        if ( img.channels == 1 ) {
          img.data[img.width*y+x] = row[x];
        } else if (img.channels == 3) {
          img.data[offset*3] = row[3*x];
          img.data[offset*3+1] = row[3*x+1];
          img.data[offset*3+2] = row[3*x+2];
        }

      }
    }
  }
  for (int y = 0; y < img.height; y++)
    free(rowPointers[y]);
  free(rowPointers);
  return img;
}
void writePNG(Image& img, std::string filename) {
  png_byte colorType;
  if(!img.isValid()){
    cout << "ERR: Image not valid" << endl;
    return;
  }
  if(img.bits != 8){
    cout << "writing  != 8 bit images currently not implemented" << endl;
    return;
  }
  if(img.channels == 1){
    colorType = PNG_COLOR_TYPE_GRAY;
  }
  else if(img.channels == 3){
    colorType = PNG_COLOR_TYPE_RGB;
  }
  png_byte bitDepth = img.bits;;

  png_structp pngPtr;
  png_infop infoPtr;
  png_bytep * rowPointers;

  //Copy Image matrix such that we can do a conservative resize to its correct size
  FILE *fp = fopen(filename.c_str(), "wb");
  if (!fp)
    cout << "ERR: File"  << filename << " could not be opened for writing" << endl;

  // init
  pngPtr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!pngPtr)
    cout << "ERR: png_create_write_struct failed" << endl;

  infoPtr = png_create_info_struct(pngPtr);
  if (!infoPtr)
    cout << "ERR: png_create_info_struct failed" << endl;

  if (setjmp(png_jmpbuf(pngPtr)))
    cout << "ERR: Error during init_io" << endl;

  png_init_io(pngPtr, fp);

  // write header
  if (setjmp(png_jmpbuf(pngPtr)))
    cout << "ERR: Error during writing header" << endl;

  png_set_IHDR(pngPtr, infoPtr, img.width, img.height,
      bitDepth, colorType, PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(pngPtr, infoPtr);

  // write bytes
  if (setjmp(png_jmpbuf(pngPtr)))
    cout << "ERR: Error during writing bytes" << endl;

  //Allocate row pointers
  rowPointers = (png_bytep*) malloc(sizeof(png_bytep) * img.height);
  for (int y = 0; y < img.height; y++)
    rowPointers[y] = (png_byte*) malloc(png_get_rowbytes(pngPtr, infoPtr));

  for (int y = 0; y < img.height; y++) {
    png_byte* row = rowPointers[y];
    for (int x = 0; x < img.width; x++) {
      int offset;
      if(img.channels == 1){
        offset = y * img.width + x;
        row[x * img.channels] = img.data[offset];
      }else if(img.channels ==  3){
        offset = (y * img.width + x) * 3;
        row[x * 3] = img.data[offset];
        row[x * 3 + 1] = img.data[offset+1];
        row[x * 3 + 2] = img.data[offset+2];
      }
    }
  }

  png_write_image(pngPtr, rowPointers);

  // end write
  if (setjmp(png_jmpbuf(pngPtr)))
    cout << "ERR: Error during end of write" << endl;

  png_write_end(pngPtr, NULL);

  for (int y = 0; y < img.height; y++)
    free(rowPointers[y]);
  free(rowPointers);

  fclose(fp);
}
#endif
