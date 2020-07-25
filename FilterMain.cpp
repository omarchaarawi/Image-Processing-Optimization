#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "Filter.h"
#include "omp.h"

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

class Filter *
readFilter(string filename)
{
    ifstream input(filename.c_str());

    if ( ! input.bad() ) {
        int size = 0;
        input >> size;
        Filter *filter = new Filter(size);
        int div;
        input >> div;
        filter -> setDivisor(div);
        for (int i=0; i < size; i++) {
          for (int j=0; j < size; j++) {
        int value;
        input >> value;
        filter -> set(i,j,value);
            }
        }
        return filter;
    } else {
        cerr << "Bad input in readFilter:" << filename << endl;
        exit(-1);
    }
}


double applyFilter(class Filter *filter, cs1300bmp *input, cs1300bmp *output)
{    
    long long cycStart, cycStop;

    cycStart = rdtscll();

    output -> width = input -> width;
    output -> height = input -> height;
    
    // define height_minus_one and height_minus_one for use later in loop
    const short int height_minus_one = (input -> height) - 1;
    const short int width_minus_one = (input -> width) - 1;
    
   
    // get filter data outside of loop so we don't have to keep calling filter->get(i,j)
//     int local_filter_data[3][3];
//     short int i;
//    #pragma omp parallel for num_threads(3)
//     for(i = 0; i < 3; ++i){
//         local_filter_data[i][0] = filter -> get(i,0);
//         local_filter_data[i][1] = filter -> get(i,1);
//         local_filter_data[i][2] = filter -> get(i,2);
//     }
    short int a,b,c,d,e,f,g,h,i;
    a = filter->get(0,0);
    b = filter->get(0,1);
    c = filter->get(0,2);
    d = filter->get(1,0);
    e = filter->get(1,1);
    f = filter->get(1,2);
    g = filter->get(2,0);
    h = filter->get(2,1);
    i = filter->get(2,2);

    
    //int filter_size = filter -> getSize();   
    // called filter->getSize() rather than in loop to avoid unnecessary function call
    const short int filter_divisor = filter -> getDivisor();  // called filter->getDivisor() outside of loop to avoid unnecessary function call
    

    // run loops in parallel using omp
    // my source for this was https://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-loop.html#Loopschedules 
    #pragma omp parallel for 
    for(int plane = 0; plane < 3; ++plane) {
        for(int row = 1; row < height_minus_one; row++) {    // swapped column, row, and plane to create better memory efficiency
            
            // save ptr to row plus one and row minus one for use in inner loop
            // should improve spatial locality
            const short int* row_minus_one_ptr = &input->color[plane][row - 1][0];
            const short int* row_plus_one_ptr = &input->color[plane][row + 1][0];
            for(int col = 1; col < width_minus_one ; col+=2) {   // go two columns at a time                 
                // split outputs in two
                output -> color[plane][row][col] = 0;
                output -> color[plane][row][col+1] = 0;
                int output_temp = 0;
                int output_temp2 = 0;
                
           int x0, x1, x2, y0, y1, y2, z0, z1, z2;      // split into three processes
                int xx0, xx1, xx2, yy0, yy1, yy2, zz0, zz1, zz2;  // unrolled once

                x0 = (row_minus_one_ptr[col - 1] * a);
                xx0 = (row_minus_one_ptr[col] * a);

                x1 = (row_minus_one_ptr[col] * b);
                xx1 = (row_minus_one_ptr[col+1] * b);

                x2 = (row_minus_one_ptr[col + 1] * c);
                xx2 = (row_minus_one_ptr[col + 2] * c);

                y0 = (input -> color[plane][row][col - 1] * d);
                yy0 = (input -> color[plane][row][col] * d);

                y1 = (input -> color[plane][row][col] * e);
                yy1 = (input -> color[plane][row][col+1] * e);

                y2 = (input -> color[plane][row][col + 1] * f);
                yy2 = (input -> color[plane][row][col + 2] * f);

                z0 = (row_plus_one_ptr[col - 1] * g);
                zz0 = (row_plus_one_ptr[col] * g);

                z1 = (row_plus_one_ptr[col] * h);
                zz1 = (row_plus_one_ptr[col+1] * h);
                
                z2 = (row_plus_one_ptr[col + 1] * i);
                zz2 = (row_plus_one_ptr[col + 2] * i);
                
                output_temp += (x0 + x1 + x2 + y0 + y1 + y2 + z0 + z1 + z2);
                output_temp2 += (xx0 + xx1 + xx2 + yy0 + yy1 + yy2 + zz0 + zz1 + zz2);
        


            // only divide if filter divisor is not 1 or else it's a waste of computing power    
            if(filter_divisor != 1){
                output_temp = output_temp / filter_divisor;
                output_temp2 = output_temp2 / filter_divisor;
            }
//             output_temp = (filter_divisor == 1) ? output_temp : (output_temp/filter_divisor);
//             output_temp2 = (filter_divisor == 1) ? output_temp2 : (output_temp2/filter_divisor);
                
            // replaced if statements with faster conditionals
            output_temp = (output_temp < 0) ? 0 : output_temp;
            output_temp2 = (output_temp2 < 0) ? 0 : output_temp2;


            output_temp = (output_temp > 255) ? 255 : output_temp;
            output_temp2 = (output_temp2 > 255) ? 255 : output_temp2;


            output -> color[plane][row][col] = output_temp;
            output -> color[plane][row][col+1] = output_temp2;

          }
        }
    }
    
    cycStop = rdtscll();
    double diff = cycStop - cycStart;
    double diffPerPixel = diff / (output -> width * output -> height);
    fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
      diff, diff / (output -> width * output -> height));
    return diffPerPixel;
}