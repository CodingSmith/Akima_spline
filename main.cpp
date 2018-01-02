#include <iostream>
#include "akima_spline.h"

int main(int argc, char** argv) 
{
  std::vector<float> X(5), Y(5), Xi(4), Yi(4);
  X={1, 2, 3, 4, 5};
  Y={0.1, 0.7, 0.6, 1.1, 0.9};
  Xi={1.5, 2.5, 3.5, 4.5};
//  Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;
  
  cs::AkimaSpline spline_generator;
  spline_generator.set_points(X,Y);

  spline_generator.generate_points(Xi, Yi);

  for(size_t i = 0; i < Xi.size(); i++) {
    std::cout << Xi[i] <<" "<< Yi[i] << std::endl;
  }

}
