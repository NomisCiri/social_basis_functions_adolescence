#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <time.h>

using namespace Rcpp;
using namespace RcppParallel;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator

// [[Rcpp::depends(RcppParallel)]]
struct vanilla_ddm2w : public Worker {
  
  // input vectors/matrices to read from
  const double d_v;
  const double thres;
  const double nDT;
  const double bias;
  const double vd;
  GEN gen;
  
  // output vector to write to
  RVector<double> vecOut;
  
  // initialize from Rcpp input and output matrices/vectors (the RMatrix/RVector class
  // can be automatically converted to from the Rcpp matrix/vector type)
  
  vanilla_ddm2w(const double d_v, const double thres, const double nDT, const double bias, const double vd, NumericVector vecOut , GEN gen)
    : d_v(d_v), thres(thres), nDT(nDT), bias(bias), vd(vd), gen(gen), vecOut(vecOut) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    
    double T = 5.2, dt = 0.008, lt;
    lt = (int)(T/dt);//time divided by discrete time.
    
    for (std::size_t i = begin; i < end; i++) {
      vecOut[i] = T;
      double X = bias;
      int flag = 0;
      double cont = 0;
      double noise = 0;
      
      while (flag==0 && cont<lt) {
        
        noise=gen()*sqrt(dt);
        X = X + (d_v*vd)*dt + noise;
        
        if (X > 1) {
          flag=1;
          vecOut[i] = nDT + cont*dt;
        }
        else if (X < -1) {
          flag=1;
          vecOut[i] = -nDT -cont*dt;
        }
        cont++;
        
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector vanilla_ddm2_parallel(double d_v, double thres, double nDT, double bias, double vd, unsigned int N) {
  
  const double sd_n = 1.4;
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0,sd_n);
  GEN  gen(eng,dist);
  
  //output vector
  NumericVector vecOut(N);
  
  // create the worker
  vanilla_ddm2w vanilla_ddm2w(d_v, thres, nDT, bias, vd, vecOut, gen);
  
  // call the worker
  parallelFor(0, N, vanilla_ddm2w);
  
  return vecOut;
}

// [[Rcpp::export]]
NumericVector vanilla_ddm2_serial(double d_v, double thres, double nDT, double bias, double vd, unsigned int N) {
  
  const double sd_n = 1.4;
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0,sd_n);
  GEN  gen(eng,dist);
  
  // Output vector
  NumericVector vecOut(N);
  
  // Execute the logic for each index serially
  for (unsigned int i = 0; i < N; ++i) {
    double T = 5.2, dt = 0.008, lt;
    lt = (int)(T/dt);
    
    vecOut[i] = T;
    double X = bias;
    int flag = 0;
    double cont = 0;
    double noise = 0;
    
    while (flag==0 && cont<lt) {
      noise=gen()*sqrt(dt);
      X = X + (d_v*vd)*dt + noise;
      
      if (X > 1) {
        flag=1;
        vecOut[i] = nDT + cont*dt;
      } 
      else if (X < -1) {
        flag=1;
        vecOut[i] = -nDT -cont*dt;
      } 
      cont++;
    } 
  }
  
  return vecOut;
} 