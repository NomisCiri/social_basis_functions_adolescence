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
struct ddm2w : public Worker {
  
  // input vectors/matrices to read from
  const double d_v;
  const double d_h;
  const double d_b;
  const double thres;
  const double nDT;
  const double tIn;
  const double bias;
  const double vd;
  const double hd;
  const double bd;
  GEN gen;
  
  // output vector to write to
  RVector<double> vecOut;
  
  // initialize from Rcpp input and output matrices/vectors (the RMatrix/RVector class
  // can be automatically converted to from the Rcpp matrix/vector type)

  ddm2w(const double d_v, const double d_h, const double d_b, const double thres, const double nDT, const double tIn, const double bias, const double vd, const double hd, const double bd, NumericVector vecOut , GEN gen)
    : d_v(d_v), d_h(d_h), d_b(d_b), thres(thres), nDT(nDT), tIn(tIn), bias(bias), vd(vd), hd(hd), bd(bd), gen(gen), vecOut(vecOut) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    
    double T = 5.2, dt = 0.008, lt;
    lt = (int)(T/dt);//time divided by descrete time.
    
    // define vectors with info when to consider attribute or not to consider it (0)
    std::vector<double> vec_tsecBF(lt,1);
    std::vector<double> vec_tprimBF(lt,1);
    std::vector<double> vec_tBonus(lt,1);
    int aux = std::abs((int)(tIn/dt));
    
    // here RST (= tIN) is taken into account
    // if tIN > 0 -> primBF enters first (secBF set to 0)
    // if tIN < 0 -> secBF enters first (primBF set to 0)
    // bonus always enters with the first attribute
    if (tIn > 0) {
      for (int t=0; t<aux; t++) {
        vec_tsecBF[t] = 0;
      }
    }
    else if (tIn < 0) {
      for (int t=0; t<aux; t++) {
        vec_tprimBF[t] = 0;
      }
    }

    for (std::size_t i = begin; i < end; i++) {
      vecOut[i] = T;
      //std::vector<double> X(lt,)
      double X = bias;
      int flag = 0;
      double cont = 0;
      double noise = 0;
      
      while (flag==0 && cont<lt) {
        
        noise=gen()*sqrt(dt);
        X = X + (d_v*vd*vec_tprimBF[cont] + d_h*hd*vec_tsecBF[cont] + d_b*bd*vec_tBonus[cont])*dt + noise;
        
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
NumericVector ddm2_parallel(double d_v, double d_h, double d_b, double thres, double nDT, double tIn, double bias, double vd, double hd, double bd, unsigned int N) {
  
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
  ddm2w ddm2w(d_v, d_h, d_b, thres, nDT, tIn, bias, vd, hd, bd, vecOut, gen);
  
  // call the worker
  parallelFor(0, N, ddm2w);
  
  return vecOut;
}







// [[Rcpp::export]]
NumericVector ddm2_serial(double d_v, double d_h, double d_b, double thres, double nDT, double tIn, double bias, double vd, double hd, double bd, unsigned int N) {
  
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
    
    // Define vectors with info when to consider attribute or not to consider it (0)
    std::vector<double> vec_tsecBF(lt,1);
    std::vector<double> vec_tprimBF(lt,1);
    std::vector<double> vec_tBonus(lt,1);
    int aux = std::abs((int)(tIn/dt));
    
    // Here RST (= tIN) is taken into account
    // If tIN > 0 -> primBF enters first (secBF set to 0)
    // If tIN < 0 -> secBF enters first (primBF set to 0)
    // Bonus always enters with the first attribute
    if (tIn > 0) {
      for (int t=0; t<aux; t++) {
        vec_tsecBF[t] = 0;
      }
    }
    else if (tIn < 0) {
      for (int t=0; t<aux; t++) {
        vec_tprimBF[t] = 0;
      }
    }
    
    vecOut[i] = T;
    double X = bias;
    int flag = 0;
    double cont = 0;
    double noise = 0;
  
    while (flag==0 && cont<lt) {
      noise=gen()*sqrt(dt);
      X = X + (d_v*vd*vec_tprimBF[cont] + d_h*hd*vec_tsecBF[cont] + d_b*bd*vec_tBonus[cont])*dt + noise;
  
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

