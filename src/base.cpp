/*
 * Created on: October 1, 2017
 * Author: Hjalmar Boulouednine <hjalmar@b-nine.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <functional>
#include <Rcpp.h>

// marathon includes
#include "marathon/binary_matrix/binary_matrix.h"
#include "marathon/binary_matrix/metrics.h"
#include "marathon/binary_matrix/random_generator.h"
#include "marathon/binary_matrix/sampling_engine.h"

// fixed margin
#include "marathon/binary_matrix/fixed_margin/random_generator_mcmc.h"
#include "marathon/binary_matrix/fixed_margin/random_generator_exact.h"
#include "marathon/binary_matrix/fixed_margin/instance.h"
#include "marathon/binary_matrix/fixed_margin/switch_chain.h"
#include "marathon/binary_matrix/fixed_margin/curveball.h"
#include "marathon/binary_matrix/fixed_margin/realize.h"

// interval margin
#include "marathon/binary_matrix/interval_margin/random_generator_mcmc.h"
#include "marathon/binary_matrix/interval_margin/random_generator_exact.h"
#include "marathon/binary_matrix/interval_margin/instance.h"
#include "marathon/binary_matrix/interval_margin/simple_chain_dynamic.h"
#include "marathon/binary_matrix/interval_margin/realize.h"


Rcpp::NumericMatrix toMatrix(const marathon::binary_matrix::BinaryMatrix& mtr)
{
  Rcpp::NumericMatrix X(mtr.getNumRows(), mtr.getNumCols());
  for (size_t i = 0; i < mtr.getNumRows(); ++i) {
    for (size_t j = 0; j < mtr.getNumCols(); ++j) {
      X(i, j) = mtr.get(i, j);
    }
  }
  return X;
}



std::unique_ptr<marathon::binary_matrix::RandomGenerator>
  construct_rg(const Rcpp::NumericVector &rowsums,
               const Rcpp::NumericVector &colsums,
               size_t steps, 
               const std::string& method)
  {
    using namespace marathon::binary_matrix;
    
    std::vector<int> v_rowsums(rowsums.begin(), rowsums.end());
    std::vector<int> v_colsums(colsums.begin(), colsums.end());
    
    fixed_margin::Instance inst(v_rowsums, v_colsums);
    
    if (method.compare("exact") == 0) {
      return std::make_unique<fixed_margin::RandomGeneratorExact>(inst);
    } else if (method.compare("ktv-switch") == 0) {
      return std::make_unique<fixed_margin::RandomGeneratorMCMC>(inst,
                                                                 fixed_margin::RandomGeneratorMCMC::classical_switch,
                                                                 steps);
    } else if (method.compare("edge-switch") == 0) {
      return std::make_unique<fixed_margin::RandomGeneratorMCMC>(inst,
                                                                 fixed_margin::RandomGeneratorMCMC::edge_switch,
                                                                 steps);
    } else if (method.compare("curveball") == 0) {
      return std::make_unique<fixed_margin::RandomGeneratorMCMC>(inst,
                                                                 fixed_margin::RandomGeneratorMCMC::curveball,
                                                                 steps);
    } else {
      std::cerr << "Unknown method: " << method << std::endl;
      return std::unique_ptr<interval_margin::RandomGeneratorMCMC>(nullptr);
    }
  }

std::unique_ptr<marathon::binary_matrix::RandomGenerator>
  construct_rg(const Rcpp::NumericVector &rowsumsL,
               const Rcpp::NumericVector &rowsumsU,
               const Rcpp::NumericVector &colsumsL,
               const Rcpp::NumericVector &colsumsU,
               size_t steps, 
               const std::string& method)
  {
    using namespace marathon::binary_matrix;
    
    std::vector<int> v_rowsumsL(rowsumsL.begin(), rowsumsL.end());
    std::vector<int> v_rowsumsU(rowsumsU.begin(), rowsumsU.end());
    std::vector<int> v_colsumsL(colsumsL.begin(), colsumsL.end());
    std::vector<int> v_colsumsU(colsumsU.begin(), colsumsU.end());
    
    interval_margin::Instance inst(v_rowsumsL, v_rowsumsU, v_colsumsL, v_colsumsU);
    
    if (method.compare("exact") == 0) {
      return std::make_unique<interval_margin::RandomGeneratorExact>(inst);
    } else if (method.compare("simple") == 0) {
      return std::make_unique<interval_margin::RandomGeneratorMCMC>(inst,
                                                                    interval_margin::RandomGeneratorMCMC::simple,
                                                                    steps);
    } else if (method.compare("informed") == 0) {
      return std::make_unique<interval_margin::RandomGeneratorMCMC>(inst,
                                                                    interval_margin::RandomGeneratorMCMC::informed,
                                                                    steps);
    } else {
      std::cerr << "Unknown method: " << method << std::endl;
      return std::unique_ptr<interval_margin::RandomGeneratorMCMC>(nullptr);
    }
  }

/**
 * Create N random binary matrices with prescribed row and column sums.
 */

// [[Rcpp::export]]
std::vector<Rcpp::NumericMatrix>
  CppSampleBinaryMatricesInterval(const Rcpp::NumericVector &rowsumsL,
                                  const Rcpp::NumericVector &rowsumsU,
                                  const Rcpp::NumericVector &colsumsL,
                                  const Rcpp::NumericVector &colsumsU,
                                  int N, int steps,
                                  const std::string& method)
  {
    auto rg = construct_rg(rowsumsL,
                           rowsumsU,
                           colsumsL,
                           colsumsU,
                           steps,
                           method);
    
    std::vector<Rcpp::NumericMatrix> res;
    marathon::binary_matrix::SamplingEngine se(*rg);
    
    auto samples = se.sample(N);
    
    for (size_t i = 0; i < N; i++) {
      res.emplace_back(toMatrix(samples[i]));
    }
    
    return res;
  }

/**
 * Create N random binary matrices with prescribed row and column sums.
 */

// [[Rcpp::export]]
std::vector<Rcpp::NumericMatrix>
  CppSampleBinaryMatricesFixed(const Rcpp::NumericVector &rowsums,
                               const Rcpp::NumericVector &colsums,
                               int N, int steps,
                               const std::string& method)
  {
    
    auto rg = construct_rg(rowsums,
                           colsums,
                           steps,
                           method);
    
    std::vector<Rcpp::NumericMatrix> res;
    marathon::binary_matrix::SamplingEngine se(*rg);
    
    auto samples = se.sample(N);
    
    for (size_t i = 0; i < N; i++) {
      res.emplace_back(toMatrix(samples[i]));
    }
    
    return res;
  }

// [[Rcpp::export]]
bool CppIsRealizableInterval(const Rcpp::NumericVector &rowsumsL,
                             const Rcpp::NumericVector &rowsumsU,
                             const Rcpp::NumericVector &colsumsL,
                             const Rcpp::NumericVector &colsumsU)
{
  using namespace marathon::binary_matrix;
  
  std::vector<int> v_rowsumsL(rowsumsL.begin(), rowsumsL.end());
  std::vector<int> v_rowsumsU(rowsumsU.begin(), rowsumsU.end());
  std::vector<int> v_colsumsL(colsumsL.begin(), colsumsL.end());
  std::vector<int> v_colsumsU(colsumsU.begin(), colsumsU.end());
  
  interval_margin::Instance inst(v_rowsumsL, v_rowsumsU, v_colsumsL, v_colsumsU);
  return interval_margin::isRealizable(inst);
}

// [[Rcpp::export]]
bool CppIsRealizableFixed(const Rcpp::NumericVector &rowsums,
                          const Rcpp::NumericVector &colsums)
{
  using namespace marathon::binary_matrix;
  
  std::vector<int> v_rowsums(rowsums.begin(), rowsums.end());
  std::vector<int> v_colsums(colsums.begin(), colsums.end());
  
  fixed_margin::Instance inst(v_rowsums, v_colsums);
  return fixed_margin::isRealizable(inst);
}
