//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008-2015 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <example_compute.h>
#include <example_likelihood.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>

void compute(const QUESO::FullEnvironment& env) {
  // Step 1 of 5: Instantiate the parameter space
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix>
    paramSpace(env, "param_", 2, NULL);

  // Step 2 of 5: Instantiate the parameter domain
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  paramMins.cwSet(-INFINITY);
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMaxs.cwSet( INFINITY);
  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_",paramSpace,paramMins,paramMaxs);

  // Step 3 of 5: Instantiate the likelihood function object
  QUESO::GslVector meanVecG(paramSpace.zeroVector());
  meanVecG[0] = 3.5;
  meanVecG[1] =  0;

  QUESO::GslVector meanVecR(paramSpace.zeroVector());
  meanVecR[0] = -5.;
  meanVecR[1] =  0;

  QUESO::GslMatrix covMatG(paramSpace.zeroVector());
  covMatG(0,0) = 0.15; covMatG(0,1) = 0.;
  covMatG(1,0) = 0.; covMatG(1,1) = 0.15;

  QUESO::GslMatrix covMatR(paramSpace.zeroVector());
  covMatR(0,0) = 0.25; covMatR(0,1) = 0.;
  covMatR(1,0) = 0.; covMatR(1,1) = 0.81;

  QUESO::GslMatrix invCovMatG(paramSpace.zeroVector());
  invCovMatG(0,0) = 6.66666666666667; invCovMatG(0,1) = 0.;
  invCovMatG(1,0) = 0.; invCovMatG(1,1) = 6.66666666666667;

  QUESO::GslMatrix invCovMatR(paramSpace.zeroVector());
  invCovMatR(0,0) = 4.; invCovMatR(0,1) = 0.;
  invCovMatR(1,0) = 0.; invCovMatR(1,1) = 1.23456790123457;

  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.meanVecG = &meanVecG;
  likelihoodRoutine_Data.covMatG  = &covMatG;
  likelihoodRoutine_Data.meanVecR = &meanVecR;
  likelihoodRoutine_Data.covMatR  = &covMatR;
  likelihoodRoutine_Data.invCovMatG = &invCovMatG;
  likelihoodRoutine_Data.invCovMatR  = &invCovMatR;

  QUESO::GenericScalarFunction<QUESO::GslVector,QUESO::GslMatrix>
    likelihoodFunctionObj("like_",
                          paramDomain,
                          likelihoodRoutine,
                          (void *) &likelihoodRoutine_Data,
                          true); // routine computes [ln(function)]

  // Step 4 of 5: Instantiate the inverse problem
  QUESO::UniformVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    priorRv("prior_", paramDomain);
  QUESO::GenericVectorRV<QUESO::GslVector,QUESO::GslMatrix>
    postRv("post_", paramSpace);
  QUESO::StatisticalInverseProblem<QUESO::GslVector,QUESO::GslMatrix>
    ip("", NULL, priorRv, likelihoodFunctionObj, postRv);

  // Step 5 of 5: Solve the inverse problem
  double scal = .1;
  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.1;
  paramInitials[1] = -1.4;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = scal*8.; proposalCovMatrix(0,1) = scal*4.;
  proposalCovMatrix(1,0) = scal*4.; proposalCovMatrix(1,1) = scal*16.;

  ip.solveWithBayesMLSampling();

  return;
}
