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
  QUESO::GslVector thetaVector(paramSpace.zeroVector());
  thetaVector[0] = -5.12;
  thetaVector[1] = 5.12;

  likelihoodRoutine_DataType likelihoodRoutine_Data;
  likelihoodRoutine_Data.thetaVector = &thetaVector;


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
  double scaleFac = 1;

  QUESO::GslVector paramInitials(paramSpace.zeroVector());
  paramInitials[0] = 0.1;
  paramInitials[1] = 1.4;

  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());
  proposalCovMatrix(0,0) = scaleFac*8.; proposalCovMatrix(0,1) = scaleFac*4.;
  proposalCovMatrix(1,0) = scaleFac*4.; proposalCovMatrix(1,1) = scaleFac*16.;

  ip.solveWithBayesMetropolisHastings(NULL,paramInitials, &proposalCovMatrix);

  return;
}
