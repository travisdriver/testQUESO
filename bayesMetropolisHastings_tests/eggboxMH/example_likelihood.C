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

#include <example_likelihood.h>
#include <cmath>
#include <limits>

static unsigned int likelihoodCounter = 0;

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{

  likelihoodCounter ++;
  // Logic just to avoid warnings from INTEL compiler
  const QUESO::GslVector* aux1 = paramDirection;
  if (aux1) {};
  aux1 = gradVector;
  aux1 = hessianEffect;
  QUESO::GslMatrix* aux2 = hessianMatrix;
  if (aux2) {};

  // Just checking: the user, at the application level, expects
  // vector 'paramValues' to have size 2.
  UQ_FATAL_TEST_MACRO(paramValues.sizeGlobal() != 2,
                      QUESO::UQ_UNAVAILABLE_RANK,
                      "likelihoodRoutine()",
                      "paramValues vector does not have size 2");

  // Actual code
  //
  // This code exemplifies multiple Metropolis-Hastings solvers, each calling this likelihood
  // routine. In this simple example, only node 0 in each subenvironment does the job even
  // though there might be more than one node per sub-environment. In a more realistic
  // situation, if the user is asking for multiple nodes per subenvironment, then the model
  // code in the likelihood routines might really demand more than one node. Here we use
  // 'env.subRank()' only. A realistic application might want to use either 'env.subComm()'
  // or 'env.subComm().Comm()'.

  double result = 0.;

  const QUESO::BaseEnvironment& env = paramValues.env();
  if (env.subRank() == 0) {
    const QUESO::GslVector& thetaVector =
      *((likelihoodRoutine_DataType *) functionDataPtr)->thetaVector;

    double x = paramValues[0];
    double y = paramValues[1];

    double logL = pow(( 2 + cos(x/2) * cos(y/2) ), 5);

    if (logL == 0.)
    {
      logL = -std::numeric_limits<double>::max();
    }
    else
    {
      // Do nothing
    }

    result = logL;
  }
  else {
    // Do nothing
  }

  return -0.5 * result;
}
