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

#include <bimodal_likelihood.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

int likelihoodCounter = 0;

double likelihoodRoutine(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  const void*             functionDataPtr,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect)
{
  likelihoodCounter++;

  if (paramDirection  ||
      functionDataPtr ||
      gradVector      ||
      hessianMatrix   ||
      hessianEffect) {}; // just to remove compiler warning

    
  const QUESO::GslVector& meanVecG =
  	  *((likelihoodRoutine_DataType *) functionDataPtr)->meanVecG;
  const QUESO::GslVector& meanVecR =
  	  *((likelihoodRoutine_DataType *) functionDataPtr)->meanVecR;
  const QUESO::GslMatrix& covMatG  =
  	  *((likelihoodRoutine_DataType *) functionDataPtr)->covMatG;
  const QUESO::GslMatrix& covMatR  =
  	  *((likelihoodRoutine_DataType *) functionDataPtr)->covMatR;
  const QUESO::GslMatrix& invCovMatG  =
  	  *((likelihoodRoutine_DataType *) functionDataPtr)->invCovMatG;
  const QUESO::GslMatrix& invCovMatR  =
  	  *((likelihoodRoutine_DataType *) functionDataPtr)->invCovMatR;
                                                                                          
  double sqrtDetInvCovMatR = sqrt(4.93827160493827);
  double sqrtDetInvCovMatG = sqrt(44.4444444444445);
  double invsqrt2pi = 0.398942280401432;
  double nd = 2;
  double normFac = 0.6;
  double a = 0.7; double b = 1.5;
                                                                                          
  QUESO::GslVector diffVecG(paramValues - meanVecG);
  double gaussian = pow(invsqrt2pi,nd) * sqrtDetInvCovMatG \
  				             * exp( -0.5 * scalarProduct( diffVecG , invCovMatG.multiply(diffVecG) ) );
  
  QUESO::GslVector ParamR = paramValues;
  ParamR[1] = a * ParamR[1];
  ParamR[0] = ParamR[0]/a - b * (pow(ParamR[1], 2) + pow(a, 2));
  QUESO::GslVector diffVecR(ParamR - meanVecR);
  double rosenbrock = pow(invsqrt2pi, nd) * sqrtDetInvCovMatG \
  					             * exp( -0.5 * scalarProduct( diffVecR , invCovMatR.multiply(diffVecR) ) );
  
  double getLogProbGR = rosenbrock + normFac * gaussian;
  if (getLogProbGR==0.)
  {
  	  getLogProbGR = -std::numeric_limits<double>::max();
  }
  else
  {
  	  getLogProbGR = log( rosenbrock + normFac * gaussian );
  }

  double resultValue = getLogProbGR;
  //resultValue = log(resultValue);

  //if (resultValue == INFINITY) {
    //std::cerr << "WARNING In likelihoodRoutine"
    //          << ", fullRank "       << paramValues.env().fullRank()
    //          << ", subEnvironment " << paramValues.env().subId()
    //          << ", subRank "        << paramValues.env().subRank()
    //          << ", inter0Rank "     << paramValues.env().inter0Rank()
    //          << ": x = "            << x
    //          << ", z1 = "           << z1
    //          << ", z2 = "           << z2
    //          << ", resultValue = "  << resultValue
    //          << std::endl;
    //resultValue = 1040.;
  //}
  //if (resultValue == -INFINITY){
  //	resultValue = -1040.;
  //}


  double returnValue = resultValue;

  
  if (paramValues.env().exceptionalCircumstance()) {
    if ((paramValues.env().subDisplayFile()       ) &&
        (paramValues.env().displayVerbosity() > 0)) { // detailed output debug
      *paramValues.env().subDisplayFile() << "Leaving likelihood function"
                                          << ": paramValues = " << paramValues
                                          << ", returnValue = " << returnValue
                                          << std::endl;
    }
  }
  
  return returnValue;
}
