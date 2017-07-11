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

    
    double returnValue = 0.;
    double prod = 0.;
    prod = cos(paramValues[0]/2) * cos(paramValues[1]/2) * cos(paramValues[2]/2) * cos(paramValues[3]/2); //* cos(paramValues[4]/2);

    double logL = pow(( 2 + prod ), 5);
    //logL = log(logL);

    double resultValue = -2*logL;

    if (resultValue == INFINITY) {
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
      resultValue = 1040.;
    }
	//if (resultValue == -INFINITY) {
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
      //resultValue = -1040.;
    //}	

    returnValue = -.5*resultValue;

  
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
