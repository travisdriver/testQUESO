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

   
    const QUESO::GslVector& meanVector1 =
      *((likelihoodRoutine_DataType *) functionDataPtr)->meanVector1;
    const QUESO::GslVector& meanVector2 =
      *((likelihoodRoutine_DataType *) functionDataPtr)->meanVector2;
    const QUESO::GslVector& meanVector3 =
      *((likelihoodRoutine_DataType *) functionDataPtr)->meanVector3;
    const QUESO::GslVector& meanVector4 =
      *((likelihoodRoutine_DataType *) functionDataPtr)->meanVector4;
                                                                          
    double r = 1;
    double w = 0.1;
    double A = .9;
                                                                          
    QUESO::GslVector diffVec1(paramValues - meanVector1);
    double logLike1 = A - pow((diffVec1.norm2() - r),2) / (2 * pow(w,2));
    logLike1 = exp(logLike1);
                                                                          
    QUESO::GslVector diffVec2(paramValues - meanVector2);
    double logLike2 = A - pow((diffVec2.norm2() - r),2) / (2 * pow(w,2));
    logLike2 = exp(logLike2);  	

    QUESO::GslVector diffVec3(paramValues - meanVector3);
    double logLike3 = A - pow((diffVec3.norm2() - r),2) / (2 * pow(w,2));
    logLike3 = exp(logLike3);

    QUESO::GslVector diffVec4(paramValues - meanVector4);
    double logLike4 = A - pow((diffVec4.norm2() - r),2) / (2 * pow(w,2));
    logLike4 = exp(logLike4);
                                                                          
    double logL = log(logLike1 + logLike2 + logLike3 + logLike4);
	double resultValue = logL;

    if (resultValue == -INFINITY) {
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
      resultValue = -1040.;
    }
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
