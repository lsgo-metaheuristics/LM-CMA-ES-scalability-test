#pragma once

#include "UnitFunction1D.h"

//
//  DoubleDip.java
//  FractalFunctions
//
//  Created by Cara MacNish on 25/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.

/// <summary>
/// DoubleDip is a base function that uses a segment of a sextic polynomial. It is similar in shape
/// to DoubleCosine, but considerably faster. It results in a continuously differentiable surface. </summary>
/// <seealso cref= UnitFunction1D
/// @author <seealso cref="<a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>"/>, University of Western Australia
/// @version 1.0RC1, 7th Nov 2007
/// <br>For the latest version and additional information see the
/// <seealso cref="<a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>"/> </seealso>
class DoubleDip : public UnitFunction1D
{

public:
    DoubleDip();

    DoubleDip(double centre, double scale);

    virtual double getValue(double point) override;

    virtual double twist(double x, double y) override;

};
