//
//  DoubleDip.java
//  FractalFunctions
//
//  Created by Cara MacNish on 25/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.
//
//  Translated in C++ by Giuseppe A. Trunfio on Nov 10, 2013

/**
 * DoubleDip is a base function that uses a segment of a sextic polynomial. It is similar in shape
 * to DoubleCosine, but considerably faster. It results in a continuously differentiable surface.
 * @see UnitFunction1D
 * @author {@link <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>}, University of Western Australia
 * @version 1.0RC1, 7th Nov 2007
 * <br>For the latest version and additional information see the
 * {@link <a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>}
 */

#include "DoubleDip.h"
#include <math.h>       /* fmod */

DoubleDip::DoubleDip() : UnitFunction1D()
{
}

DoubleDip::DoubleDip(double centre, double scale) : UnitFunction1D(centre, scale)
{
}

double DoubleDip::getValue(double point)
{
    double depth = 0;
    double xs;
    double x = (point - centre) / scale;
    if (x > -0.5 && x < 0.5)
    {
        xs = 4*x*x; // scale to -1 to 1 -> (2x)^2
        depth = (-96*xs*xs*xs + 193*xs*xs - 98*xs + 1) * scale; // much quicker to use * than Math.pow;
    }
    return depth;
}

double DoubleDip::twist(double x, double y)
{
    double dx = 0;
    y = fmod(y, 1);

    double ys = y*y;
    if (y > 0) // twisting quartic
    {
        dx = 4*(ys*ys - 2*ys*y + ys);
    }
    else // twisting quartic
    {
        dx = 4*(ys*ys + 2*ys*y + ys);
    }
    return dx;
}
