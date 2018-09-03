//
//  UnitFunction1D.java
//  FractalFunctions
//
//  Created by Cara MacNish on 18/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.
//
//  Translated in C++ by Giuseppe A. Trunfio on Nov 10, 2013

/**
 * Each fractal function is constructed from base functions or "unit functions" at various scales
 * (a little like wavelets) chosen from appropriate probability distributions to preserve
 * self-similarity.
 *<p>
 * The name <i>unit function</i> comes from the fact that the non-zero portion of each function
 * is constrained to a unit hypercube, centred at the origin. In the 1-d case this is therefore the
 * interval (-0.5, 0.5).
 *<p>
 * This is the superclass of all 1-d unit functions. If you are writing your own unit
 * function you only need to call the constructors from your own constructors using
 * super(), and provide an implementation for getValue(point).
 *
 * @author {@link <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>}, University of Western Australia
 * @version 1.0RC1, 7th Nov 2007
 * <br>For the latest version and additional information see the
 * {@link <a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>}
 */


#include "UnitFunction1D.h"

UnitFunction1D::UnitFunction1D()
{
    InitializeInstanceFields();
}

UnitFunction1D::UnitFunction1D(double centre, double scale)
{
    InitializeInstanceFields();
    setCentre(centre);
    setScale(scale);
}

void UnitFunction1D::setParams(double centre, double scale)
{
    this->centre = centre;
    this->scale = scale;
}

void UnitFunction1D::setCentre(double centre)
{
    this->centre = centre;
}

double UnitFunction1D::getCentre()
{
    return centre;
}

void UnitFunction1D::setScale(double scale)
{
    this->scale = scale;
}

double UnitFunction1D::getScale()
{
    return scale;
}

double UnitFunction1D::twist(double x, double y)
{
    return 0;
}


void UnitFunction1D::InitializeInstanceFields()
{
    centre = 0;
    scale = 1;
}
