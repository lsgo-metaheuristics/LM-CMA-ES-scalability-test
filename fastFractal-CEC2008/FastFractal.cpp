//
//  FastFractal.java
//  FractalFunctions
//
//  Created by Cara MacNish on 19/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.
//
//  Translated in C++ by Giuseppe A. Trunfio on Nov 10, 2013

/**
 * This is the top level class (called by the user) for generating fast multidimensional
 * fractal functions.
 *<p>
 * Each function (or "landscape") is made from a base function (or <i>unit function</i>) at various scales.
 * For the fast multidimensional functions these are subclasses of {@link UnitFunction1D}.
 *<p>
 * The choice of base function determines a <i>class</i> of fractal functions.
 * Each instance within that class is then determined by three parameters:
 * the fractal depth, the density, and a sequence number.
 *<p>
 * For the motivation behind the fractal functions see:
 * MacNish, C., Towards Unbiased Benchmarking of Evolutionary and Hybrid Algorithms for Real-valued
 * Optimisation</a>, <i><a href="http://www.tandf.co.uk/journals/titles/09540091.asp">Connection Science</a></i>,
 * Vol. 19, No. 4, December 2007. Or visit <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>'s website.
 * @see UnitFunction1D
 * @author {@link <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>}, University of Western Australia
 * @version 1.0RC1, 7th Nov 2007
 * <br>For the latest version and additional information see the
 * {@link <a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>}
 */


#include "FastFractal.h"



FastFractal::FastFractal(UnitFunction1D* uf, int fractalDepth, int density, long index, int d)
{
    dimensions = d;
    unitFunction = uf;
    ff = new FractalFunction1D(unitFunction, fractalDepth, density, index);
}

FastFractal::~FastFractal()
{
    delete ff;
}

double FastFractal::evaluate(double *point)
{

    double depth = 0;
    double x, lastx, dx;
    ff->setIndex((6*dimensions - 1) + 1);
    lastx = point[dimensions - 1];
    for (int i = 0; i < dimensions; i++)
    {
        ff->setIndex(6*i + 1); // spread to small "prime-ish" seeds for diversity
        x = point[i];
        dx = unitFunction->twist(x,lastx);
        depth = depth + ff->evaluate(x + dx); // "twist" and evaluate
        lastx = x;
    }
    return depth;
}
