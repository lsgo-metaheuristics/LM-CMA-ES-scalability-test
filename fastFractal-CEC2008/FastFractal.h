#pragma once

#include "UnitFunction1D.h"
#include "FractalFunction1D.h"
#include <string>
#include <stdexcept>

//
//  FastFractal.java
//  FractalFunctions
//
//  Created by Cara MacNish on 19/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.

/// <summary>
/// This is the top level class (called by the user) for generating fast multidimensional
/// fractal functions.
/// <para>
/// Each function (or "landscape") is made from a base function (or <i>unit function</i>) at various scales.
/// For the fast multidimensional functions these are subclasses of <seealso cref="UnitFunction1D"/>.
/// </para>
/// <para>
/// The choice of base function determines a <i>class</i> of fractal functions.
/// Each instance within that class is then determined by three parameters:
/// the fractal depth, the density, and a sequence number.
/// </para>
/// <para>
/// For the motivation behind the fractal functions see:
/// MacNish, C., Towards Unbiased Benchmarking of Evolutionary and Hybrid Algorithms for Real-valued
/// Optimisation</a>, <i><a href="http://www.tandf.co.uk/journals/titles/09540091.asp">Connection Science</a></i>,
/// Vol. 19, No. 4, December 2007. Or visit <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>'s website.
/// </para>
/// </summary>
/// <seealso cref= UnitFunction1D
/// @author <seealso cref="<a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>"/>, University of Western Australia
/// @version 1.0RC1, 7th Nov 2007
/// <br>For the latest version and additional information see the
/// <seealso cref="<a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>"/> </seealso>
class FastFractal
{

private:
    UnitFunction1D *unitFunction;
    int dimensions;
    FractalFunction1D *ff;

    /// <summary>
    /// Create a fast fractal function generator. </summary>
    /// <param name="unitFunctionName"> the name of the base function for the generator. It must match the
    /// class name of a subclass of <seealso cref="UnitFunction1D"/>. </param>
    /// <param name="fractalDepth"> recursive depth of fractal - each increment adds detail at half the scale
    /// (double the resolution).
    /// Must be between 1 and 48 (the maximum supported by IEEE 64-bit floating point resolution). </param>
    /// <param name="density"> average number of base functions per unit area at each resolution </param>
    /// <param name="index"> the sequence number of this surface (for the given fractal depth and density) </param>
    /// <param name="dimensions"> number of dimensions (free variables) of the parameter space </param>
public:
    FastFractal(UnitFunction1D *unitFunction, int fractalDepth, int density, long index, int dimensions);
    ~FastFractal();

/// <summary>
/// Evaluate the function at the given point. </summary>
/// <param name="point"> the point to evaluate. The size of the array must match the dimension of the problem.
/// point[0] is co-ordinate x1, point[1] is co-ordinate x2, ..., point [D-1] is co-ordinate xD, where
/// D is the dimension. </param>
/// <returns> the value </returns>
    virtual double evaluate(double *point);
};

