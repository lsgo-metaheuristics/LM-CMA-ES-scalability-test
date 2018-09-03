#pragma once

#include "RanTable.h"
#include "UnitFunction1D.h"
#include <cmath>

//
//  FractalFunction1D.java
//  FractalFunctions
//
//  Created by Cara MacNish on 18/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.

/// <summary>
/// This class is called by <seealso cref="FastFractal"/> to evaluate the components for each dimension to
/// evaluate fast very high-dimensional fractal functions (landscapes).
/// <para>
/// It can also be used directly to generate 1-dimensional fractal landscapes.
/// @author <seealso cref="<a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>"/>, University of Western Australia
/// @version 1.0RC1, 7th Nov 2007
/// <br>For the latest version and additional information see the
/// <seealso cref="<a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>"/>
/// </para>
/// </summary>
class FractalFunction1D
{

    // random tables
public:
    static const int DOUBLE_TABLE_SIZE = 0x3fff + 1; // 16384
    static const int INT_TABLE_SIZE = 0xff + 1; //   256
private:
    RanTable *ranTable;

    // declare constructor parameters with some default values, these ones chosen for speed
    // for very high dimensional, use higher values for challenging problem in one dimension only
    int fractalDepth; // maximum recursive depth of 40 suggested for 64-bit architecture
    int density;
    long long index;
    UnitFunction1D *unitFunction;

    // Constructors

    /// <summary>
    /// Create a new 1D fast fractal function generator. </summary>
    /// <param name="unitFunction"> the base function for this generator </param>
    /// <param name="fractalDepth"> recursive depth of fractal - each increment adds detail at half the scale
    /// (double the resolution).
    /// Must be between 1 and 2^64 although in practice maximum supported by IEEE 64-bit floating point
    /// is in the low 40s. Recommend maximum of 40. </param>
    /// <param name="density"> average number of base functions per unit area at each resolution </param>
    /// <param name="index"> the sequence number of this surface (for the given fractal depth and density) </param>
public:
    FractalFunction1D(UnitFunction1D *unitFunction, int fractalDepth, int density, long index);

    /// <summary>
    /// Create a new 1D fast fractal function generator using default values (fractal depth = 3). </summary>
    /// <param name="unitFunction"> the base function for this generator </param>
    /// <param name="density"> average number of base functions per unit area at each resolution </param>
    /// <param name="index"> the sequence number of this surface (for the given fractal depth and density) </param>
    FractalFunction1D(UnitFunction1D *unitFunction, int density, long index);

    /// <summary>
    /// Create a new 1D fast fractal function generator using default values (fractal depth = 3, density = 1). </summary>
    /// <param name="unitFunction"> the base function for this generator </param>
    /// <param name="index"> the sequence number of this surface (for the given fractal depth and density) </param>
    FractalFunction1D(UnitFunction1D *unitFunction, long index);
    ~FractalFunction1D();

    /// <summary>
    /// Create a new 1D fast fractal function generator using default values
    /// (fractal depth = 3, density = 1, index = 1). </summary>
    /// <param name="unitFunction"> the base function for this generator </param>
    FractalFunction1D(UnitFunction1D *unitFunction);

    /// <summary>
    /// Create a new 1D fast fractal function generator using default values
    /// (unitFunction = DoubleDip, fractal depth = 3, density = 1, index = 1).
    /// </summary>
    FractalFunction1D();


    /// <summary>
    /// Create a new generator in the same series by resetting the index (faster than creating a new object). </summary>
    /// <param name="index"> the new index (sequence number). </param>
    virtual void setIndex(long index);


    /// <summary>
    /// Evaluate the function at the given co-ordinate. </summary>
    /// <param name="x"> the point at which to evaluate </param>
    /// <returns> the value at that point </returns>
    virtual double evaluate(double x);


private:
    double getDepthLocal(double x, int recDepth, long seed, long span);


    double getDepthWRTSquare(double x, long square, int recDepth, long seed, long span, double scale);


private:
    void InitializeInstanceFields();
};
