//
//  FractalFunction1D.java
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
 * This class is called by {@link FastFractal} to evaluate the components for each dimension to
 * evaluate fast very high-dimensional fractal functions (landscapes).
 * <p>
 * It can also be used directly to generate 1-dimensional fractal landscapes.
 * @author {@link <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>}, University of Western Australia
 * @version 1.0RC1, 7th Nov 2007
 * <br>For the latest version and additional information see the
 * {@link <a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>}
 */


#include "FractalFunction1D.h"
#include "DoubleDip.h"

FractalFunction1D::FractalFunction1D(UnitFunction1D *unitFunction, int fractalDepth, int density, long index)
{
    InitializeInstanceFields();
    this->unitFunction = unitFunction;
    this->fractalDepth = fractalDepth;
    this->density = density;
    this->index = index;
    ranTable = new RanTable(DOUBLE_TABLE_SIZE, INT_TABLE_SIZE, density, index);
}

FractalFunction1D::~FractalFunction1D()
{
    delete ranTable;
    delete unitFunction;
}

FractalFunction1D::FractalFunction1D(UnitFunction1D *unitFunction, int density, long index)
{
    InitializeInstanceFields();
    this->unitFunction = unitFunction;
    this->density = density;
    this->index = index;
}

FractalFunction1D::FractalFunction1D(UnitFunction1D *unitFunction, long index)
{
    InitializeInstanceFields();
    this->unitFunction = unitFunction;
    this->index = index;
}

FractalFunction1D::FractalFunction1D(UnitFunction1D *unitFunction)
{
    InitializeInstanceFields();
    this->unitFunction = unitFunction;
}

FractalFunction1D::FractalFunction1D()
{
    InitializeInstanceFields();
    this->unitFunction = new DoubleDip();
}

void FractalFunction1D::setIndex(long index)
{
    this->index = index;
    ranTable->setSeed(index);
}

double FractalFunction1D::evaluate(double x)
{
    x = fmod(x, 1);

    // note in Java -4.3%1 is -0.3 not 0.7 ie Matlab 'rem' function not 'mod'
    if (x <= 0) // 0 must move to 1, or will be in wrong "square"
    {
        x = x + 1;
    }
    if (fractalDepth < 1) // check for valid depth argument, should never fire
    {
        return 0;
    }
    else // start recursion
    {
        return getDepthLocal(x, 1, index, 1);
    }
}

double FractalFunction1D::getDepthLocal(double x, int recDepth, long seed, long span)
{
    double depth = 0;
    double scale = 1.0 / span;
    long  square = static_cast<long>(ceil(x*span));
    long  newSeed, square1;
    double x1;
    // get contribution from each of the 3 relevant squares...
    for (int offset = -1; offset < 2; offset++)
    {
        x1 = x;
        square1 = square + offset;
        if (square1 == 0) // wrap to r.h.s.
        {
            square1 = span;
            x1 = x1 + 1;
        }
        else if (square1 > span) // wrap to l.h.s.
        {
            square1 = 1;
            x1 = x1 - 1;
        }
        depth = depth + getDepthWRTSquare(x1, square1, recDepth, seed, span, scale); // accumulate contributions
    }
    // now fire recursion to next level down...
    if (recDepth < fractalDepth)
    {
        newSeed = (span + seed) & (DOUBLE_TABLE_SIZE-1); // unique seeds up to random table size
        long newSpan = span << 1; // newSpan = 2^(recDepth-1), bit shift faster
        depth = depth + getDepthLocal(x,recDepth + 1,newSeed,newSpan); // recur to next level
    }
    return depth;
}

double FractalFunction1D::getDepthWRTSquare(double x, long square, int recDepth, long seed, long span, double scale)
{
    double depth = 0;
    long  squareSeed = (square-1); // unique seed for this square
    long  localSeed = (seed + squareSeed) & (DOUBLE_TABLE_SIZE-1); // unique seed for square at depth up to table size
    ranTable->setSeed(localSeed); // apply seed for this square
    int numUnits = ranTable->nextInteger(); // choose number of unit functions from uniform dist whose average is 'density'
    for (int i = 1; i <= numUnits; i++) // get contribution from each
    {
        double diameter = 1 / (2 - ranTable->nextDouble()) * scale; // get diameter from quadratic distribution
        double centre = (square - ranTable->nextDouble()) * scale; // get centre from uniform distribution
        if ((x - centre)*(x - centre) < diameter*diameter / 4) // save making unnecessary call if unit function is too far to affect point
        {
            unitFunction->setCentre(centre); // faster to set individually (believe it or not)
            unitFunction->setScale(diameter);
            depth = depth + unitFunction->getValue(x); // add this unit function's contribution
        }
    }
    return depth;
}

void FractalFunction1D::InitializeInstanceFields()
{
    fractalDepth = 3;
    density = 1;
    index = 1;
}
