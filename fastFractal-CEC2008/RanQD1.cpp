//
//  RanQD1.java
//  FractalFunctions
//
//  Created by Cara MacNish on 20/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.
//
//  Translated in C++ by Giuseppe A. Trunfio on Nov 10, 2013

/**
 * This class implements the (repeatable) "quick and dirty" portable pseudo-random
 * generator with modulus 2<sup>32</sup> in
 * W. H. Press et al, "Numerical Recipes in C: The Art of Scientific Computing",
 * Cambridge University Press, 2nd Ed., 1992, with values recommended by Knuth.
 * It is the fastest generator recommended in the above tome.
 * <p>
 * The original relies on the fact that in C on a 32-bit machine, multiplying two unsigned long ints
 * returns the lower 32 bits of the 64 bit product. Since we are using Java (no unsigned ints)
 * and a 64-bit architecture, we mimic this using 64 bit longs with bit-masking instead.
 * @author {@link <a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>}, University of Western Australia
 * @version 1.0RC1, 7th Nov 2007
 * <br>For the latest version and additional information see the
 * {@link <a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>}
 */

#include "RanQD1.h"

const double RanQD1::MAX_INT = 4294967295.0;

RanQD1::RanQD1(long long seed)
{
    InitializeInstanceFields();
    idum = seed;
    nextLong(); // one multiple to scatter seeds
}

void RanQD1::setSeed(long long seed)
{
    idum = seed;
    nextLong(); // one multiple to scatter seeds
}

long long RanQD1::nextLong()
{
    idum = (A * idum + C) & MASK;
    return idum;
}

double RanQD1::nextDouble()
{
    return nextLong() / MAX_INT;
}

int RanQD1::nextInt(int min, int max)
{
    return min + static_cast<int>(floor(nextDouble()*(max - min + 1)));
}

void RanQD1::InitializeInstanceFields()
{
    idum = 0;
}
