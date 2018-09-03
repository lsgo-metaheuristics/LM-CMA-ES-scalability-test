#pragma once

#include <cmath>

//
//  RanQD1.java
//  FractalFunctions
//
//  Created by Cara MacNish on 20/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.

/// <summary>
/// This class implements the (repeatable) "quick and dirty" portable pseudo-random
/// generator with modulus 2<sup>32</sup> in
/// W. H. Press et al, "Numerical Recipes in C: The Art of Scientific Computing",
/// Cambridge University Press, 2nd Ed., 1992, with values recommended by Knuth.
/// It is the fastest generator recommended in the above tome.
/// <para>
/// The original relies on the fact that in C on a 32-bit machine, multiplying two unsigned long ints
/// returns the lower 32 bits of the 64 bit product. Since we are using Java (no unsigned ints)
/// and a 64-bit architecture, we mimic this using 64 bit longs with bit-masking instead.
/// @author <seealso cref="<a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>"/>, University of Western Australia
/// @version 1.0RC1, 7th Nov 2007
/// <br>For the latest version and additional information see the
/// <seealso cref="<a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>"/>
/// </para>
/// </summary>
class RanQD1
{

public:
    static const long long MASK = 0xffffffffLL; // lower order 32 bits of long
    static const double MAX_INT; // 2^32-1 (MASK as a double)
    static const long long A = 1664525LL; // suggested by Knuth
    static const long long C = 1013904223LL; // suggested by Lewis

private:
    long long idum;

    /// <summary>
    /// Create a new pseudo-random generator. </summary>
    /// <param name="seed"> the seed </param>
public:
    RanQD1(long long seed);

    /// <summary>
    /// Reset the seed (quicker than creating a new instance). </summary>
    /// <param name="seed"> the seed </param>
    virtual void setSeed(long long seed);

    /// <summary>
    /// Get the next long.
    /// </summary>
    virtual long long nextLong();

    /// <summary>
    /// Get the next double.
    /// </summary>
    virtual double nextDouble();

    /// <summary>
    /// Get an integer with equal probability from [min, max] inclusive.
    /// </summary>
    virtual int nextInt(int min, int max);

private:
    void InitializeInstanceFields();
};
