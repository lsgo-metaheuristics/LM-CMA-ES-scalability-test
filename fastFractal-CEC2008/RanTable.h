#pragma once

#include "RanQD1.h"

//
//  RanTable.java
//  FractalFunctions
//
//  Created by Cara MacNish on 23/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.

/// <summary>
/// For very high dimensional functions (eg 5M evaluations of a 1000 dimension function), re-seeding
/// random generators, generating the random numbers and mapping them into the right range on the fly
/// is too slow, even for <seealso cref="RanQD1"/>. 
/// <para>
/// This companion class pre-generates the random numbers needed by <seealso cref="FractalFunction1D"/> 
/// (in turn called by <seealso cref="FastFractal"/>) and stores them in lookup tables. Changing the index
/// simply becomes resetting the array index.
/// </para>
/// <para>
/// This results in a small loss of statistical diversity for a considerable gain in speed, even 
/// for quite large tables (eg 16K) since it is only done once for any given function.
/// </para>
/// <para>
/// The user does not need to access this class directly, all values are set by its owner, such
/// as <seealso cref="FractalFunction1D"/>.
/// @author <seealso cref="<a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>"/>, University of Western Australia
/// @version 1.0RC1, 7th Nov 2007
/// <br>For the latest version and additional information see the
/// <seealso cref="<a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>"/>
/// </para>
/// </summary>
class RanTable
{
  private:
  double *doubleTable;
  int doubleTableIndex;
  int doubleTableSize;
  int *intTable;
  int intTableIndex;
  int intTableSize;

  public:
  RanQD1 *ran;

  /// <summary>
  /// Create the tables and populate with appropriate random values. </summary>
  /// <param name="doubleTableSize"> size of the table for double values </param>
  /// <param name="intTableSize"> size of the table for int values </param>
  /// <param name="aveInt"> average int value, this will be set from the density </param>
  /// <param name="index"> the index or seed for the tables </param>
  RanTable(int doubleTableSize, int intTableSize, int aveInt, long long index);

  ~RanTable();

  /// <summary>
  /// Reset the seed.
  /// </summary>
  virtual void setSeed(long long seed);

  /// <summary>
  /// Get the next double.
  /// </summary>
  virtual double nextDouble();

  /// <summary>
  /// Get the next int
  /// </summary>
  virtual int nextInteger();

private:
	void InitializeInstanceFields();
};
