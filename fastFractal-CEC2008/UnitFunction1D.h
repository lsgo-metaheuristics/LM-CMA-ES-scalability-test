#pragma once

#include <string>

//
//  UnitFunction1D.java
//  FractalFunctions
//
//  Created by Cara MacNish on 18/10/07.
//  Copyright 2007 CSSE, UWA. All rights reserved.
//
//  This source is distributed under GPL3.0. See ../index.html
//  for important information on modifying and distributing.

/// <summary>
/// Each fractal function is constructed from base functions or "unit functions" at various scales
/// (a little like wavelets) chosen from appropriate probability distributions to preserve
/// self-similarity.
/// <para>
/// The name <i>unit function</i> comes from the fact that the non-zero portion of each function
/// is constrained to a unit hypercube, centred at the origin. In the 1-d case this is therefore the
/// interval (-0.5, 0.5).
/// </para>
/// <para>
/// This is the superclass of all 1-d unit functions. If you are writing your own unit
/// function you only need to call the constructors from your own constructors using
/// super(), and provide an implementation for getValue(point).
///
/// @author <seealso cref="<a href="http://www.csse.uwa.edu.au/~cara/">Cara MacNish</a>"/>, University of Western Australia
/// @version 1.0RC1, 7th Nov 2007
/// <br>For the latest version and additional information see the
/// <seealso cref="<a href="http://www.cs.bham.ac.uk/research/projects/ecb/">Birmingham Repository</a>"/>
/// </para>
/// </summary>
class UnitFunction1D
{

protected:
    double centre;
    double scale;


    /// <summary>
    /// Construct a default unit function whose centre and scale will be set later - only called
    /// by subclasses.
    /// </summary>
public:
    UnitFunction1D();

    /// <summary>
    /// Construct a unit function - only called by subclasses. </summary>
    /// <param name="centre"> the x-axis value to which the centre of this unitfunction is mapped </param>
    /// <param name="scale"> the factor by which this unit function is scaled </param>
    UnitFunction1D(double centre, double scale);

    /// <summary>
    /// Set the location and scale at which this unit function is applied. </summary>
    /// <param name="centre"> the x-axis value to which the centre of this unitfunction is mapped </param>
    /// <param name="scale"> the factor by which this unit function is scaled </param>
    virtual void setParams(double centre, double scale);

    /// <summary>
    /// Set the location at which this unit function is applied. </summary>
    /// <param name="centre"> the x-axis value to which the centre of this unitfunction is mapped </param>
    virtual void setCentre(double centre);

    /// <summary>
    /// Get the location of this unit function. </summary>
    /// <returns> the x-axis value to which the centre of this unitfunction is mapped </returns>
    virtual double getCentre();

    /// <summary>
    /// Set the scale at which this unit function is applied. </summary>
    /// <param name="scale"> the factor by which this unit function is scaled </param>
    virtual void setScale(double scale);

    /// <summary>
    /// Get the scale of this unit function. </summary>
    /// <returns> the factor by which this unit function is scaled </returns>
    virtual double getScale();


    /// <summary>
    /// Evalutate this unit function at the given x-value. </summary>
    /// <param name="point"> the point at which this function is evaluated </param>
    /// <returns> the value </returns>
    virtual double getValue(double point) = 0;


    virtual double twist(double x, double y);



private:
    void InitializeInstanceFields();
};
