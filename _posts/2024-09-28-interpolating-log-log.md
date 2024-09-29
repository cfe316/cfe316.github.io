---
layout: post
title: "Interpolating in log-log space is often nice"
author: "Jacob Schwartz"
categories: journal
tags: [interpolation,math,mistakes]
image: interpolating_loglog_data.svg
use_math: true
---

In nuclear and atomic physics it's common to have data spanning many order of magnitude, like cross sections or reaction rates.
If we want to plot, integrate, or generally do calculations on the data, we usually need an interpolation scheme to generate values between the existing data points.
The interpolation scheme one chooses affects the accuracy of the output data, but it can also have significant and *surprising* effects on important properties of the resulting function. This could lead to issues with numerical solvers or optimizers downstream.

In this post I demonstrate the effects of four interpolation methods on some strictly positive data that spans multiple order of magnitude in _x_ and _y_, specially, data on the cross section of the D-T fusion reaction as a function of energy. The interpolating functions are shown on a log-log plot.

* Standard 1D linear interpolation (blue) is the simplest, but the sharp transitions can cause issues with numerical integrators, solvers, or gradient-based optimization schemes. In log-log space these linear segments look like rounded stair steps.
* Standard cubic interpolation (orange) produces a smooth curve, but the resulting interpolating curve has two _new_ bad features: first, even though the x and y data each are strictly increasing functions (in this region of the plot) with positive slopes, the resulting curve has a *negative slope* near 2 keV, in the center of the plot. This could trip up a gradient-based optimizer. Even worse, the function starts to oscillate between positive and **negative values** below 1 keV, even though all the data points have positive $y$.

The stair-step shape and the unexpected negative values and slope are eliminated by performing the interpolation in log-log space.

## Doing interpolation in log-log space

$$ y_\mathrm{new}=\exp\left(\mathrm{interp}_{(\log(x), \log(y))}(\log(x_\mathrm{new})) \right)$$

1. Take the log of the _x_ and _y_ data points. 
2. Take the log of the desired $x_\mathrm{new}$.
3. Do interpolation (linear, cubic, or otherwise) on these log-space values to get $\log(y_\mathrm{new})$.
4. Exponentiate to get the final $y_\mathrm{new}$.

This only works on data with strictly positive *x* and *y*, since $\log$ is undefined for zero or negative arguments.

* The green curve shows linear interpolation in log-log space. This is a very good fit to the data when viewed in log-log space, but the sharp transitions between segments can still cause issues.
* The red surve is cubic interpolation in log-log space. This is the nicest function! It'll have good derivatives everywhere.

## Wrapping up, lessons learned

Standard (linear-space) interpolation is probably fine for data all on the same order of magnitude,
but for strictly positive data spanning **many** orders of magnitude, 
cubic interpolation in log-log space yields nice, smooth functions.

Finally, always plot your functions! I did not expect standard cubic interpolation to have the problems with ringing and negative values.

### Wait, what about back in linear space

Does the data which was interpolated in log-log space look weird when you plot it in linear space? Nope, it's still nice and smooth.

#### The data
If you'd like to try your own interpolation experiments, I got the data from the [National Nuclear Data Center](https://www.nndc.bnl.gov/sigma/getInterpreted.jsp?evalid=19788&mf=3&mt=50),
hosted by Brookhaven National Lab.
The first few points shown here are

```
beam energy / eV,      σ/barns
             100,   2.0469e-56
             200,   7.4327e-39
             300,   4.0555e-31
             400,  1.58620e-26
             500,   2.1029e-23
             600,   4.1704e-21
             700,   2.5141e-19
             800,   6.7828e-18
             900,   1.0321e-16
            1000,   1.0268e-15
            2000,   2.3063e-10
            3000,    4.9369e-8
            4000,    1.1604e-6
            5000,    9.7767e-6
            6000, 0.0000464530
            7000,   0.00015443
            8000,   0.00040383
            9000,  0.000890740
           10000,    0.0017326
```
