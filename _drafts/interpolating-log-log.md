---
layout: post
title: "Interpolating in log-log space is often nice"
author: "Jacob Schwartz"
categories: journal
tags: [interpolation,math,mistakes]
image: interpolating_loglog_data.svg
use_math: true
---

In nuclear and atomic physics it's common to have data spanning many order of magnitude, for things like cross sections or reaction rates.
If we want to plot, integrate, or generally do calculations on the data, we often need an interpolation scheme to generate values between the existing data points.
The interpolation scheme one chooses affects the accuracy of the output data; but it can also have significant and surprising effects on important properties of the resulting function.

In this post I demonstrate the effects of four interpolation methods on some data that spans multiple order of magnitude in X and Y, specially, data on the cross section of the D-T fusion reaction as a function of energy. The interpolating functions are shown on a log-log plot.

* Standard 1D linear interpolation (blue) is the simplest, but the sharp transitions can cause issues with some numerical integrators, solvers, or gradient-based optimization schemes. In log-log space these linear segments look like rounded stair steps.
* Standard cubic interpolation (orange) produces a smooth curve, but the resulting interpolating curve has two new bad features: first, even though the x and y data each are strictly increasing functions (in this region of the plot) with positive slopes, the resulting curve has a *negative slope* near 2 keV. This could trip up a gradient-based optimizer. Even worse, the function starts to oscillate between positive and *negative values* below 1 keV, even though all the data points have positive $y$.

The stair-step shape and the unexpected negative values and negative slope are solved by performing the interpolation in log-log space. What I mean by this is
