---
layout: post
title: "Polar form for an ellipse offset from the origin"
author: "Jacob Schwartz"
categories: journal
tags: [algebra]
image: polar_ellipse.png
use_math: true
---

I recently needed a polar-form equation for an ellipse, where the center of the ellipse is offset from the origin. I couldn't easily find such an equation, so I derived it and am posting it here. (It's easy to find expressions for ellipses where the *focus* is at the origin.)

## Result

The radius as a function of angle is

$$
\rho(\theta) = \frac{b^2 x \cos (\theta ) + a^2 y \sin (\theta )+a b \sqrt{\left(a^2-x^2\right) \sin ^2(\theta )+\left(b^2-y^2\right) \cos ^2(\theta )+2 x y \sin (\theta ) \cos (\theta )}}{a^2 \sin ^2(\theta )+b^2 \cos ^2(\theta )}
$$

where $a$ is the horizontal semi-major axis, $b$ is the vertical semi-major axis, and the center of the ellipse is $(x, y)$.

In Mathematica:
```
d[a_, b_, x_, 
  y_, \[Theta]_] := (a^2 y Sin[\[Theta]] + b^2 x Cos[\[Theta]] + 
   a b  Sqrt[
    Cos[\[Theta]]^2 (b^2 - y^2) + 
     2 x y Sin[\[Theta]] Cos[\[Theta]] + (a^2 - 
        x^2) Sin[\[Theta]]^2])/(b^2 Cos[\[Theta]]^2 + 
   a^2 Sin[\[Theta]]^2)
```
### Limitations

This equation only works where the origin is inside the ellipse: otherwise there are some angles for which there is no radius to the ellipse, and two solutions for the radius elsewhere.

## Hints for the derivation

Solve $\tan \theta = (b \sin t + y) / (a \cos t + x)$ for t. Insert the result for $t$ into
$$d^2 = (a \cos(t) + x)^2 + (b \sin(t) + y)^2$$.

