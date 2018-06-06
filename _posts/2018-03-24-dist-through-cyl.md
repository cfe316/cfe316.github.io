---
layout: post
title: "Average distance a particle starting inside an infinite tube travels through the wall of the tube."
author: "Jacob Schwartz"
categories: blog
tags: [integrals]
image: distance_through_cyl.svg
use_math: true
---

A particle is emitted from a uniform random location inside an infinite tube and travels in a uniform random direction. What is the average distance of tube wall it must travel through to escape?

## Solution

\\[ \overline{d} = \frac{2}{3} \left(\left(r-r^3\right) K\left(\frac{1}{r^2}\right)+\left(r^3+r\right) E\left(\frac{1}{r^2}\right)-2\right) \tag{1}\label{eq:sol} \\]

where the inner radius of the tube is $1$, the outer radius is $r$, $K$ and $E$ are the complete [elliptic integrals](https://en.wikipedia.org/wiki/Elliptic_integral) of the [first kind](http://mathworld.wolfram.com/CompleteEllipticIntegraloftheFirstKind.html) and [second kind](http://mathworld.wolfram.com/CompleteEllipticIntegraloftheSecondKind.html), respectively.

### Asymptotics
#### For thin-wall tubes
In the limit of a tube with outer radius $1+\epsilon$, the expression becomes : \\[ \lim_{r \to 1+\epsilon} \overline{d} \to 2 \epsilon + \frac{1}{4} \epsilon ^2 (2 \log (\epsilon )+5-6 \log (2))\\]
#### For thick-wall tubes
In the limit of a thick wall tube, $r \gg 1 $, the leading terms are: \\[ \lim_{r \to \infty} \overline{d} \to \frac{\pi}{2} r - \frac{4}{3} \\]

## Explanation

### Calculate distance through wall as a function of angle and $r'$
Parameterize the particle's track by $t$. Write an expression for the radius $R$ at a given $t$. Solve for t when $R =1$ and when $R = r$. The difference in $t$ is the track distance $d$ through the cylinder.
\\[ R=\sqrt{(r'+t \sin \theta \cos \phi)^2+t^2 \sin
   ^2\theta \sin ^2\phi} \\]
Solve for $t$
\\[ t(R) = \csc \theta \left(\sqrt{R^2+r'^2 \sin^2 \phi}- r' \cos\phi\right) \\]
\\[ d = t(r) - t(1) \\]

\\[ d(\phi, \theta, r', r) = \csc\theta\left(
\sqrt{r^2 - r'^2 \sin^2\phi} - \sqrt{1 - r'^2 \sin^2\phi}\right) \\]
### Set up the integral

Integrate over all directions $(\theta, \phi)$ with geometric factor $\sin\theta$and integrate over all locations in a circle $(r', \phi')$ with geometric factor $r'$. The integral is normalized by total steradians $4\pi$ and circle area $\pi$.

\\[ \frac{1}{4\pi} \frac{1}{\pi} \int_0^{2\pi} \int_0^1 \int_0^{2\pi} \int_0^{\pi} d(\phi, \theta, r', r) \; r' \sin\theta \;d \phi \; d \theta \; dr' \; d\phi' \\]

In Mathematica:
```
Integrate[
 1/\[Pi] r 1/(4 \[Pi])
   d[\[Phi], \[Theta], rp, r] Sin[\[Theta]], {\[Theta], 
  0, \[Pi]}, {\[Phi], 0, 2 \[Pi]}, {rp, 0, 1}, {\[Phi]\[Phi], 0, 
  2 \[Pi]}, Assumptions -> r > 1]
```

The integral evaluates to Equation $(\ref{eq:sol})$.

## An integral avoiding near-polar $\theta$

It may be desired to ignore angles near the poles, where particles would travel a long distance through the tube walls.

If averaging over only the remaining angles is desired, drop the $\frac{1}{4 \pi}$ normalization from the start of the integral and instead do

\\[ \frac{\frac{1}{\pi} \int_0^{2\pi} \int_0^1 \int_0^{2\pi} \int_{\theta_1}^{\theta_2} d(\phi, \theta, r', r) \; r' \sin\theta \;d \phi \; d \theta \; dr' \; d\phi'}{\int_0^{2\pi}\int_{\theta_1}^{\theta_2}\sin\theta\; d \theta \; d \phi} \\]

## Better: an integral avoiding d over a certain distance

Of course, zenith angle $\theta$ is not the only variable that affects distance through the wall. If you wish to more accurately average over only short distances, Mathematica's `ImplicitRegion` function comes in handy. For example:

```
Ω = ImplicitRegion[ d[ph, th, rp, 1.1] <= 0.2, { {ph, 0, 2 Pi}, {rp, 0, 1}, {th, 0, Pi}}];

NIntegrate[2 rp Sin[th] d[ph, th, rp, 1.1], {ph, rp, th} \[Element] Ω]/ NIntegrate[2 rp Sin[th], {ph, rp, th} \[Element] Ω]
```
