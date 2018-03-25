---
layout: post
title: "Probability of a radioactive particle decaying <em>after</em> passing through a surface"
author: "Jacob Schwartz"
categories: blog
tags: [integrals]
image: chance_of_radioactive_particle_escaping.svg
use_math: true
---

A radioactive particle with mean lifetime $\tau$ is a distance $z$ from a planar surface in 3D space. It begins traveling in a constant, uniform random direction at constant velocity $v$. What is the chance it passes through the plane before decaying?

Let the angle of $v$ from the normal to the surface be $\theta$ and the azimuthal angle $\phi$.

The chance that the particle has not decayed at time $t$ is $e^{-t/\tau}$.

The length to the plane along the direction $\vec{v}$ is $z/\cos\theta$. Thus the chance that it escapes is

\\[ e^{-z/v\tau\cos\theta} \\]

Let $b$ be the normalized distance $b=z/v \tau$.

We’d like to average over all possible directions of travel.
\\[ \frac{1}{4\pi} \int_0^{2 \pi} \int_0^{\pi/2} e^{-b/\cos\theta} \sin \theta \; d\theta \; d\phi \\]

Notice that the $\theta$ integral is only to $\pi/2$: if the particle travels downward it will never reach the plane, so contributes nothing.

The answer is \\[\boxed{\frac{1}{2}Ei_2(b)}\\] where $Ei_2$ is an [exponential]( http://mathworld.wolfram.com/ExponentialIntegral.html) [integral](https://en.wikipedia.org/wiki/Exponential_integral) function.

## In order to do the integral:

First do the $\phi$ integral.

\\[\frac{1}{2} \int_0^{\pi/2} e^{-b/\cos\theta} \sin \theta \; d\theta\\]

(From this point on we drop the leading $1 / 2$ for brevity.)

Let $x=\cos\theta$, $dx = - \sin\theta d\theta$.

\\[ -\int_1^0 e^{-b/x}\;dx =  \int_0^1 e^{-b/x}\;dx \\]

Let $y=1/x$, $dy=-1/x^2 dx$ or $dy/y^2 = -dx$

\\[ -\int_\infty^1 e^{-b y} y^{-2} \;dy =  \int_1^\infty e^{-b y} y^{-2} \;dy \\]

This is the definition of the second exponential integral Ei$_2(b)$.

Putting back the $1 / 2$,

\\[\boxed{\frac{1}{2}Ei_2(b)}\\]
### More about $Ei$
The definition of the exponential integral is

\\[ E_n(z) = \int_0^\infty e^{-z t} t^{-n} \; d t \\] In [Mathematica](https://reference.wolfram.com/language/ref/ExpIntegralE.html) it is `ExpIntegralE[n, z]`.

It happens that $Ei_2(b) = e^{-b} - b \Gamma(0,b)$ where $\Gamma(s,z)$ is the Upper Incomplete Gamma function.
