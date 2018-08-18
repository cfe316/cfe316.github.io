---
layout: post
title: "Particles passing through an infinite cylinder"
author: "Jacob Schwartz"
categories: journal
tags: [integrals]
image: cards.jpg
use_math: true
---

Part of my thesis work involves modeling a cylindrical plasma, 
and examining what happens when neutral molecules hit it.
If the plasma density is low enough, the neutrals might pass straight through.
If the density is higher, they could be elastically scattered by plasma particles.
Or, if the density and temperature are high, they could be ionized.
If they are ionized, at what radius does that typically occur?
What is the distribution of ionization radii as a function of the incoming particle's velocities and the mean free time before ionization?

This post will attempt to lay the groundwork for examining the above issues,
by calculating the average distance that particles encountering a transparent, infinite cylinder pass through it before leaving.

## Solution
The average distance $$\overline{d}$$ that particles travel through the transparent cylinder of radius 1 is

$$\overline{d} = 2$$


We can perform this in two coordinate systems: first, a more complicated calculation in spherical coordinates, and then in a hybrid of cartesian and spherical.

# Method 1

# Method 2

A second calculation is a bit simpler because it can more easily be broken up into factors.

$$ \overline{d} = \frac{\pi}{2} * \frac{4}{\pi} = 2 $$

The first factor is the average distance the particles travel through the cylinder, perpendicular to the axis. The second factor is the average distance multiplier due to particles' velocity along the axis.

## First factor, perpendicular to the axis

See Figure 1. A particle encountering the cylinder does so with a particular impact parameter $b$. We want to find $\overline{h}$, the average distance that the particle goes through the cylinder.

{% include figure.html url="through_circle.svg" 
caption="Figure 1: Particles passing through a circular slice of a cylinder." %}

Since the particle bath is uniform and isotropic, the distribution of $b$ is uniform between -1 and 1. 
Find $h$ as a function of $b$: using the formula for a circle in cartesian coordinates, $h = 2 \sqrt{1 - b^2}$. Taking the average over $h$ over all the $b$:

$$ \overline{h} =  \frac{\int_{-1}^1 2\sqrt{1- b^2} \; db}{\int_{-1}^{1} \;db} = \frac{\pi}{2} $$

The average distance is just the area of the circle divided by the diameter.

## Second factor, axial

If a particle has a velocity along the axis it will travel an extra distance in that direction. We need to find the distribution of axial velocities among particles that hit the cylinder. We have assumed an isotropic velocity distribution; $$ f(\theta, \phi) = \frac{1}{4 \pi} \sin\theta $$.

Particles that will encounter the cylinder do so proportional to their velocity perpendicular to it, so the distribution of particles that hit the cylinder has an extra factor of $\sin\theta$.

$$f_\text{hit}(\theta) \propto \sin^2 \theta$$

{% include figure.html url="axial_factor_outlined.svg" 
caption="Figure 2: Particles passing through a length ." %}
See Figure 2.
The $d$ for a particle moving through the cylinder at an angle $\theta$ is multiplied by $1/\sin\theta$ compared to if it had $\theta = \pi/2$.

So, the average axial factor is 

$$ \overline{a} = \frac{\int_{0}^\pi \frac{1}{\sin\theta} \sin^2 \theta \; d \theta}{\int_0^\pi \sin^2 \theta \; d \theta} = \frac{4}{\pi} $$

The solution is the product of the perpendicular distance and axial factor:

\\[ \overline{d} = \overline{h} * \overline{a} = \frac{\pi}{2} \frac{4}{\pi} = 2 \tag{1} \\]
