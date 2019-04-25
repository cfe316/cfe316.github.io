---
layout: post
title: "Launching particles with a Maxwellian-influx distribution"
author: "Jacob Schwartz"
categories: blog
tags: [particles]
image: Maxwellian-Laplacian.png
use_math: true
---

Part of my thesis is about simulating gases inside a box: their flow and how they interact with a plasma that is also inside a box. In order to do this, I'm using a FORTRAN code that's been developed by generations of graduate students. 
I ran a simple test case in which I should end up with a gas of a uniform temperature and density inside my box, but to my horror the temperature was lower and the density higher than what I expected. 
In my test case, I model the boundaries box as surfaces beyond which is a resevoir of a stationary gas in thermal equilibrium.
The code works by launching particles from the walls, as if they had come from that reservoir, and tracking their passage through the box.
I realized there was a problem with the part of the code that assigns particles an initial velocity starting from the wall.

## How NOT to launch particles from the wall.

A stationary gas in thermal equilibrium is isotropic: particles have an equal chance of moving in all directions.
In Figure 1, points represent the velocity vector of such particles. Each of $v_x$, $v_y$, and $v_z$ are drawn from a Gaussian distribution with a mean $\mu=0$ and standard deviation $\sigma = 1$. The mean of the square of any velocity component is 1, and since they are independent, the mean of the full $v^2$ is 3.

{% include figure.html url="3dMaxwellian.png" 
caption="Figure 1: The velocities in a gas of (dimensionless) temperature $kT/m = 1$." %}

