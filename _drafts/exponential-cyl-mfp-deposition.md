---
layout: post
title: "Radial profiles of deposition in a cylinder"
author: "Jacob Schwartz"
categories: journal
image: exponential_disk_mfp_geometry.svg
tags: [integration,mean free path]
use_math: true
---

An infinite cylinder of radius $1$ is surrounded by a uniform, isotropic, collisionless gas of particles.
When particles enter the cylinder region, they are absorbed with some mean free path $\lambda$.
I provide an expression for the radial profile of absorption (deposition) intensity $d(\rho, \lambda)$

This this the third and hopefully final post in a series.
The previous two examined 2D version of this geometry---a disk in a plane---and calculated the fraction of particles that survive crossing the disk, and the radial profile of the deposition intensity.

$$
\newcommand{\cancelcolor}[1]{\color{midnightblue}{#1}}
\newcommand{\gone}[1]{\color{midnightblue}{#1}}
\newcommand{\gtwo}[1]{\color{forestgreen}{#1}}
\newcommand{\gthree}[1]{\color{crimson}{#1}}
\newcommand{\gfour}[1]{\color{purple}{#1}}
\newcommand{\ki}[2]{\mathrm{Ki}_{#1}\left(#2\right)}
$$

## The Bickley-Naylor function

The Bickley-Naylor function is defined as

$$\ki{N}{x} \equiv \int_0^{\pi/2} e^{-x/\cos\theta}\cos^{N-1}\theta\;d\theta.$$

To match $\int_0^\pi (\ldots) e^{-1/\lambda \sin\beta}\;d\beta$ we perform some elementary manipulations: change $x$ to $1/\lambda$, replace $\theta$ with $\beta$, replace the $\cos$ with $\sin$ (equal, by symmetry) and extend from $0\to\pi/2$ to $0\to\pi$ (double, by symmetry). Thus

$$\int_0^{\pi} \frac{e^{-1/\lambda\sin\beta}}{\sin^{1-N}\beta}\;d\beta = 2\ki{N}{\frac{1}{\lambda}}.$$

For the sake of completeness


## Problem background and setup
### Background gas, rate of particles entering cylinder
Assume the background density of particles is $\gone{1}$ per unit length squared and all particles move at a speed of $\gtwo{1}$ unit length per unit time.
The one-way flux across a unit area with a horizontal normal is $\left(\int_0^\pi \sin^2\beta \; d\beta \int_{-\pi/2}^{\pi/2} \gone{1} * \gtwo{1} * \cos\phi\;d\phi\right)/ \left(\int_0^{\pi}\sin\beta \int_{-\pi}^{\pi}\;d\phi\right) = 1/4$ per unit time,
where $\phi$ is the azimuthal angle and $\beta$ is the angle from the pole. 

{% include figure.html url="sin_and_sinsq.svg" 
caption="Figure 2: Vertical angle integral.
 "%} 

In the numerator's integral over $\beta$, one factor of $\sin\beta$ is the standard Jacobian of polar coordinates, and the other accounts for how particles traveling at $\beta$ near $\pi/2$ have a higher horizontal velocity and are more likely to pass through a surface with a horizontal normal.

Since the cylinder has an area of $2\pi$ per unit length, the rate of particles entering the cylinder is $\pi/2$ per unit time.

### Setup for integration
I want to find the intensity of deposition at a given radius $\rho$.
This is the total deposition per unit length and per unit time in the cylindrical shell of radius $\rho$ and width $d\rho$, summed over all central angles $\theta$, divided by the annulus area of $2\pi\rho\,d\rho$.

The deposition intensity is

$$ i(\rho, \lambda) = \frac{\int_{\phi=0}^{2\pi}d\phi \int_{\beta=0}^{\pi}\sin^2\beta\,d\beta\, \int_{\theta=0}^{2\pi}\cancelcolor{\rho}\, d\theta\,y(\phi,\beta,\theta, \rho, \lambda)\,d\rho}{2\pi\cancelcolor{\rho} \,d\rho} $$

where $y(\phi, \beta, \theta, \rho, \lambda)$ is the deposition intensity from particles moving at horizontal angle $\phi$, angle from the pole $\beta$, at radius $\rho$ and central angle $\theta$.
The $\cancelcolor{\rho}$ in the numerator is the standard Jacobian of polar coordinates. 
This will cancel with the $\cancelcolor{\rho}$ in the denominator; the $d\rho$s also cancel.

By symmetry, and to resolve the integral over $\phi$, we note that the radial intensity from all angles $\phi$ is the same;
thus we can consider an equivalent scenario where all particles come from the same direction, but with $2\pi$ times the intensity: $y'(\beta, \theta, \rho, \lambda) = 2\pi\,y(\phi, \beta, \theta, \rho, \lambda)$.

<!---
The cover illustration shows all particles entering from the right, with $\phi=\pi$.
--->
<!---
The bad news is, I have not found a nice form for the deposition profile $d(\rho, \lambda)$.
(I've since found the 1977 paper by Michael Milgram [On the properties of collision probability integrals in annular geometry](https://pubs.aip.org/jmp/article/18/12/2456/225431/On-the-properties-of-collision-probability), which I'm working through.)
-->

{% include figure.html url="exponential_cylinder_mfp_detailed_geometry.svg" 
caption="Figure 1: Geometry detail for the decay distance as a function of $\theta$. At this stage, the integral has been rearranged so that all particles come from the right side of the disk.
 "%} 

For particles coming from a now-fixed $\phi$ and $\beta$, the deposition at the point given by $(\rho, \theta)$ is proportional to $\frac{1}{\lambda \sin\beta}\exp(-d_\mathrm{decay}/\lambda)$, where 

$$d_\mathrm{decay} = \left(\sqrt{1 - \rho^2\sin^2\theta} - \rho \cos \theta\right)/\sin\beta.$$

The leading $1/(\lambda\sin\beta)$ is the intensity of deposition from particles with mean free path $\lambda$ combined with the fact that particles travel across the cylinder at a speed proportional to $1/\sin\beta$.

Cancelling the $\cancelcolor{\rho}$, the intensity of deposition at a given radius $\rho$ is

$$
\begin{equation}
i = \frac{1}{2\pi}\int_{0}^{\pi}\,d\beta \sin^2 \beta \int_0^{2\pi}d\theta\; \frac{1}{\lambda \sin \beta} e^{-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda \sin\beta}.
 \tag{1}\label{eq:one}
\end{equation}
$$

For brevity I will define

$$
\begin{equation}
f \equiv \exp\left(-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda \sin \beta\right)
\end{equation}
$$

so Equation \eqref{eq:one} can be expressed as

$$ 
\begin{equation}
i = \int_{\beta=0}^\pi d\beta\,\sin^2\beta \frac{1}{2\pi\lambda \sin \beta} \int_0^{2\pi} f\; d\theta.
\tag{2}\label{eq:two}
\end{equation}
$$

This is essentially the same as Equation (2) of the earlier post, but inside an integral over $\beta$.

## The solution as a Taylor series

Instead, I found a Taylor series for $i$ around $\rho=0$. 
I did this by expanding $f$ in a Taylor series and integrating each term.

The deposition intensity is $I(\rho,\lambda)$ is

$$
\begin{equation}
\boxed{
d(\rho, \lambda) = \frac{2}{\pi \lambda} \sum_{\text{even } n\ge0} \frac{\rho^n}{n!} \sum_{k=0}^{n-1} \frac{u_{n,k}}{\lambda^{n-k}} \mathrm{Ki}_{2 + k -n}\left(\frac{1}{\lambda}\right)
}
\end{equation}
$$

## Plots of $d(\rho,\lambda)$

### Special cases: center and edge
