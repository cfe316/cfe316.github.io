---
layout: post
title: "Deposition profiles in a circle: a tricky Taylor series"
author: "Jacob Schwartz"
categories: journal
image: exponential_circle_mfp_traces.svg
tags: [integration,mean free path]
use_math: true
---

The previous two posts have been about collision probability integrals in slab geometry.
I've tried extending this to the geometry of an infinite cylinder, but got stuck. 
The quantity of interest to me is the radial profile of interactions coming from particles in a uniform,
isotropic gas outside of the cylinder.
When the particles enter the cylinder, they are absorbed with some mean free path $\lambda$.
This is the same geometry as in [this previous post]({% post_url 2018-11-15-average-cyl-transmission %}), but with a different question---asking where the particles stop instead of how many make it through.
$$
\newcommand{\cancelcolor}[1]{\color{midnightblue}{#1}}
\newcommand{\gone}[1]{\color{midnightblue}{#1}}
\newcommand{\gtwo}[1]{\color{forestgreen}{#1}}
\newcommand{\gthree}[1]{\color{crimson}{#1}}
\newcommand{\gfour}[1]{\color{purple}{#1}}
$$

## Bad news: this post only solves a 2D version of the problem, not 3D

### Let me explain

The bad news is, I have not found a nice form for the deposition profile $d(\rho, \lambda)$.
(I've since found the 1977 paper by Michael Milgram [On the properties of collision probability integrals in annular geometry](https://pubs.aip.org/jmp/article/18/12/2456/225431/On-the-properties-of-collision-probability), which I'm working through.)

Finding the deposition intensity at a given radius requires summing over particles coming from all possible directions. There are many ways to add them up: one could sum over all the possible $v_x$, $v_y$, and $v_z$ that hit the cylinder, for example. 
For this problem, I decompose this into sums over particles that that hit the cylinder coming from
* All possible 'impact parameters' $b$.
* All possible azimuthal angles $\phi$.
* All the polar angles $\beta$ from $0$ to $\pi$.

This blog post only covers the first two, which is hard enough.
Actually, it's mostly about the first integral, over the relevant impact factors. The second is trivial, as by symmetry the deposition from any $\phi$ is the same.
So, we have a flat, circular geometry, shown in Figure 1.


{% include figure.html url="exponential_circle_geometry.png" 
caption="Figure 1: Geometry for computing the intensity of deposition in a ring at radius $\rho$, from sources located outside the disk and which are absorbed in the disk with mean free path $\lambda$, by integrating over the central angle $\theta$." %}

Actually, rather than integrate over $b$ and then somehow distribute that deposition to different radii $\rho$, we integrate over points with the same $\rho$ but different central angles $\theta$.

Note that the figure shows only particles coming from one source angle $\phi$; I'll call the deposition intensity $\underline{i} \equiv i/2\pi$, where $i$ is the intensity from particles coming from all $\phi$ (but still in 2D) and $I$ would be the intensity in the 3D problem (including the effect of the polar angles).

The reduced intensity of the deposition at a given radius $\rho$ is

$$
\begin{equation}
\underline{i} = \frac{1}{\lambda}\frac{1}{2\pi\cancelcolor{\rho}}\int_0^{2\pi}\cancelcolor{\rho}\;d\theta\; \exp\left(-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda\right).
 \tag{1}\label{eq:one}
\end{equation}
$$

I will explain this expression.
The numerator of the quantity in the exponent has two terms.
The term $\sqrt{1-\rho^2\sin^2\theta}$ is the x-location of the right-hand edge of the circle, at the y-value $\rho \sin\theta$.
The term $\rho \cos\theta$ is the x-location of the point on the right at radius $\rho$. The difference between them is the thickness of material that particles from the direction $\phi$ need to travel through to get to the point on the ring.
The leading $1/\lambda$ is the intensity of deposition from particles with mean free path $\lambda$.

The $\color{midnightblue}{\rho}$ in the integral is the standard Jacobian for polar coordinates.
The $1/2\pi\color{midnightblue}{\rho}$ outside the integral is because $\underline{i}$ is the deposition intensity (averaged over all angles); not the total deposition at that radius. The $2\pi\color{midnightblue}{\rho}$ is the circumference of the ring with radius $\rho$.

Cancelling the $\color{midnightblue}{\rho}$, the reduced deposition intensity at a given radius is 

$$
\begin{equation}
\underline{i} = \frac{1}{\lambda}\frac{1}{2\pi}\int_0^{2\pi}\;d\theta\; \exp\left(-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda\right).
 \tag{2}\label{eq:two}
\end{equation}
$$

## As a Taylor series

I express $\underline{i}$ as a Taylor series around $\rho=0$.
Pleasingly, only the even terms are present, and all coefficients are positive.

$$
\begin{equation}
\boxed{
\;\underline{i} = \frac{e^{-1/\lambda}}{\lambda}
\sum_{n \text{ even},\;n\ge0}^{\infty}
\frac{1}{n!} \frac{\rho^n}{\lambda^n}
\sum_{k = 0}^n \frac{\lambda^k}{(n/2)!}
\sum_{c = \max(k - n/2,\,0)}^{\min(k/2, n)} t(n,k,c)
\;}

\end{equation}
$$


<!--

$$
\begin{equation}
T(n,k) \equiv \frac{\lambda^k}{(n/2)!} \sum_{c = \max(k - n/2,\,0)}^{\min(k/2, n)} t(n,k,c)
\end{equation}
$$

and where
-->

where

$$
\begin{equation}
t(n,k,c) \equiv \frac{(k-c)_c}{2^c \,\pi} \binom{k-c-1}{c} \binom{n}{2 c-2 k+n}\Gamma\left(k - c+\frac{1}{2}\right)  \Gamma
   \left(c-k+\frac{n+1}{2}\right) (2k-2c-1)\text{!!}
\end{equation}
$$


where $(x)_y$ is the [Pochhammer symbol](https://en.wikipedia.org/w/index.php?title=Pochhammer_symbol), $\binom{x}{y}$ is a binomial, $\Gamma$ is the Gamma function, and `!!` is the double factorial $(n (n-2) \ldots)$. Here Pochhammer is preferred to a ratio of Gamma functions because it is `1` rather than undefined when $x=y=0$.

**Woof.** Those are quite some coefficients!
Mechanically, it's easy to get the coefficients of the taylor series one by one, via repeated differentiation, but it was more difficult to figure out how to express them in a form such that the coefficient of the $n'$th term can be computed.

In the rest of the post I'll describe how I came up with this formula for the coefficients $\mathscr{T}_n(\lambda)$.

Here by $\mathscr{T}_n(\lambda)$ I mean the polynomial in $\lambda$ such that

$$
\begin{equation}
\underline{i} = \frac{e^{-1/\lambda}}{\lambda}
\sum_{n \text{ even},\;n\ge0}^{\infty}
\frac{1}{n!}\frac{\rho^n}{\lambda^n} \mathscr{T}_n(\lambda).
\end{equation}
$$

This polynomial $\mathscr{T}_n(\lambda)$ handles the interesting complex part of the result;
what's left is a constant in the front and more or less the expected parts of a Taylor series.
Dialing in even further,

$$
\mathscr{T}_n(\lambda) \equiv \sum_{k = 0}^n \lambda^k \mathscr{T}_{n,k}
$$

where

$$
\begin{equation}
\mathscr{T}_{n,k} = \frac{1}{(n/2)!}
\sum_{c = \max(k - n/2,\,0)}^{\min(k/2, n)} t(n,k,c).
\end{equation}
$$

## Prelude: Slight Mechanical Derivative

Equation \eqref{eq:two} can be expressed as

$$ \underline{i} = \frac{1}{2\pi\lambda} \int_0^{2\pi} f\; d\theta$$

where $f \equiv \exp(-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda)$.

Let's expand the integrand in a Taylor series around $\rho=0$:

$$ \underline{i} = \frac{1}{2\pi\lambda} \int_0^{2 \pi} \sum_{n=0}^{\infty} \frac{\rho^n}{n!} \left(\left. \frac{d}{d\rho^n} f\right) \right|_{\rho=0}  \; d\theta.$$

Since the integral is over $\theta$ only we can move the integral sign inside the sum:

$$ \underline{i} =  \frac{1}{2\pi\lambda}\sum_{n=0}^{\infty} \frac{\rho^n}{n!} \int_0^{2 \pi} \left(\left. \frac{d}{d\rho^n} f\right) \right|_{\rho=0}  \; d\theta.$$

The task is then to find a general form for the coefficient of $\rho^n$,

$$
\begin{equation}
a_n = \frac{1}{2\pi\lambda\,n!}\int_0^{2 \pi} \left(\left. \frac{d}{d\rho^n} f\right) \right|_{\rho=0}  \; d\theta.
 \tag{3}\label{eq:generalTaylorForm}
\end{equation}
$$

The first few are

{% include margin-note.html id="1" content="
In this notation,
$a_n = e^{-1/\lambda}\lambda^{-(n+1)}\mathscr{T}_n(\lambda)/n!$.
" %}

$$
\begin{aligned}
a_0 &= e^{-1/\lambda} / \lambda \\
a_2 &= (1 + \lambda)\, e^{-1/\lambda} / 4\,\lambda^3\\
a_4 &= (1 + 2\lambda + 3 \lambda^2 + 3 \lambda^3) \, e^{-1/\lambda} / 64\,\lambda^5\\
a_6 &= (1+3 \lambda +9 \lambda ^2+24 \lambda ^3+45 \lambda ^4+45 \lambda ^5) e^{-1/\lambda} / 2304\,\lambda^7 \\
\end{aligned}
$$

The integers in the denominator are given by $2^n((n/2)!)^2$.
The process that forms the polynomials in the numerator is more complex, and not amenable by a simple lookup in [OEIS](https://oeis.org/A002454) (A002454, with $n \to n/2$). 
This process will be discussed in the next section.


## The polynomial game

Let's look at the integral of Equation \eqref{eq:generalTaylorForm}.
Taking derivatives of $f$ yields sums of terms of the form

$$
\begin{equation}
w\left(\frac{\rho^r \sin^{2s}\theta \cos^q\theta}{\lambda^g(1 - \rho^2 \sin^2\theta)^m}\right)f
 \tag{4}\label{eq:rsqform}
\end{equation}
$$

where $w$, $r$, $s$, $q$, and $g$ are non-negative integers, and $m$ can be an integer or half-integer $(0, 1/2, 1, 3/2, \ldots)$.

For example the first derivative is
{% include margin-note.html id="2" content="
The leading $f^{-1}$ is for brevity of the displayed math.
It's not part of the integral or Taylor series process.
" %}

$$f^{-1} f' = \frac{\cos\theta}{\lambda }+\frac{\rho  \sin ^2\theta}{\lambda  (1-\rho ^2 \sin ^2\theta)^{1/2}}$$

and the second derivative is

$$
\begin{aligned}
f^{-1} f'' &= \frac{\sin ^2\theta}{\lambda  \left(1-\rho ^2 \sin ^2\theta\right)^{3/2}} + \frac{2 \rho \cos\theta \sin ^2\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)^{3/2}} - \frac{2 \rho ^3 \cos \theta \sin ^4\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)^{3/2}} \\
&+ \frac{\cos ^2\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)} - \frac{\rho ^2 \cos ^2\theta \sin ^2\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)} + \frac{\rho ^2 \sin ^4\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)}.
\end{aligned}
$$

This intermediate step can generate many terms, but the expression simplifies considerably when setting $\rho \to 0$:

$$
\begin{aligned}
\left.f'\right|_{\rho=0} &= \frac{\cos\theta}{\lambda} e^{-1/\lambda} \\
\left.f''\right|_{\rho=0} &= \left(\frac{\cos^2\theta}{\lambda^2} + \frac{\sin^2\theta}{\lambda}\right)e^{-1/\lambda}.
\end{aligned}
$$

Integrating from 0 to $2\pi$, these evaluate to

$$
\begin{aligned}
b_1 \equiv \int \left.f'\right|_{\rho=0} &= 0 \\
b_2 \equiv \int \left.f''\right|_{\rho=0} &= e^{-1/\lambda}\frac{\pi}{\lambda^2}\left(1 + \lambda\right) = a_2(2\pi\lambda\,2!) .
\end{aligned}
$$


Examining the derivative process more closely, let's take the derivative with respect to $\rho$ of a general one of the Expression $\eqref{eq:rsqform}$-type terms.
The existing $\sin^{2s}\cos^q/\lambda^g$ are independent of $\rho$, so I will ignore them.

$$
\begin{aligned}
f^{-1}\frac{d}{d\rho} \frac{\rho^r}{\left(1 - \rho^2 \sin^2 \right)^m} f &=
\frac{\rho^{r+1}\,\sin^2}{\lambda (1 - \rho^2 \sin^2)^{m+1/2}} 
+ \frac{\rho^r \cos}{\lambda (1 - \rho^2 \sin^2)^{m}} \\
&+ 2m\frac{\rho^{r+1}\sin^2}{(1 - \rho^2 \sin^2)^{m+1}}
+ r\frac{\rho^{r-1}}{(1 - \rho^2 \sin^2)^{m}}.
\end{aligned}
$$

Thinking of a term with some combination of $\{r, m, s, q, g\}$ as a point in a 5-d space, the derivative operation on a point yields up to four new points with associated weights $1$, $1$, $2m$, and $r$:

```
The four derivative steps:
    {r+1, m+1/2, s+1, q  , g+1}   # A
+   {r  , m    , s  , q+1, g+1}   # B
+ 2m{r+1, m+1  , s+1, q  , g }    # C
+  r{r-1, m    , s  , q  , g }    # D
```

{% include margin-note.html id="3" content="
The base $f$ has $r=m=s=q=g=0$.
" %}

I say _up to_ four new points because if $m=0$ (as in the base $f$) or if $r=0$ (also in the base $f$, and along paths where $r$ increased and then decreased to zero) then those terms will not be generated. I've labeled the four new points, or steps, as $A, B, C, D$.

For convenience I define 

$$
b_n \equiv 2\pi\lambda \,n!\,a_n = \int_0^{2\pi}\left(\left. \frac{d}{d\rho^n} f\right) \right|_{\rho=0}  \; d\theta.
$$

We will find $b_n$ by taking $n$ steps along all the possible derivative paths, or, the valid combinations of $A$, $B$, $C$, $D$ of length $n$; then setting $\rho$ to zero, and finally integrating over $\theta$ to end up with a polynomial in terms of $1/\lambda$.
The coefficient of each power of $\lambda$ is formed from the sum of the weights of all the valid paths that generate $\lambda$'s of that power.

### A very general start

One initial way to think about constructing the terms of the n'th derivative of $f$ is something like

$$
\begin{equation}
\frac{d }{d\rho^n }f = \sum_{r,m,s,q,g} \left(\frac{\rho^r \sin^{2s}\theta \cos^q\theta}{\lambda^g(1 - \rho^2 \sin^2\theta)^m}\right)f
\end{equation}
$$

### Restrictions on useful paths

#### Useful terms have $r=0$. 
The fact that we set $\rho$ to zero before integrating means that the only terms left to integrate will ones with the power $r=0$. 
Steps $A$ and $C$ both increase $r$, and step $D$ decreases it; $B$ leaves it unchanged.
So, a term will remain only if the sum of the number of $A$s and $C$s equals the number of $D$s along a path.

$$\boxed{a + c -d = 0}$$

We use this restriction in the subsection immediately below.

### Weights of paths to get to a point

#### Factor from the order of (A|C)'s vs D's
As discussed above, the number of $D$s must equal the number of $A$s and $C$s.
The *order* of the $A$s, $C$s, and $D$s contributes an extra weight to path because of the coefficient $r$ whenever step $D$ is taken.
For example, $AAADDD$ has a weight of $3∗2∗1=6$ but $ADADAD$ has a weight of $1∗1∗1=1$.
{% include margin-note.html id="3" content="
$A$ and $C$ are equivalent here as they both increment $r$.
" %}
As mentioned before, not all paths are valid: all the prefix substrings must have the number of ($A$'s or $C$'s) greater than or equal to the number of $D$'s, or else the weight will crash down to zero (from differentiating a constant).

For paths of length n=4

```
AADD : weight of 2
ADAD : weight of 1
ADDA : invalid
DAAD : invalid
DADA : invalid
DDAA : invalid
```

I found experimentally that the sum of the weights of all these paths is 

$$ \gone{(2(a+c)-1)!! }$$

where `!!` is the double factorial.

### Factor from order of A's and C's
In the above process, the order of $A$s and $C$s does not matter.
But, there is also a weight for each path from the coefficients $2m$, which depends *only* on the order of the $A$s and $C$s.
Each $A$ increases $m$ by $1/2$ and each $C$ multiplies the weight by $2m$ and increments $m$ by $1$.

I simulated this process and found that the pattern of the sum of the weights over all permutations of $A$ and $C$ here was basically the same as the [Bessel numbers](https://oeis.org/A100861), but with the indices shifted. This yields a factor of 

$$\gtwo{\frac{(a + 2c - 1)!}{2^c (a-1)!c!}}.$$

### Factor from insertion of B's
The steps $B$ can be inserted anywhere in the above sequences without changing the outcome of the $r$-weights or $m$-weights.
This contributes a factor

$$\gthree{\binom{2(a+c) + b-1}{b-1}}.$$


### Restriction
I noticed experimentally that the $n=\text{odd}$ derivatives are always zero, because they have only $q=\text{odd}$ powers of cosine. 
But the $n$-even derivatives _also_ have only even powers of cosine.
This is because if $n = a + b + c + d$ is even, and $(a+c) = d$ then $a + c + d$ is even, so $b$ (which contributes a power of cosine) is even.

Therefore before integrating over $\theta$, $b_n$ is a sum of terms $\cos^n, \cos^{n-2}\sin^{2}, \cos^{n-4}\sin^4, \ldots \sin^{n}$.

### Thing
$A$ and $C$ both bring a factor of $\sin^2$, while $A$ and $B$ bring a factor of $1/\lambda$. Also $a + c = d$.
For any given $\lambda^{k-n}$ *and* combination of $\cos$ and $\sin$ there is a limited range of $c$ which end up contributing.
This leads to the restrictions on the sum to form the term with power $\lambda^{k-n}$, where $c$ ranges from $\max(k - n/2, 0)$ to $\min(k/2, n)$. 

### Thing
Finally, integrate. The integral of $\sin^{2s}\theta\cos^q\theta$ contributes a factor

$$\gfour{\frac{2 \Gamma \left(\frac{q+1}{2}\right) \Gamma
   \left(s+\frac{1}{2}\right)}{\Gamma
   \left(\frac{q}{2}+s+1\right)}}$$

(and remember that the useful $q$ here are always even).
