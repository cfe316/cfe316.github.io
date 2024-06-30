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
isotropic gas outside of the cylinder, where the particles
are absorbed with some mean free path $\lambda$ when they enter the cylinder. 
This is the same geometry as in [this earlier post]({% post_url 2018-11-15-average-cyl-transmission %}), but with a different question---asking where the particles stop instead of how many make it through.
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
The $\sqrt{1-\rho^2\sin^2\theta}$ is the x-location of the right-hand edge of the circle, at the y-value $\rho \sin\theta$.
The $\rho \cos\theta$ is the x-location of the point on the right at radius $\rho$. The difference between them is the thickness of material that particles from the direction $\phi$ need to travel through to get to the point on the ring.
The leading $1/\lambda$ is the intensity of deposition from particles with mean free path $\lambda$.

The $\color{midnightblue}{\rho}$ in the integral is the standard Jacobian for polar coordinates.
The $1/2\pi\color{midnightblue}{\rho}$ outside the integral is because $\underline{i}$ is the deposition intensity (averaged over all angles); not the total deposition at that radius. The $2\pi\color{midnightblue}{\rho}$ is the circumference of the ring with radius $\rho$.

Cancelling the $\color{midnightblue}{\rho}$, the reduced deposition intensity at a given radius is 

$$
\begin{equation}
\underline{i} = \frac{1}{2 \pi \lambda}\int_0^{2\pi}\;d\theta\; \exp\left(-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda\right).
 \tag{2}\label{eq:two}
\end{equation}
$$

For brevity I will define

$$
\begin{equation}
f \equiv \exp\left(-\left(\sqrt{1 - \rho^2 \sin^2\theta} - \rho \cos\theta\right)/\lambda\right)
\end{equation}
$$

so Equation \eqref{eq:two} can be expressed as

$$ \underline{i} = \frac{1}{2\pi\lambda} \int_0^{2\pi} f\; d\theta.$$

Unfortunately, we cannot evaluate the integral in terms of elementary functions, or
as far as I can tell, even in terms of complicated non-elementary functions like the Meijer-G function.

Instead, I will find a Taylor series for $\underline{i}$ around $\rho=0$. 
I will do this by expanding $f$ in a Taylor series and integrating each term.

## Outline

I will express $f$ as a Taylor series:

$$ f = \sum_{n=0}^\infty j_n(\lambda, \theta) \frac{\rho^n}{n!} $$

where $j_n(\lambda, \theta)$ are functions of $\theta$ and $\lambda$.
By the standard Taylor series construction,

$$ j_n \equiv \left(\left.\frac{d}{d\rho^n}f\right)\right|_{\rho=0} .$$

Then 

$$\underline{i} = \frac{1}{2\pi\lambda} \int d\theta \sum_n j_n(\lambda,\theta) \frac{\rho^n}{n!}.$$

I reverse the order of the sum and integral (and move left the $\rho^n/n!$) to get a Taylor series for $\underline{i}$,

$$\underline{i} = \frac{1}{2\pi\lambda} \sum_n\frac{\rho^n}{n!} l_n(\lambda)$$

where I'm defining

$$ l_n(\lambda) \equiv \int_0^{2\pi} j_{n}(\lambda,\theta)\;d\theta .$$

I will show that each $j_n(\lambda, \theta)$ is a sum of terms of the form

$$ w\, e^{-1/\lambda} \frac{\sin^{2s}\theta\cos^b\theta}{\lambda^g}$$

where $w$, $s$, $b$, and $g$ are non-negative integers, and further, $b$ is even.

Then, $l_n$ is a polynomial in $1/\lambda$,

$$ l_n(\lambda) = e^{-1/\lambda} \pi \sum_{k=0}^{n-1} \frac{u_{n,k}}{\lambda^{n-k}}$$

where $u_{n,k}$ are ratios of integers.
The coefficients $w$ and then $u_{n,k}$ will be found by carefully analyzing the 
derivative-taking process.
It will turn out that $u_{n,k}$ is best expressed as a sum over a third index, which I'll call $c$.


## The derivative-taking game

Taking repeated derivatives of $f$ with respect to $\rho$ yields sums of terms of the form

$$
\begin{equation}
w\left(\frac{\rho^r \sin^{2s}\theta \cos^b\theta}{\lambda^g(1 - \rho^2 \sin^2\theta)^m}\right)f
\tag{3}\label{eq:rsqform}
\end{equation}
$$

where $w$, $r$, $s$, $b$, and $g$ are non-negative integers, and $m$ can be an integer or half-integer $(0, 1/2, 1, 3/2, \ldots)$.

For example the first derivative is
{% include margin-note.html id="1" content="
The leading $f^{-1}$ is for brevity.
It's not part of the integral or Taylor series process.
" %}

$$f^{-1} f' =
\frac{\rho  \sin ^2\theta}{\lambda  (1-\rho ^2 \sin ^2\theta)^{1/2}} +
\frac{\cos\theta}{\lambda }
$$

and the second derivative is

$$
\begin{aligned}
f^{-1}f'' &=
\frac{\rho ^2 \sin ^4\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)}
+\frac{2 \rho  \cos \theta \sin ^2\theta}{\lambda ^2 \left(1-\rho ^2 \sin ^2\theta\right)^{1/2}} \\
&+ \frac{\rho ^2 \sin ^4\theta}{\lambda  \left(1-\rho ^2 \sin ^2(\theta)\right)^{3/2}}
  +\frac{\sin ^2\theta}{\lambda  \left(1-\rho ^2 \sin ^2\theta\right)^{1/2}}
+\frac{\cos ^2\theta}{\lambda ^2} \\
\end{aligned}
$$

Examining the derivative process more closely, let's take the derivative with respect to $\rho$ of a general one of the Expression $\eqref{eq:rsqform}$-type terms.
The existing $w \sin^{2s}\cos^b/\lambda^g$ are independent of $\rho$, so I will ignore them.

$$
\begin{equation}
\begin{aligned}
f^{-1}\frac{d}{d\rho} \frac{\rho^r}{\left(1 - \rho^2 \sin^2 \right)^m} f &=
\frac{\rho^{r+1}\,\sin^2}{\lambda (1 - \rho^2 \sin^2)^{m+1/2}} 
+ \frac{\rho^r \cos}{\lambda (1 - \rho^2 \sin^2)^{m}} \\
&+ 2m\frac{\rho^{r+1}\sin^2}{(1 - \rho^2 \sin^2)^{m+1}}
+ r\frac{\rho^{r-1}}{(1 - \rho^2 \sin^2)^{m}}.
\end{aligned}
\tag{4}\label{eq:foursteps}
\end{equation}
$$

Thinking of a term with some combination of $\{r, m, s, b, g\}$ as a point in a 5-d space, the derivative operation on a point with initial weight $w=1$ yields up to four new points with associated weights $1$, $1$, $2m$, and $r$:

```
The four derivative steps:
    {r+1, m+1/2, s+1, b  , g+1}   # A
+   {r  , m    , s  , b+1, g+1}   # B
+ 2m{r+1, m+1  , s+1, b  , g }    # C
+  r{r-1, m    , s  , b  , g }    # D

# These are in the same order as the terms in Eq (4).
```

{% include margin-note.html id="2" content="
The base $f$ has $r=m=s=b=g=0$.
" %}

I say _up to_ four new points because if $m=0$ (as in the base $f$) or if $r=0$ (also in the base $f$, and along paths where $r$ increased and then decreased to zero) then those terms will not be generated. I've labeled the four new points, or steps, as $A, B, C, D$.

Here's a digram that shows this in action:


{% include figure.html url="exponential_circle_derivatives_tree.svg" 
caption="Figure 2: Each term is generated by one of the four derivative steps.
Any combination of steps like the set $\{A,B\}$ will contribute to the same term: $AB$ and $BA$ end up in the same place, though they might contribute different weights. 
Not all terms have descendents for all steps: taking steps $C$ or $D$ from the root $1$ or $\cos/\lambda$ yields $0$.
(I've dropped the final factor $f$ for all terms here.)
 "%} 

Using a computer algebra system like Mathematica, it's not difficult to take derivatives, then set $\rho$ to 0, then integrate over $\theta$ to find the polynomial $l_n(\lambda)$.[^mechanical]
But, it would be still be interesting to come up with a closed-form expression for $l_n(\lambda)$.

[^mechanical]: Mathematica snippet to generate $l_n(\lambda)$
    ```
    f = Exp[-(Sqrt[1 - ρ^2 Sin[θ]^2] - ρ Cos[ρ])/λ];
    der[n_] := (D[f, {ρ, n}] // Expand)
    j[n_] := (der[n]) /. ρ -> 0
    l[n_] := Integrate[j[n], {θ, 0, 2 π}]
    ```

In order to do this, let's make some observations about the patterns of derivatives, or the patterns of steps, as in the diagram.

### Only r=0 terms survive in $j_n$
The $n'th$ derivative is the sum of all permutations of length $n$ of $\{A,B,C,D\}$.

$$
\begin{equation}
n = a + b + c + d
\tag{5}\label{eq:nabcd}
\end{equation}
$$

After taking the derivatives, the next step in forming $j_n$ is to set $\rho \to 0$.
Only terms with the exponent $r=0$ (in $\rho^r$) will survive.
The root term $1$ has $r=0$; steps $A$ and $C$ both increment $r$ and step $D$ decremenets it.
Therefore, any surviving term will have

$$
\begin{equation}
a + c = d
\tag{6}\label{eq:acd}
\end{equation}
$$

where $a$ is the number of steps $A$, and so on.

### Only even-$n$ $l_n(\lambda)$ are nonzero.

Step $B$ contributes a factor $\cos\theta$ to the numerator.
When this is integrated from $0$ to $2\pi$, terms with odd powers of cosine will be cancelled out; only terms with even $b$ will survive.
As shown above, surviving terms have $a + c = d$, so $a + c + d$ is even.
If $n$ is odd, $(a + c +d) + b$ must be odd as well, so $b$ will be odd.
Conversely, if $n$ is even $b$ will be even. 

Therefore $l_n(\lambda)$ will be zero for all odd $n$.
This means the Taylor series of $\underline{i}$ only has even powers of $\rho$!

Also, all the terms contributing to even $l_n$ will have $\sin^{2s}\cos^{b}$ with even $b$.

### Restrictions on valid $c$ for given $n, k$

The polynomial $l_n(\lambda)$ is a sum of terms with different powers of $1/\lambda^{n-k}$.
For a given $n$ and $k$ there are restrictions on $c$, the number of steps $C$ which can be taken to contribute to a term with $1/\lambda^{n-k}$.

For example, by inspection of the table of steps above, only $A$ and $B$ increase $g$, that is, contribute a power of $1/\lambda$.
Thus,

$$
\begin{equation}
n-k = a+b.
\tag{7}\label{eq:nkab}
\end{equation}
$$

Using Equations \eqref{eq:nabcd}, \eqref{eq:acd}, and \eqref{eq:nkab}, together with the inequalities

{% include margin-note.html id="3" content="
The restriction on $b$ is less harsh than on $a, c, d$.
Some of these 8 inequalities are redundant.
" %}

$$
\begin{aligned}
0 &\le a \le n/2 \\
0 &\le b \le n  \\
0 &\le c \le n/2  \\
0 &\le d \le n/2
\end{aligned}
$$

we find that for a given $(n,k)$,

$$ \max(0, k-n/2) \le c \le k/2.$$

From the equations listed above we can express $a$, $b$, and $d$ in terms of $n$, $k$, and $c$:

{% include margin-note.html id="4" content="
We could choose to put things in terms of $n$, $k$ and $b$
but that comes out messier.
" %}

$$
\begin{aligned}
a &= k - 2c \\
b &= n - 2k + 2c \\
d &= k - c.
\end{aligned}
\tag{8}
$$

### Weights of paths to get to a term

#### Factor from the order of (A|C)'s vs D's
From Equation \eqref{eq:acd}, $a + c = d$.
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
ADDA : weight of 0
DAAD : weight of 0
DADA : weight of 0
DDAA : weight of 0
```


I found experimentally[^acd] that the sum of the weights of all these paths matches
[OEIS A001147](https://oeis.org/A001147), which is

$$ (2(a+c)-1)!! = \gone{(2k-2c-1)!!}$$

where `!!` is the double factorial.

[^acd]: Mathematica snippet for the (A|C,D) process
    ```
    o = {1, 0};
    ac = {#[[1]], #[[2]] + 1} &;
    d = {Times @@ #, #[[2]] - 1} &;
    p[i_] := Permutations[Flatten[{ac, d}~Table~i], {2 i}]
    sim[i_] := Total[#[[1]] & /@ Map[(Composition @@ #)@o &, p[i]]]
    sim /@ Range[0,6]
```

### Factor from order of A's and C's
In the above process, the order of $A$s and $C$s does not matter.
But, there is also a weight for each path from the coefficients $2m$, which depends *only* on the order of the $A$s and $C$s.
Each $A$ increases $m$ by $1/2$ and each $C$ multiplies the weight by $2m$ and increments $m$ by $1$.

I simulated[^aAndC] this process and found that the pattern of the sum of the weights over all permutations of $A$ and $C$ here was basically the same as the [Bessel numbers](https://oeis.org/A100861), but with the indices shifted. This yields a factor of 

$$\frac{(a + 2c - 1)!}{2^c (a-1)!\,c!} =
\frac{\Gamma(k)}{2^c\, c!\, \Gamma(k-2c)}=
\gtwo{\frac{(k-2c)_{2c}}{2^c\,c!}}
$$

where $(k-2c)_{2c}$ is a [Pochhammer symbol](https://en.wikipedia.org/w/index.php?title=Pochhammer_symbol). The latter form is preferred over the ratio of two Gamma functions to avoid an intederminate expression when $k=c=0$.

[^aAndC]: Snippet for the A,C process
    ```
    o = {1, 0};
    A = {#[[1]], #[[2]] + 1/2} &;
    Cc = {2 Times @@ #, #[[2]] + 1} &;
    perms[a_, c_] := Permutations[Join[A~Table~a, Cc~Table~c], {a + c}];
    sim[a_, c_] := Total[#[[1]] & /@ Map[(Composition @@ #)@o &, perms[a, c]]];
    ```

### Factor from insertion of B's
The steps $B$ can be inserted anywhere in the above sequences without changing the outcome of the $r$-weights or $m$-weights.
Thus the number of paths is multiplied by

$$\binom{2(a+c) + b}{b} =
\gthree{\binom{n}{n+2c-2k}} 
$$

where $\binom{x}{y}$ is a binomial.

### Harmonic factors

Steps $A$ and $C$ contribute factors of $\sin^2\theta$ while step $B$ contributes $\cos\theta$.

Thus the harmonic part will be 

$$\sin^{2(a+c)}\theta \cos^b\theta =
\sin^{2k-2c}\theta \cos^{n-2k+2c}\theta
.$$

### Setting $\rho$ to zero

When setting $\rho$ to 0 each term like that in \eqref{eq:rsqform} becomes

$$
w \frac{\sin^{2s}\theta \cos^b \theta}{\lambda^g}e^{-1/\lambda} 
 = w e^{-1/\lambda} \frac{\sin^{2k-2c}\theta \cos^{n-2k+2c}\theta}{\lambda^{n-k}}
$$

## Construct j

Therefore

$$ 
\begin{aligned}
j_n(\lambda,\theta) =
e^{-1/\lambda}
\sum_{k=0}^n 
\frac{1}{\lambda^{n-k}}
\sum_{c=\max(0,k-n/2)}^{k/2} &
\gone{(2k-2c-1)!!}
\gtwo{\frac{(k-2c)_{2c}}{2^c\,c!}} \\
&
\gthree{\binom{n}{n+2c-2k}} 
\sin^{2k-2c}\theta \cos^{n-2k+2c}\theta \\
\end{aligned}
$$

## Integrate to get $l_n(\lambda)$

For even $b$,

$$ \int_0^{2\pi}
\sin^{2s}\theta \cos^b\theta
= \frac{2\, \Gamma \left(\frac{b+1}{2}\right) \Gamma
   \left(s+\frac{1}{2}\right)}{\Gamma
   \left(\frac{b}{2}+s+1\right)}
$$

Here $s = a + c = k-c$ and $b=n-2k+2c$. Therefore

$$
\begin{aligned}
\frac{2\, \Gamma \left(\frac{b+1}{2}\right) \Gamma
   \left(s+\frac{1}{2}\right)}{\Gamma
   \left(\frac{b}{2}+s+1\right)}
& \to
\frac{2 \Gamma \left(\frac{1}{2}-c+k\right)
\Gamma \left(\frac{1}{2}+c-k+\frac{n}{2}\right)}
{\Gamma \left(1+\frac{n}{2}\right)}
 \\
& =
\gfour{
\frac{2 \Gamma \left(\frac{1}{2}-c+k\right)
\Gamma \left(\frac{1}{2}+c-k+\frac{n}{2}\right)}
{(n/2)!}
}
\end{aligned}
$$

and

$$

\boxed{
\begin{aligned}


l_n(\lambda) &=
\frac{\gfour{2} e^{-1/\lambda}} {\gfour{(n/2)!}}
\sum_{k=0}^{n-1}
\frac{1}{\lambda^{n-k}}
\sum_{c=\max(0,\,k-n/2)}^{k/2}
\gone{(2k-2c-1)!!} \\
&
\gtwo{\frac{(k-2c)_{2c}}{2^c\,c!}}
\gthree{\binom{n}{n+2c-2k}} 

\gfour{\Gamma \left(\frac{1}{2}-c+k\right)
\Gamma \left(\frac{1}{2}+c-k+\frac{n}{2}\right)}.

\end{aligned}
}
$$

I've moved the $\gfour{2/(n/2)!}$ out to the front because it does not depend on $k$ or $c$.

The first few $l_n$ are

$$
\begin{aligned}
l_0 &= 2 e^{-1/\lambda } \pi  \\
l_2 &= \frac{e^{-1/\lambda } \pi  (1+\lambda )}{\lambda ^2} \\
l_4 &= \frac{3 e^{-1/\lambda } \pi  \left(1+2 \lambda +3 \lambda ^2+3 \lambda ^3\right)}{4 \lambda ^4} \\
l_6 &= \frac{5 e^{-1/\lambda } \pi  \left(1+3 \lambda +9 \lambda ^2+24 \lambda ^3+45 \lambda ^4+45 \lambda ^5\right)}{8 \lambda ^6} \\
l_8 &= \frac{35 e^{-1/\lambda } \pi  \left(1+4 \lambda +18 \lambda ^2+78 \lambda ^3+285 \lambda ^4+810 \lambda ^5+1575 \lambda ^6+1575 \lambda ^7\right)}{64 \lambda ^8}
\end{aligned}
$$

For completeness,

$$

\begin{aligned}


u_{n,k} &=
\frac{\gfour{2}} {\pi \gfour{(n/2)!}}
\sum_{c=\max(0,\,k-n/2)}^{k/2}
\gone{(2k-2c-1)!!} \\
&
\gtwo{\frac{(k-2c)_{2c}}{2^c\,c!}}
\gthree{\binom{n}{n+2c-2k}} 

\gfour{\Gamma \left(\frac{1}{2}-c+k\right)
\Gamma \left(\frac{1}{2}+c-k+\frac{n}{2}\right)}.

\end{aligned}
$$

## Finally $\underline{i}$

$$

\boxed{
\begin{aligned}

\underline{i} &=
\frac{e^{-1/\lambda}}{\pi}
\sum_{\text{even}\;n \ge 0}
\frac{\rho^n}{n! \gfour{(n/2)!}}
\sum_{k=0}^{n-1}
\frac{1}{\lambda^{n-k+1}}
\sum_{c=\max(0,\,k-n/2)}^{k/2}
\gone{(2k-2c-1)!!} \\
&
\gtwo{\frac{(k-2c)_{2c}}{2^c\,c!}}
\gthree{\binom{n}{n+2c-2k}} 

\gfour{\Gamma \left(\frac{1}{2}-c+k\right)
\Gamma \left(\frac{1}{2}+c-k+\frac{n}{2}\right)}.

\end{aligned}
}
$$

Putting $l_n(\lambda)$ back into the formula for $\underline{i}$ I was able to combine the $\lambda$ with the others in the denominator and cancel a $2$.

## The shape of $\underline{i}$

The cover image shows $f/\lambda$ as a function of $\theta$ for $\lambda=1/2$ (blue) and $\lambda=2$ (green).
Figure 3 shows the $\underline{i}$ that result from integrating over $\theta$.

{% include figure.html url="exponential_circle_mfp_plot.png" 
caption="Figure 3: 
Plot of $\underline{i}$ for $\lambda=(1/2)$ and $2$. For $\lambda=(1/2)$ I also show
the partial sums from Taylor series of order 2, 4, 8, and 16.
The Taylor series converges much more quickly for higher $\lambda$.
 "%} 
