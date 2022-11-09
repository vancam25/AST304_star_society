# <p align="center">AST 304 Project 2 Write Up</p>

### <p align="center">Abram Anderson, Hannah Sullivan, Steven VanCamp, Sanskriti Verma</p>
<br>

## Contents
1. [Parameter Space Analysis](https://github.com/vancam25/AST304_star_society/new/main/Project%20%232#parameter-space-analysis) 
2. [Mass-Radius Table](https://github.com/vancam25/AST304_star_society/new/main/Project%20%232#mass-radius-table)
3. [Comparison To Observations](https://github.com/vancam25/AST304_star_society/new/main/Project%20%232#comparison-to-observations)

## Parameter Space Analysis
Before completing the bulk of the analysis for this project we analyzed the parameter space for the
integration input parameters `xi` and `eta`. By doing this we can verify that the integration is insensitive
to their precise value.

The parameter `xi` or $\xi$ is used to determine the stepsize $h$ for the integration loop. It is updated
during the integration loop following Eq.(1).

$$ H_{r} = 4\pi^{3}\rho $$

$$ H_{p} = \frac{4\pi^{4}\rho}{Gm} $$

$$ h={\xi}min(H_{r},H_{p}) \tag{1} $$

Where $\rho$ is the stellar density at the current mass cordinate. $m$ is the current mass cordinate. 

Eq.(1) is chosen to be this so that the step size is always less than the scale for which the radius or pressure 
changes appreciably.

The parameter `eta` or $\eta$ is used to determine the upper bound of the integration so that we do not
get a negative value for the pressure. Physically $\eta$ is used to chose the 'edge' of the star. We do this
with Eq.(2)

$$ P<\eta P_{c} \tag{2} $$

Where $P$ is the stellar pressure at the current mass coordinate. $P_{c}$ is the stars central pressure.

Parameter Space Exploration for $\xi$ and `eta`



The range where the integration is insentive to the choice of $\xi$ is $0\lessapprox\xi\lessapprox0.5$. 
So we can safely choose $\xi$ to be $0.1$. The range where the integration is insensitive to the 
choice of $\eta$ is $0\lessapprox\eta\lessapprox0.1^{-5}$$

## Mass-Radius Table

## Comparison To Observations
