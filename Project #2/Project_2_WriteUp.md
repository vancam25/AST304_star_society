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

### <p align="center">Figure 1 Parameter Space Exploration for $\xi$ and $\eta$</p>
![ParameterAnalysis](https://user-images.githubusercontent.com/90724182/201415012-922f0855-50c8-4f4b-b1aa-6a544421f95f.svg)

The range where the integration is insentive to the choice of $\xi$ is $0\lessapprox\xi\lessapprox0.5$. 
So we can safely choose $\xi$ to be $0.1$. The range where the integration is insensitive to the 
choice of $\eta$ is $0\lessapprox\eta\lessapprox0.1^{-5}$$

## Mass-Radius Table

Now that we know what a good choice for $\xi$ and $\eta$ are we can show the relation between white dwarf mass and radius. This is shown on Table (1).

We found this relation by supplying our routine with a desired total mass $M_{want}$. From this we can make an estimate for $P_{c}$ using Eq.(3), which comes from virial relations.

$$ P_{c}^{guess} = (\frac{G^{5}}{K_{e}^{4}})(M\mu_{e}^{2})^{\frac{10}{3}} \tag{3} $$

Where $\mu_{3}$ is the nucleon-electron ration.

After making this guess we run the integration loop supplying our found values for $\xi$, $\eta$, and $P_{c}$. We choose a maximum number of integration steps to be 100000 since this allows the time for the integration to converge and our computers were able to handle it. A 'core' mass is also supplied, and is chosen to be a small fraction ($10^{-4}$) of the desired mass, this gives a lower bound for the integration. 

With this the integration runs, until the condition specified by Eq.(2) is met and the 'edge' of the white dwarfs atmosphere is reached. The results of these runs are shown in Table (1).

### <p align="center">Table 1 Mass-Radius Relationship for White Dwarves</p>
![mass_radius_relationship](https://user-images.githubusercontent.com/90724182/201414979-e0103d17-3546-4680-bb1a-e40103b21dfe.png)

The final column in these tables, $\rho_{c}/[3M/(3*\pi*R^{3}]$, has the same value for all masses. Its the ratio of central density to overall density. This value dose not seem excessively large for a compact object.

## Comparison To Observations
<br>

### <p align="center">Figure 2 Mass and Radius Comparison to Observations</p>
![Mass_Radius_Plot](https://user-images.githubusercontent.com/90724182/201415029-96bf1ed4-ef61-417e-bd23-cb41dfda6ab0.svg)

Our data is within the range shown by the other data sets, HST & FUSE. However the points at lower masses $<=0.2$ are difficult to include this statement since there is no supporting data from HST or FUSE.

