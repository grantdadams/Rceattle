---
title: "CEATLLE Walkthrough"
author: "Grant Adams"
date: "2018-08-02"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{CEATTLE Walkthrough}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: ceattle_bibliography.bib
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/dependent/fisheries-research.csl
---

CEATTLE is short for Climate-Enhanced, Age-based model with Temperature-specific Trophic Linkages and Energetics, which is a multi-species age-structured assessment model developed for groundfish in the Bering Sea, USA by Holsman et al. [-@Holsman2015]. To incorporate the impacts of climate, the model includes temperature-dependent von Bertalanffy weight-at-age functions (VBGF) and temperature-specific, bioenergetics-based predation interactions. Inputs of the model include U.S. National Marine Fisheries Service Alaska Fisheries Science Center (AFSC) survey and fishery data. Outputs include historical estimates of predation mortality, fishing mortality, biomass, recruitment, etc.

`Rceattle` is an `R` package designed to implement the CEATTLE model using Template Model Builder (`TMB`; @Kristensen2015). Rceattle is structured similar to the original manuscript in terms of modularization. Seperate functions (i.e. modules) estimate retrospective temperature- and size-specific predator rations, prey preference, and weight-at- age. These are then used as inputs to the CEATTLE model to evaluate how predation mortality, recruitment, and survival of three target species change under historical climate conditions and harvest rates. Estimates from CEATTLE in `Rceattle` can then be forward projected to derive estimates of unfished biomass and fishing mortality under future climate conditions and various harvest scenarios.  

## Equations 


**Table 1.** Population dynamics equations for species $i$ and age $j$ in each simulation year $y$. BT indicates the AFSC bottom trawl survey and EIT represents the echo-integrated acoustic- trawl survey. For all parameter definititions see Table 3.
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Definition </th>
   <th style="text-align:left;"> Equation </th>
   <th style="text-align:left;">  </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Recruitment </td>
   <td style="text-align:left;"> $N_{i1,y} = R_{i,y} = R_{0,i}e^{\tau_{i,y}}$ </td>
   <td style="text-align:left;"> T1.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Initial abundance </td>
   <td style="text-align:left;"> $N_{ij,1} = \left \{ \begin{array}{} R_{0,i}e^{ \left( -j M1_{ij} \right)}N_{0,ij} \\ R_{0,i}e^{ \left( -j M1_{ij} \right)}N_{0,i,A_i} / \left(1 - e^{\left( -j M1_{i,A_i} \right)} \right) \end{array} \right.$ </td>
   <td style="text-align:left;"> T1.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Numbers at age </td>
   <td style="text-align:left;"> $N_{i, j+1, y+1} = N_{ij,y} e^{-Z_{ij,y}}$ </td>
   <td style="text-align:left;"> T1.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $N_{i, A_i, y+1} = N_{i, A_{i}-1,y} e^{-Z_{i, A_{i}-1,y}} + N_{i, A_{i},y} e^{-Z_{i, A_{i},y}}$ </td>
   <td style="text-align:left;"> ... </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Catch </td>
   <td style="text-align:left;"> $C_{ij,y} = \frac{F_ij,y} {Z_{ij,y}} \left( 1 - e^{-Z_{ij,y}} \right) N_{ij,y}$ </td>
   <td style="text-align:left;"> T1.4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total yield (kg) </td>
   <td style="text-align:left;"> $Y_{i,y} = \sum_j^{A_i} \frac{F_ij,y} {Z_{ij,y}} \left( 1 - e^{-Z_{ij,y}} \right) N_{ij,y} W_{ij,y}$ </td>
   <td style="text-align:left;"> T1.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biomass at age (kg) </td>
   <td style="text-align:left;"> $B_{ij,y} = N_{ij,y} W_{ij,y}$ </td>
   <td style="text-align:left;"> T1.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Spawning biomass at age (kg) </td>
   <td style="text-align:left;"> $SSB_{ij,y} = B_{ij,y} \rho_{ij}$ </td>
   <td style="text-align:left;"> T1.7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total mortality at age </td>
   <td style="text-align:left;"> $Z_{ij,y} = M1_{ij} + M2_{ij,y} + F_{ij,y}$ </td>
   <td style="text-align:left;"> T1.8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fishing mortality at age </td>
   <td style="text-align:left;"> $F_{ij,y} = F_{0,i}e^{\epsilon_{i,y} s^f_{ij}}$ </td>
   <td style="text-align:left;"> T1.9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Weight at age (kg) </td>
   <td style="text-align:left;"> $W_{ij,y} = W_{\infty,iy} \left(1 - e^{\left( -K_i \left( 1 - d_{i,y} \right) \left( j - t_{0,i} \right) \right)} \right) ^ {\frac{1} {1-d_{i,y }}}$ </td>
   <td style="text-align:left;"> T1.10a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $d_{i,y} = e^{\left( \alpha_{d,i,y} + \alpha 0_{d,i} + \beta_{d,i} T_y \right)}$ </td>
   <td style="text-align:left;"> T1.10b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $W_{\infty,iy} = \left( \frac{H_i}{K_i} \right)^{1/\left(1-d_{i,y} \right)}$ </td>
   <td style="text-align:left;"> T1.10c </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BT survey biomass (kg) </td>
   <td style="text-align:left;"> $\hat{\beta}^s_{i,y} = \sum_j^{A_i} \left( N_{ij,y} e^{-0.5 Z_{ij,y}} W_{ij,y} S^s_{ij} \right)$ </td>
   <td style="text-align:left;"> T1.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EIT survey biomass (kg) </td>
   <td style="text-align:left;"> $\hat{\beta}^{eit}_{y} = \sum_j^{A_1} \left( N_{1j,y} e^{-0.5 Z_{1j,y}} W_{1j,y} S^{eit}_{1j} q_{1}^{eit}  \right)$ </td>
   <td style="text-align:left;"> T1.12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fishery age composition </td>
   <td style="text-align:left;"> $\hat{O}^f_{ij,y} = \frac{c_{ij,y}} {\sum_j C_{ij,y}}$ </td>
   <td style="text-align:left;"> T1.13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BT survey age composition </td>
   <td style="text-align:left;"> $\hat{O}^{s}_{ij,y} = \frac{N_{ij,y} e^{-0.5 Z_{ij,y}}  s^s_{ij}} {\sum_j \left(N_{ij,y} e^{-0.5 Z_{ij,y}}  s^s_{ij} \right)}$ </td>
   <td style="text-align:left;"> T1.14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EIT survey age composition </td>
   <td style="text-align:left;"> $\hat{O}^{eit}_{1j,y} = \frac{N_{1j,y} e^{-0.5 Z_{1j,y}}  s^{eit}_{1j} q^{eit}_{1}} {\sum_j \left( N_{1j,y} e^{-0.5 Z_{1j,y}}  s^{eit}_{1j} q^{eit}_{1} \right)}$ </td>
   <td style="text-align:left;"> T1.15 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BT selectivity </td>
   <td style="text-align:left;"> $s_{ij}^s = \frac{1} {1 + e^{\left(-b_i^s * j-a_i^s \right)}}$ </td>
   <td style="text-align:left;"> T1.16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fishery selectivity </td>
   <td style="text-align:left;"> $s^f_{ij} = \left \{ \begin{array}{} e^{\eta_ij} ~~ j \leq A_{\eta , i} \\ e^{\eta_{i,A_{\eta,i}}} ~~ j \gt A_{\eta,i} \end{array} \right.$ </td>
   <td style="text-align:left;"> T1.17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proportion females </td>
   <td style="text-align:left;"> $\omega_{ij} = \frac{e^{-j M_{fem}}} { e^{-j M_{fem}} + e^{j M_{male}}}$ </td>
   <td style="text-align:left;"> T1.18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proportion of mature females </td>
   <td style="text-align:left;"> $\rho_{ij} = \omega_{ij} \phi_{ij}$ </td>
   <td style="text-align:left;"> T1.19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Weight at age (kg) </td>
   <td style="text-align:left;"> $W_{ij,y} = W^{fem}_{ij,y} \omega_{ij} + \left(1 - \omega_{ij} \right) W^{male}_{ij,y}$ </td>
   <td style="text-align:left;"> T1.20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residual natural mortality </td>
   <td style="text-align:left;"> $M1_{ij} = M^{fem}_{ij} \omega_{ij} + \left(1 - \omega_{ij} \right) M^{male}_{ij}$ </td>
   <td style="text-align:left;"> T1.21 </td>
  </tr>
</tbody>
</table>
  
  
**Table 2.** Predation mortality equations for predators $p$ of age $a$, and prey $i$ of age $j$
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Definition </th>
   <th style="text-align:left;"> Equation </th>
   <th style="text-align:left;">  </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Predation mortality </td>
   <td style="text-align:left;"> $M2_{ij,y} = \sum_{pa} \left( \frac{N_{pa,y} \delta_{pa,y} {S}_{paij}} {\sum_{ij} \left({S}_{paij} B_{ij,y} \right) + B_p^{other} \left(1 - \sum_{ij} \left({S}_{paij} \right) \right)} \right)$ </td>
   <td style="text-align:left;"> T2.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Predator-prey suitability </td>
   <td style="text-align:left;"> $\hat{S}_{paij} = \frac{1}{n_y} \sum_y \left( \frac{ \frac{\bar{U}_{paij}} {B_{ij,y}} } {\sum_{ij} \left( \frac{\bar{U}_{paij}} {B_{ij,y}} \right) + \frac{1 + \sum_{ij} \bar{U}_{paij}} {B_p^{other}} } \right)$ </td>
   <td style="text-align:left;"> T2.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mean gravimetric diet proportion </td>
   <td style="text-align:left;"> $\bar{U}_{paij} = \frac{\sum_y {U_{paij,y}}} {n_y}$ </td>
   <td style="text-align:left;"> T2.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Individual specific ration ($kg~kg^{-1}~yr^{-1}$) </td>
   <td style="text-align:left;"> $\delta_{pa,y} = \varphi_p \alpha_{\delta} W_{pa,y}^{\left(1+\beta_{\delta}\right)} f\left(T_y \right)_p$ </td>
   <td style="text-align:left;"> T2.4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Temperature scaling algorithim </td>
   <td style="text-align:left;"> $f \left(T_y \right)_p = V^X e^{\left(X \left(1-V \right) \right)}$ </td>
   <td style="text-align:left;"> T2.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $V = \left( T^{cm}_p - T_y \right)/ \left( T^{cm}_p - T^{co}_p \right)$ </td>
   <td style="text-align:left;"> T2.5a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $X = \left(Z^2 \left(1+ \left( 1 + 40 /Y \right)^{0.5} \right)^2 \right)/400$ </td>
   <td style="text-align:left;"> T2.5b </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $Z = ln\left(Q^c_p \right) \left( T^{cm}_p - T^{co}_p \right)$ </td>
   <td style="text-align:left;"> T2.5c </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ... </td>
   <td style="text-align:left;"> $Y = ln\left(Q^c_p \right) \left( T^{cm}_p - T^{co}_p + 2\right)$ </td>
   <td style="text-align:left;"> T2.5d </td>
  </tr>
</tbody>
</table>
  
  
**Table 3.** Parameter definitions.  
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Parameter </th>
   <th style="text-align:left;"> Definition </th>
   <th style="text-align:left;"> Type </th>
   <th style="text-align:left;"> Model Object </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Year </td>
   <td style="text-align:left;"> $y$ </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> i </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Predator </td>
   <td style="text-align:left;"> $p$ </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Predator age (years) </td>
   <td style="text-align:left;"> $a$ </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Prey </td>
   <td style="text-align:left;"> $i$ </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> k </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Prey age (years) </td>
   <td style="text-align:left;"> $j$ </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of prey species </td>
   <td style="text-align:left;"> $n_i$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nspp </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of predator species </td>
   <td style="text-align:left;"> $n_p$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of prey ages </td>
   <td style="text-align:left;"> $A_i$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nages </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of predator ages </td>
   <td style="text-align:left;"> $A_p$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of simulation years </td>
   <td style="text-align:left;"> $n_y$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Start year </td>
   <td style="text-align:left;"> $y_0$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> styr </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Annual relative foraging rate ($d~yr^{-1}$) </td>
   <td style="text-align:left;"> $\hat{\varphi}_p$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Intercept of the allometic maximum consumption function ($g~g^{-1}~yr^{-1}$) </td>
   <td style="text-align:left;"> $\alpha_{\delta}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> aLW </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Allometric slope of maximum consumption </td>
   <td style="text-align:left;"> $\beta_{\delta}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> bLW </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Consumption maximum physiological temperature (°C) </td>
   <td style="text-align:left;"> $T^{cm}_p$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> Tcm </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Consumption optimum physiological temperature (°C) </td>
   <td style="text-align:left;"> $T^{co}_p$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> Tco </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Max consumption parameter </td>
   <td style="text-align:left;"> $Q^c_p$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> Qc </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mean recruitment </td>
   <td style="text-align:left;"> $R_{0,i}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Annual recruitment deviation </td>
   <td style="text-align:left;"> $\tau_{i,y}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;"> rec_dev </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Initial abundance </td>
   <td style="text-align:left;"> $N_{0,ij}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mean fishing mortality </td>
   <td style="text-align:left;"> $F_{0,i}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Anuual fishing mortality deviation </td>
   <td style="text-align:left;"> $\epsilon_{i,y}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fishery age selectivity coefficient </td>
   <td style="text-align:left;"> $\eta_{ij}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Survey age selectivity slope </td>
   <td style="text-align:left;"> $b^s_i$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Survey age selectivity limit </td>
   <td style="text-align:left;"> $a^s_i$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VBGF allometric slope of consumption </td>
   <td style="text-align:left;"> $d_{i,y}$ </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> d </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VBGF max asymptotic weight (kg) </td>
   <td style="text-align:left;"> $W_{\infty,iy}$ </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> Winf </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proportion of mature females at age </td>
   <td style="text-align:left;"> $\rho_{ij}$ </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residual natural mortality </td>
   <td style="text-align:left;"> $M1_{ij}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> M1_base </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Intercept for VBGF $d$ parameter </td>
   <td style="text-align:left;"> $\alpha0_{d,i}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Annual intercept for VBGF $d$ parameter </td>
   <td style="text-align:left;"> $\alpha_{d,i,y}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> log_mean_d </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Temperature covariate for VBGF $d$ parameter </td>
   <td style="text-align:left;"> $\beta_{d,i}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> Tcoef </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VBGF energy loss constant ($kg~kg^{-1}~yr^{-1}$) </td>
   <td style="text-align:left;"> $K_i$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> logK </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VBGF assimilation constant ($kg~kg^{-1}~yr^{-1}$) </td>
   <td style="text-align:left;"> $H_i$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> logH </td>
  </tr>
  <tr>
   <td style="text-align:left;"> VBGF age when $W_{ij,y} = 0$ (years) </td>
   <td style="text-align:left;"> $t_{0,i}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> t0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EIT survey selectivity </td>
   <td style="text-align:left;"> $S^{eit}_{1j}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Female natural mortality </td>
   <td style="text-align:left;"> $M_{i}^{fem}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Male natural mortality </td>
   <td style="text-align:left;"> $M_{i}^{male}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Female proportion of population </td>
   <td style="text-align:left;"> $\omega_{ij}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age-specific maturity proportions </td>
   <td style="text-align:left;"> $\varphi_{ij}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> pmature </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed total yield ($kg$) </td>
   <td style="text-align:left;"> $C^*_{i,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> tc_biom_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed fishery age comp. </td>
   <td style="text-align:left;"> $O^f_{ij,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> fsh_age_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed BT age comp. </td>
   <td style="text-align:left;"> $O^{s}_{ij,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> srv_age_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed EIT age comp. </td>
   <td style="text-align:left;"> $O^{eit}_{ij,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> obs_eit_age </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed BT survey biomass (kg) </td>
   <td style="text-align:left;"> $\beta^s_{i,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> srv_bio </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed EIT survey biomass (kg) </td>
   <td style="text-align:left;"> $\beta^{eit}_y$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> obs_eit </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Bottom temperature (°C) </td>
   <td style="text-align:left;"> $T_y$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> TempC </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gravimetric proportion of prey in predator stomach </td>
   <td style="text-align:left;"> $U_{paij,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biomass of other prey (kg) </td>
   <td style="text-align:left;"> $B^{other}_p$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> other_food </td>
  </tr>
  <tr>
   <td style="text-align:left;"> **Not in table 3** </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Annual survey biomass error </td>
   <td style="text-align:left;"> $CV_{s,i,y}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> srv_Mean_CV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years with total observed catch </td>
   <td style="text-align:left;"> $n_{Y,i}$ </td>
   <td style="text-align:left;"> M </td>
   <td style="text-align:left;"> nyrs_tc_biom_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years with total observed catch </td>
   <td style="text-align:left;"> $y_{Y,i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_tc_biom_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years in the fishery sp_age composition data </td>
   <td style="text-align:left;"> $n_{y,i}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;"> nyrs_fsh_comp </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of estimation years </td>
   <td style="text-align:left;"> $n_{y,est}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs_est </td>
  </tr>
  <tr>
   <td style="text-align:left;"> End year </td>
   <td style="text-align:left;"> $y_{n_y}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> endyr </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years in the fishery age composition data </td>
   <td style="text-align:left;"> $n_{y_{O^f},i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs_fsh_comp </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years in the fishery age composition data </td>
   <td style="text-align:left;"> y_{O^f,i} </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_fsh_comp </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Method of calculating fishery age </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> fsh_age_type </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of fishery age bins </td>
   <td style="text-align:left;"> $n_{f,bin}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> fsh_age_bins </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years with weight-at-age data </td>
   <td style="text-align:left;"> $n_{y,W,i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs_wt_at_age </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years with weight-at-age data </td>
   <td style="text-align:left;"> $y_{W,i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_wt_at_age </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Weight-at-age data </td>
   <td style="text-align:left;"> $W_i$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> wt </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years in the BT survey data </td>
   <td style="text-align:left;"> $n_{y_{\beta^s},i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs_srv_biom </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years in the BT survey data </td>
   <td style="text-align:left;"> $y_{\beta^s,i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_srv_biom </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BT survey standard error </td>
   <td style="text-align:left;"> $\sigma_{s,i}$ </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> srv_biom_se </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years in the BT survey age or length composition data </td>
   <td style="text-align:left;"> $n_{y_{O^s},i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs_srv_age </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years in the BT  survey age composition data </td>
   <td style="text-align:left;"> $y_{O^s,i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_srv_age </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Method of calculating BT survey age type (age or length) </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> srv_age_type </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of BT survey age bins </td>
   <td style="text-align:left;"> $n_{s,bin, i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> srv_age_bins </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample size for BT survey age composition multinomial </td>
   <td style="text-align:left;"> $n_{O^s, i}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> srv_age_n </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed survey BT size compositions </td>
   <td style="text-align:left;"> $O^{s}_{ij}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> srv_age_sizes </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age transition matrix </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> age_trans_matrix </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years in the EIT survey data </td>
   <td style="text-align:left;"> $n_{y_{\beta^{eit}}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> n_eit </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years in the BT survey data </td>
   <td style="text-align:left;"> $y_{\beta^{eit}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_eit </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample size for EIT survey age composition multinomial </td>
   <td style="text-align:left;"> $n_{O^{eit}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> eit_age_n </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Number of years in the EIT selectivtiy data </td>
   <td style="text-align:left;"> $n_{y_{S^{eit}}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> nyrs_eit_sel </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Years in the BT selectivity data </td>
   <td style="text-align:left;"> $y_{S^{eit}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> yrs_eit_sel </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sample size for EIT survey age composition multinomial </td>
   <td style="text-align:left;"> $n_{O^{eit}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> eit_sel </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Sex specific mortality and weight-at-age: 1 for combined, 2: for seperate </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> mf_type </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Proportion </td>
   <td style="text-align:left;"> $n_{O^{eit}}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> propMorF </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed catch-at-age </td>
   <td style="text-align:left;"> $C_{ij,y}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> obs_catch </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Estimated catch-at-age </td>
   <td style="text-align:left;"> $\hat{C}_{ij,y}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;"> obs_catch_hat </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed total catch </td>
   <td style="text-align:left;"> $C_{i,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> tc_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Estimated total catch </td>
   <td style="text-align:left;"> $\hat{C}_{i,y}$ </td>
   <td style="text-align:left;"> D </td>
   <td style="text-align:left;"> tc_hat </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Estimated total yield </td>
   <td style="text-align:left;"> $B_{ij,y} = N_{ij,y} W_{ij,y}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;"> tc_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Observed fishery age composition </td>
   <td style="text-align:left;"> ${O}^f_{ij,y}$ </td>
   <td style="text-align:left;"> I </td>
   <td style="text-align:left;"> fsh_age_obs </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Estimated fishery age composition </td>
   <td style="text-align:left;"> $\hat{O}^f_{ij,y}$ </td>
   <td style="text-align:left;"> E </td>
   <td style="text-align:left;"> fsh_age_hat </td>
  </tr>
</tbody>
</table>
  
  
**Table 4.** Components of the likelihood function for each species $i$ of age $j$ in year $y$.  
<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Description </th>
   <th style="text-align:left;"> Equation </th>
   <th style="text-align:left;"> Data source </th>
   <th style="text-align:left;">  </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> **Data components** </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BT survey biomass </td>
   <td style="text-align:left;"> $\sum_i \sum_y \frac{\left[ ln \left(\beta^s_{i,y} \right) - ln \left(\hat{\beta}^s_{i,y} \right) \right]^2} {2 \sigma^2_{s,i}}$ </td>
   <td style="text-align:left;"> NMFS annual EBS BT survey </td>
   <td style="text-align:left;"> T4.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BT survey age composition </td>
   <td style="text-align:left;"> $-\sum_i n_i \sum_y \sum_j \left( O^s_{ij,y} + v \right) ~ ln \left(\hat{O}^s_{ij,y} + v \right)$ </td>
   <td style="text-align:left;"> NMFS annual EBS BT survey </td>
   <td style="text-align:left;"> T4.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EIT survey biomass </td>
   <td style="text-align:left;"> $\sum_i \sum_y \frac{\left[ ln \left(\beta^{eit}_{y} \right) - ln \left(\hat{\beta}^{eit}_{y} \right) \right]^2} {2 \sigma^2_{eit}},~\sigma_{eit}=0.2$ </td>
   <td style="text-align:left;"> Pollock acoustic trawl survey </td>
   <td style="text-align:left;"> T4.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EIT survey age composition </td>
   <td style="text-align:left;"> $- n \sum_y \sum_j \left( O^{eit}_{1j,y} + v \right) ~ ln \left(\hat{O}^{eit}_{1j,y} + v \right)$ </td>
   <td style="text-align:left;"> Pollock acoustic trawl survey </td>
   <td style="text-align:left;"> T4.4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total catch </td>
   <td style="text-align:left;"> $\sum_i \sum_y \frac{\left[ ln \left(C^*_{i,y} \right) - ln \left(\hat{C}^*_{i,y} \right) \right]^2} {2 \sigma^2_{C}},~\sigma_C = 0.05$ </td>
   <td style="text-align:left;"> Fishery observer data </td>
   <td style="text-align:left;"> T4.5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fishery age composition </td>
   <td style="text-align:left;"> $-\sum_i n_i \sum_y \sum_j \left( O^f_{ij,y} + v \right) ~ ln \left(\hat{O}^f_{ij,y} + v \right)$ </td>
   <td style="text-align:left;"> Fishery observer data </td>
   <td style="text-align:left;"> T4.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> **Penalties** </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fishery selectivity </td>
   <td style="text-align:left;"> $\sum_i \sum_j^{A_i - 1} \chi \left[ ln \left(\frac{\eta^f_{ij}} {\eta^f_{ij+1}} \right) - ln \left(\frac{\eta^f_{ij+1}} {\eta^f_{ij+2}} \right) \right]^2 , ~ \chi = \left \{  \begin{array}{} 20,~if \eta^f_{ij} \gt \eta^f_{ij+1} \\ 0,~if \eta^f_{ij} \leq \eta^f_{ij+1} \end{array} \right.$ </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> T4.7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> **Priors** </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> $\sum_i \sum_j \left( \tau_{i,y} \right)^2$ </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> T4.8 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> $\sum_i \sum_j \left( N_{0,ij} \right)^2$ </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> T4.9 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> $\sum_i \sum_j \left( \varepsilon_{i,y} \right)^2$ </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> T4.10 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> $v = 0.001$ </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;">  </td>
  </tr>
</tbody>
</table>

## Errata:
Table 3: Female natural mortality is missing an F under the "Type" column.
Table 4 legend: "function" is misspelled.
Table 4. T.4.6 C is missing the hat in the likelihood equation.
Table 4. T4.7 the upper value of chi, $\eta^f_{ji+1}$ should be $\eta^f_{ij+1}$
Table 4. T4.2 Data source NFMS

## Thoughts
Table 3. Change C* to Y and in table 2 Y to hat{Y}.

## Questions
What is propMorF?

# References
