---
title: "Rceattle"
author: "Grant Adams, Kirstin Holsman"
date: "8/2/2018"
output:
  pdf_document: default
  word_document: default
csl: https://raw.githubusercontent.com/citation-style-language/styles/master/dependent/fisheries-research.csl
bibliography: ceattle_bibliography.bib
vignette: |
  %\VignetteIndexEntry{CEATTLE Walkthrough} %\VignetteEngine{knitr::rmarkdown} %\VignetteEncoding{UTF-8}
---

CEATTLE is short for Climate-Enhanced, Age-based model with Temperature-specific Trophic Linkages and Energetics, which is a multi-species age-structured assessment model developed for groundfish in the Bering Sea, USA by Holsman et al. [-@Holsman2015]. To incorporate the impacts of climate, the model includes temperature-dependent von Bertalanffy weight-at-age functions (VBGF) and temperature-specific, bioenergetics-based predation interactions. Inputs of the model include U.S. National Marine Fisheries Service Alaska Fisheries Science Center (AFSC) survey and fishery data. Outputs include historical estimates of predation mortality, fishing mortality, biomass, recruitment, etc.

`Rceattle` is an `R` package designed to implement the CEATTLE model using Template Model Builder (`TMB`; @Kristensen2015). Rceattle is structured similar to the original manuscript in terms of modularization. Seperate functions (i.e. modules) estimate retrospective temperature- and size-specific predator rations, prey preference, and weight-at- age. These are then used as inputs to the CEATTLE model to evaluate how predation mortality, recruitment, and survival of three target species change under historical climate conditions and harvest rates. Estimates from CEATTLE in `Rceattle` can then be forward projected to derive estimates of unfished biomass and fishing mortality under future climate conditions and various harvest scenarios.

Currently `RCeattle` takes CEATTLE ADMB based inputs and estimates them in `TMB`. To run `Rceattle` a user selects a control file located in a data directory with \code{.dat} files used in ADMB. The model will then initialize assuming $0$ for all parameters unless a parameter list is provided or `"ceattle.par"` or `"ceattle_par.std"` is selected. In the later case, a .par or .std file provided from a previously estimated ADMB model must by included in the folder prior to the data directory. Currently, .dat files and ADMB outputs are included in the `\data` folder.

For example, to run the 2017 single species assessment, the following can be specified:

```r
ss_run_no_re <- Rceattle(data_list = NULL, ctlFilename = "asmnt2017_0A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_SS_Files/dat files/", inits = NULL, debug = FALSE, random_rec = FALSE, niter = 3, file_name = "Report/2017_single_species_assessment")
```

For the 2017 multispecies model starting from ADMB parameters, the following can be specified:

```r
ms_run_no_re <- Rceattle(data_list = NULL, ctlFilename = "asmnt2017_2A_corrected", TMBfilename = "CEATTLE_BSAI_MS_v01", dat_dir =  "data/BS_MS_Files/dat files/", inits = "ceattle.par", debug = FALSE, random_rec = FALSE, niter = 3, file_name = "Report/2017_multi_species_assessment")
```

## Equations 


**Table 1.** Population dynamics equations for species $i$ and age $j$ in each simulation year $y$. BT indicates the AFSC bottom trawl survey and EIT represents the echo-integrated acoustic- trawl survey. For all parameter definititions see Table 3.

Definition                     Equation                                                                                                                                                                                                                  
-----------------------------  -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  -------
Recruitment                    $N_{i1,y} = R_{i,y} = R_{0,i }e^{\tau_{i,y}}$                                                                                                                                                                      T1.1   
Initial abundance              $N_{ij,1} = \left \{ \begin{array}{} R_{0,i}e^{ \left( -j M1_{ij} \right)}N_{0,ij} \\ R_{0,i}e^{ \left( -j M1_{ij} \right)}N_{0,i,A_i} / \left(1 - e^{\left( -j M1_{i,A_i} \right)} \right) \end{array} \right.$   T1.2   
Numbers at age                 $N_{i, j+1, y+1} = N_{ij,y} e^{-Z_{ij,y}}$                                                                                                                                                                         T1.3   
...                            $N_{i, A_i, y+1} = N_{i, A_{i}-1,y} e^{-Z_{i, A_{i}-1,y}} + N_{i, A_{i},y} e^{-Z_{i, A_{i},y}}$                                                                                                                    ...    
Catch                          $C_{ij,y} = \frac{F_ij,y} {Z_{ij,y}} \left( 1 - e^{-Z_{ij,y}} \right) N_{ij,y}$                                                                                                                                    T1.4   
Total yield (kg)               $Y_{i,y} = \sum_j^{A_i} \frac{F_ij,y} {Z_{ij,y}} \left( 1 - e^{-Z_{ij,y}} \right) N_{ij,y} W_{ij,y}$                                                                                                               T1.5   
Biomass at age (kg)            $B_{ij,y} = N_{ij,y} W_{ij,y}$                                                                                                                                                                                     T1.6   
Spawning biomass at age (kg)   $SSB_{ij,y} = B_{ij,y} \rho_{ij}$                                                                                                                                                                                  T1.7   
Total mortality at age         $Z_{ij,y} = M1_{ij} + M2_{ij,y} + F_{ij,y}$                                                                                                                                                                        T1.8   
Fishing mortality at age       $F_{ij,y} = F_{0,i}e^{\epsilon_{i,y} s^f_{ij}}$                                                                                                                                                                    T1.9   
Weight at age (kg)             $W_{ij,y} = W_{\infty,iy} \left(1 - e^{\left( -K_i \left( 1 - d_{i,y} \right) \left( j - t_{0,i} \right) \right)} \right) ^ {\frac{1} {1-d_{i,y }}}$                                                               T1.10a 
...                            $d_{i,y} = e^{\left( \alpha_{d,i,y} + \alpha 0_{d,i} + \beta_{d,i} T_y \right)}$                                                                                                                                   T1.10b 
...                            $W_{\infty,iy} = \left( \frac{H_i}{K_i} \right)^{1/\left(1-d_{i,y} \right)}$                                                                                                                                       T1.10c 
BT survey biomass (kg)         $\hat{\beta}^s_{i,y} = \sum_j^{A_i} \left( N_{ij,y} e^{-0.5 Z_{ij,y}} W_{ij,y} S^s_{ij} \right)$                                                                                                                   T1.11  
EIT survey biomass (kg)        $\hat{\beta}^{eit}_{y} = \sum_j^{A_1} \left( N_{1j,y} e^{-0.5 Z_{1j,y}} W_{1j,y} S^{eit}_{1j} q_{1}^{eit}  \right)$                                                                                                T1.12  
Fishery age composition        $\hat{O}^f_{ij,y} = \frac{c_{ij,y}} {\sum_j C_{ij,y}}$                                                                                                                                                             T1.13  
BT survey age composition      $\hat{O}^{s}_{ij,y} = \frac{N_{ij,y} e^{-0.5 Z_{ij,y}}  s^s_{ij}} {\sum_j \left(N_{ij,y} e^{-0.5 Z_{ij,y}}  s^s_{ij} \right)}$                                                                                     T1.14  
EIT survey age composition     $\hat{O}^{eit}_{1j,y} = \frac{N_{1j,y} e^{-0.5 Z_{1j,y}}  s^{eit}_{1j} q^{eit}_{1}} {\sum_j \left( N_{1j,y} e^{-0.5 Z_{1j,y}}  s^{eit}_{1j} q^{eit}_{1} \right)}$                                                  T1.15  
BT selectivity                 $s_{ij}^s = \frac{1} {1 + e^{\left(-b_i^s * j-a_i^s \right)}}$                                                                                                                                                     T1.16  
Fishery selectivity            $s^f_{ij} = \left \{ \begin{array}{} e^{\eta_ij} ~~ j \leq A_{ \eta , i} \\ e^{\eta_{i,A_{ \eta ,i}}} ~~ j \gt A_{ \eta,i} \end{array} \right.$                                                                    T1.17  
Proportion females             $\omega_{ij} = \frac{e^{-j M_{fem}}} { e^{-j M_{fem}} + e^{j M_{male}}}$                                                                                                                                           T1.18  
Proportion of mature females   $\rho_{ij} = \omega_{ij} \phi_{ij}$                                                                                                                                                                                T1.19  
Weight at age (kg)             $W_{ij,y} = W^{fem}_{ij,y} \omega_{ij} + \left(1 - \omega_{ij} \right) W^{male}_{ij,y}$                                                                                                                            T1.20  
Residual natural mortality     $M1_{ij} = M^{fem}_{ij} \omega_{ij} + \left(1 - \omega_{ij} \right) M^{male}_{ij}$                                                                                                                                 T1.21  
  
  
**Table 2.** Predation mortality equations for predators $p$ of age $a$, and prey $i$ of age $j$

Definition                                          Equation                                                                                                                                                                                                                  
--------------------------------------------------  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  ------
Predation mortality                                 $M2_{ij,y} = \sum_{pa} \left( \frac{N_{pa,y} \delta_{pa,y} {S}_{paij}} {\sum_{ij} \left({S}_{paij} B_{ij,y} \right) + B_p^{other} \left(1 - \sum_{ij} \left({S}_{paij} \right) \right)} \right)$                    T2.1  
Predator-prey suitability                           $\hat{S}_{paij} = \frac{1}{n_y} \sum_y \left( \frac{ \frac{\bar{U}_{paij}} {B_{ij,y}} } {\sum_{ij} \left( \frac{\bar{U}_{paij}} {B_{ij,y}} \right) + \frac{1 + \sum_{ij} \bar{U}_{paij}} {B_p^{other}} } \right)$   T2.2  
Mean gravimetric diet proportion                    $\bar{U}_{paij} = \frac{\sum_y {U_{paij,y}}} {n_y}$                                                                                                                                                                 T2.3  
Individual specific ration ($kg~kg^{-1}~yr^{-1}$)   $\delta_{pa,y} = \varphi_p \alpha_{\delta} W_{pa,y}^{\left(1+\beta_{\delta}\right)} f\left(T_y \right)_p$                                                                                                           T2.4  
Temperature scaling algorithim                      $f \left(T_y \right)_p = V^X e^{\left(X \left(1-V \right) \right)}$                                                                                                                                                 T2.5  
...                                                 $V = \left( T^{cm}_p - T_y \right)/ \left( T^{cm}_p - T^{co}_p \right)$                                                                                                                                             T2.5a 
...                                                 $X = \left(Z^2 \left(1+ \left( 1 + 40 /Y \right)^{0.5} \right)^2 \right)/400$                                                                                                                                       T2.5b 
...                                                 $Z = ln\left(Q^c_p \right) \left( T^{cm}_p - T^{co}_p \right)$                                                                                                                                                      T2.5c 
...                                                 $Y = ln\left(Q^c_p \right) \left( T^{cm}_p - T^{co}_p + 2\right)$                                                                                                                                                   T2.5d 
  
  
**Table 3.** Parameter definitions.  

Parameter                                                                      Definition                       Type   Model Object     
-----------------------------------------------------------------------------  -------------------------------  -----  -----------------
Year                                                                           $y$                              M      i                
Predator                                                                       $p$                              M                       
Predator age (years)                                                           $a$                              M                       
Prey                                                                           $i$                              M      k                
Prey age (years)                                                               $j$                              M                       
Number of prey species                                                         $n_i$                            I      nspp             
Number of predator species                                                     $n_p$                            I                       
Number of prey ages                                                            $A_i$                            I      nages            
Number of predator ages                                                        $A_p$                            I                       
Number of simulation years                                                     $n_y$                            I      nyrs             
Start year                                                                     $y_0$                            I      styr             
Annual relative foraging rate ($d~yr^{-1}$)                                    $\hat{\varphi}_p$                I                       
Intercept of the allometic maximum consumption function ($g~g^{-1}~yr^{-1}$)   $\alpha_{\delta}$                I      aLW              
Allometric slope of maximum consumption                                        $\beta_{\delta}$                 I      bLW              
Consumption maximum physiological temperature (°C)                             $T^{cm}_p$                       I      Tcm              
Consumption optimum physiological temperature (°C)                             $T^{co}_p$                       I      Tco              
Max consumption parameter                                                      $Q^c_p$                          I      Qc               
Mean recruitment                                                               $R_{0,i}$                        E                       
Annual recruitment deviation                                                   $\tau_{i,y}$                     E      rec_dev          
Initial abundance                                                              $N_{0,ij}$                       E                       
Mean fishing mortality                                                         $F_{0,i}$                        E                       
Anuual fishing mortality deviation                                             $\epsilon_{i,y}$                 E                       
Fishery age selectivity coefficient                                            $\eta_{ij}$                      E                       
Survey age selectivity slope                                                   $b^s_i$                          E                       
Survey age selectivity limit                                                   $a^s_i$                          E                       
VBGF allometric slope of consumption                                           $d_{i,y}$                        P      d                
VBGF max asymptotic weight (kg)                                                $W_{\infty,iy}$                  P      Winf             
Proportion of mature females at age                                            $\rho_{ij}$                      P                       
Residual natural mortality                                                     $M1_{ij}$                        F      M1_base          
Intercept for VBGF $d$ parameter                                               $\alpha0_{d,i}$                  F                       
Annual intercept for VBGF $d$ parameter                                        $\alpha_{d,i,y}$                 F      log_mean_d       
Temperature covariate for VBGF $d$ parameter                                   $\beta_{d,i}$                    F      Tcoef            
VBGF energy loss constant ($kg~kg^{-1}~yr^{-1}$)                               $K_i$                            F      logK             
VBGF assimilation constant ($kg~kg^{-1}~yr^{-1}$)                              $H_i$                            F      logH             
VBGF age when $W_{ij,y} = 0$ (years)                                           $t_{0,i}$                        F      t0               
EIT survey selectivity                                                         $S^{eit}_{1j}$                   F                       
Female natural mortality                                                       $M_{i}^{fem}$                    F                       
Male natural mortality                                                         $M_{i}^{male}$                   F                       
Female proportion of population                                                $\omega_{ij}$                    F                       
Age-specific maturity proportions                                              $\varphi_{ij}$                   F      pmature          
Observed total yield ($kg$)                                                    $C^*_{i,y}$                      D      tc_biom_obs      
Observed fishery age comp.                                                     $O^f_{ij,y}$                     D      fsh_age_obs      
Observed BT age comp.                                                          $O^{s}_{ij,y}$                   D      srv_age_obs      
Observed EIT age comp.                                                         $O^{eit}_{ij,y}$                 D      obs_eit_age      
Observed BT survey biomass (kg)                                                $\beta^s_{i,y}$                  D      srv_bio          
Observed EIT survey biomass (kg)                                               $\beta^{eit}_y$                  D      obs_eit          
Bottom temperature (°C)                                                        $T_y$                            D      TempC            
Gravimetric proportion of prey in predator stomach                             $U_{paij,y}$                     D                       
Biomass of other prey (kg)                                                     $B^{other}_p$                    D      other_food       
**Not in table 3**                                                                                                                      
Annual survey biomass error                                                    $CV_{s,i,y}$                     F      srv_Mean_CV      
Number of years with total observed catch                                      $n_{Y,i}$                        M      nyrs_tc_biom_obs 
Years with total observed catch                                                $y_{Y,i}$                        I      yrs_tc_biom_obs  
Number of years in the fishery sp_age composition data                         $n_{y,i}$                        E      nyrs_fsh_comp    
Number of estimation years                                                     $n_{y,est}$                      I      nyrs_est         
End year                                                                       $y_{n_y}$                        I      endyr            
Number of years in the fishery age composition data                            $n_{y_{O^f},i}$                  I      nyrs_fsh_comp    
Years in the fishery age composition data                                      $y_{O^f,i}$                      I      yrs_fsh_comp     
Method of calculating fishery age                                                                               I      fsh_age_type     
Number of fishery age bins                                                     $n_{f,bin}$                      I      fsh_age_bins     
Number of years with weight-at-age data                                        $n_{y,W,i}$                      I      nyrs_wt_at_age   
Years with weight-at-age data                                                  $y_{W,i}$                        I      yrs_wt_at_age    
Weight-at-age data                                                             $W_i$                            D      wt               
Number of years in the BT survey data                                          $n_{y_{\beta^s},i}$              I      nyrs_srv_biom    
Years in the BT survey data                                                    $y_{\beta^s,i}$                  I      yrs_srv_biom     
BT survey standard error                                                       $\sigma_{s,i}$                   F      srv_biom_se      
Number of years in the BT survey age or length composition data                $n_{y_{O^s},i}$                  I      nyrs_srv_age     
Years in the BT  survey age composition data                                   $y_{O^s,i}$                      I      yrs_srv_age      
Method of calculating BT survey age type (age or length)                                                        I      srv_age_type     
Number of BT survey age bins                                                   $n_{s,bin, i}$                   I      srv_age_bins     
Sample size for BT survey age composition multinomial                          $n_{O^s, i}$                     I      srv_age_n        
Observed survey BT size compositions                                           $O^{s}_{ij}$                     I      srv_age_sizes    
Age transition matrix                                                                                           I      age_trans_matrix 
Number of years in the EIT survey data                                         $n_{y_{\beta^{eit}}}$            I      n_eit            
Years in the BT survey data                                                    $y_{\beta^{eit}}$                I      yrs_eit          
Sample size for EIT survey age composition multinomial                         $n_{O^{eit}}$                    I      eit_age_n        
Number of years in the EIT selectivtiy data                                    $n_{y_{S^{eit}}}$                I      nyrs_eit_sel     
Years in the BT selectivity data                                               $y_{S^{eit}}$                    I      yrs_eit_sel      
Sample size for EIT survey age composition multinomial                         $n_{O^{eit}}$                    I      eit_sel          
Sex specific mortality and weight-at-age: 1 for combined, 2: for seperate                                       I      mf_type          
Proportion                                                                     $n_{O^{eit}}$                    I      propMorF         
Observed catch-at-age                                                          $C_{ij,y}$                       I      obs_catch        
Estimated catch-at-age                                                         $\hat{C}_{ij,y}$                 E      obs_catch_hat    
Observed total catch                                                           $C_{i,y}$                        D      tc_obs           
Estimated total catch                                                          $\hat{C}_{i,y}$                  D      tc_hat           
Estimated total yield                                                          $B_{ij,y} = N_{ij,y} W_{ij,y}$   E      tc_obs           
Observed fishery age composition                                               ${O}^f_{ij,y}$                   I      fsh_age_obs      
Estimated fishery age composition                                              $\hat{O}^f_{ij,y}$               E      fsh_age_hat      
  
  
**Table 4.** Components of the likelihood function for each species $i$ of age $j$ in year $y$.  

Description                  Equation                                                                                                                                                                                                                                                                                            Data source                           
---------------------------  --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  ------------------------------  ------
**Data components**                                                                                                                                                                                                                                                                                                                                                    
BT survey biomass            $\sum_i \sum_y \frac{\left[ ln \left(\beta^s_{i,y} \right) - ln \left(\hat{\beta}^s_{i,y} \right) \right]^2} {2 \sigma^2_{s,i}}$                                                                                                                                                                    NMFS annual EBS BT survey       T4.1  
BT survey age composition    $-\sum_i n_i \sum_y \sum_j \left( O^s_{ij,y} + v \right) ~ ln \left(\hat{O}^s_{ij,y} + v \right)$                                                                                                                                                                                                   NMFS annual EBS BT survey       T4.2  
EIT survey biomass           $\sum_i \sum_y \frac{\left[ ln \left(\beta^{eit}_{y} \right) - ln \left(\hat{\beta}^{eit}_{y} \right) \right]^2} {2 \sigma^2_{eit}},~\sigma_{eit}=0.2$                                                                                                                                              Pollock acoustic trawl survey   T4.3  
EIT survey age composition   $- n \sum_y \sum_j \left( O^{eit}_{1j,y} + v \right) ~ ln \left(\hat{O}^{eit}_{1j,y} + v \right)$                                                                                                                                                                                                   Pollock acoustic trawl survey   T4.4  
Total catch                  $\sum_i \sum_y \frac{\left[ ln \left(C^*_{i,y} \right) - ln \left(\hat{C}^*_{i,y} \right) \right]^2} {2 \sigma^2_{C}},~\sigma_C = 0.05$                                                                                                                                                             Fishery observer data           T4.5  
Fishery age composition      $-\sum_i n_i \sum_y \sum_j \left( O^f_{ij,y} + v \right) ~ ln \left(\hat{O}^f_{ij,y} + v \right)$                                                                                                                                                                                                   Fishery observer data           T4.6  
**Penalties**                                                                                                                                                                                                                                                                                                                                                          
Fishery selectivity          $\sum_i \sum_j^{A_i - 1} \chi \left[ ln \left(\frac{\eta^f_{ij}} {\eta^f_{ij+1}} \right) - ln \left(\frac{\eta^f_{ij+1}} {\eta^f_{ij+2}} \right) \right]^2 , ~ \chi = \left \{  \begin{array}{} 20,~if \eta^f_{ij} \gt \eta^f_{ij+1} \\ 0,~if \eta^f_{ij} \leq \eta^f_{ij+1} \end{array} \right.$                                   T4.7  
**Priors**                                                                                                                                                                                                                                                                                                                                                             
                             $\sum_i \sum_j \left( \tau_{i,y} \right)^2$                                                                                                                                                                                                                                                                                         T4.8  
                             $\sum_i \sum_j \left( N_{0,ij} \right)^2$                                                                                                                                                                                                                                                                                           T4.9  
                             $\sum_i \sum_j \left( \varepsilon_{i,y} \right)^2$                                                                                                                                                                                                                                                                                  T4.10 
                             $v = 0.001$                                                                                                                                                                                                                                                                                                                               

## Errata:
Table 3: Female natural mortality is missing an F under the "Type" column.
Table 4 legend: "function" is misspelled.
Table 4. T.4.6 C is missing the hat in the likelihood equation.
Table 4. T4.7 the upper value of chi, $\eta^f_{ji+1}$ should be $\eta^f_{ij+1}$
Table 4. T4.2 Data source NFMS

## Thoughts
Table 3. Change C* to Y and in table 2 Y to hat{Y}.

# References
