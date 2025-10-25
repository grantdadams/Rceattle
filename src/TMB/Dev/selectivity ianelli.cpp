//////////////////////////////////////////////////////////////////////////////////
  //--If time changes turned on then do 2nd differencing
//  for curvature penalty on subsequent years, otherwise only first year matters
sel_like_dev.initialize();  // Initialize penalty vector to zeros

if (active(sel_coffs_fsh))  // Check if fishery selectivity coefficients are being estimated
{
  if (active(sel_devs_fsh))  // Check if time-varying selectivity deviations are being estimated
  {
    // Penalty term 1: Overall penalty on the magnitude of selectivity deviations
    // This improves estimability by shrinking deviations toward zero when not informed by data
    // Scaled by ctrl_flag(10) parameter and divided by number of fishery groups
    sel_like_dev(1) += ctrl_flag(10)/group_num_fsh * norm2(sel_devs_fsh);  // norm2() calculates sum of squared values

    // Penalty term 2: "Smoothness" penalty on the base selectivity curve (in the start year)
    // Uses second differences (approximating second derivatives) to penalize curves with high curvature
    // This encourages smoother selectivity patterns
    // Scaled by ctrl_flag(11) and divided by number of change years to maintain consistent penalty weight
    sel_like_dev(1) += ctrl_flag(11)/nch_fsh * norm2(first_difference(first_difference(log_sel_fsh(styr))));

    // Loop through each selectivity change year to add penalties
    for (i=1; i<=nch_fsh; i++)
    {
      // Penalty term 3: "Smoothness" penalty for each change year's selectivity curve
        // Again uses second differences to penalize high curvature
        sel_like_dev(1) += ctrl_flag(11)/nch_fsh * norm2(first_difference(first_difference(log_sel_fsh(yrs_ch_fsh(i)))));

        // Penalty term 4: Year-to-year stability penalty
        // Penalizes large changes between a change year and the preceding year
        // Only applies if we're not in the start year (which has no preceding year)
      // This is a "random walk" penalty weighted by a variance parameter (sel_ch_sig_fsh)
      if (yrs_ch_fsh(i) != styr)
        sel_like_dev(1) += norm2(log_sel_fsh(yrs_ch_fsh(i)-1) - log_sel_fsh(yrs_ch_fsh(i))) /
          (2 * sel_ch_sig_fsh(i) * sel_ch_sig_fsh(i));  // Division by 2*sigma² is standard for normal likelihood
    }
  }
  else  // If no time-varying deviations, only apply smoothness penalty to the base selectivity
  {
    sel_like_dev(1) += ctrl_flag(11) * norm2(first_difference(first_difference(log_sel_fsh(styr))));
  }
}





// Function to compute time-varying fishery selectivity using a coefficient approach
FUNCTION dvar_matrix compute_fsh_selectivity(const int nsel,          // Number of age classes with unique selectivity coefficients
                                             const int stsel,          // Start year for selectivity
                                             dvariable& avgsel,        // Will store average selectivity (returned by reference)
                                             const dvar_vector& coffs, // Vector of selectivity coefficients
                                             const dvar_matrix& sel_devs, // Matrix of selectivity deviations
                                             int gn)                   // Gap between years (appears unused)
{
  // Coefficient selectivity, with deviations...
  RETURN_ARRAYS_INCREMENT();  // ADMB memory management function

  // Create matrix to store log-selectivity values by year and age
  dvar_matrix log_sel(styr,endyr_r,1,nages);  // Dimensions: years (styr to endyr_r) × ages (1 to nages)
  log_sel.initialize();  // Initialize matrix to zeros

  // Calculate mean of exponentiated coefficients and store the log
  avgsel = log(mean(mfexp(coffs)));

  // Set selectivity for start year:
  // For ages 1 to nsel: use the coefficients directly
  log_sel(stsel)(1,nsel) = coffs;
  // For ages nsel+1 to nages: use the last coefficient (flat selectivity for older ages)
  log_sel(stsel)(nsel+1,nages) = coffs(nsel);

  int ii;
  // Normalize the selectivity in start year to have mean = 1 (log mean = 0)
  log_sel(stsel) -= log(mean(exp(log_sel(stsel))));

  ii = 1;  // Counter for change years

  // Loop through all years to apply time-varying changes to selectivity
  for (i=stsel; i<endyr_r; i++)
  {
    // Check if this is a year where selectivity changes
    if (ii<=nch_fsh)  // If we haven't exceeded the number of change years
    {
      if (i==yrs_ch_fsh(ii))  // If this is a change year
      {
        // Apply selectivity deviation to the first nsel ages
        log_sel(i+1)(1,nsel) = log_sel(i)(1,nsel) + sel_devs(ii);
        // Copy the selectivity coefficient for the last specified age to all older ages
        log_sel(i+1)(nsel+1,nages) = log_sel(i+1,nsel);
        ii++;  // Move to next change year
      }
      else
        // If not a change year, copy from previous year
        log_sel(i+1) = log_sel(i);
    }
    else
      // If we've used all change years, just copy from previous year
      log_sel(i+1) = log_sel(i);

    // Normalize the selectivity in each year to have mean = 1 (log mean = 0)
    log_sel(i+1) -= log(mean(exp(log_sel(i+1))));
  }

  RETURN_ARRAYS_DECREMENT();  // ADMB memory management function
  return(log_sel);  // Return the log-selectivity matrix
}
