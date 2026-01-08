# Dimensions
nlengths = 60
lengths = 1:nlengths
nages = 10
pred_CAAL1 <- matrix(0, nlengths, nages)
pred_CAAL2 <- matrix(0, nlengths, nages)
pred_CAAL3 <- matrix(0, nlengths, nages)

# Length-based selectivity
sel_at_length  <- 1/(1+exp(-0.3 * (lengths - 35)))

# Growth matrix (length-to-age)
growth_matrix <- matrix(rnorm(nlengths*nages), nlengths, nages)
growth_matrix <- apply(growth_matrix, 2, function(x) x/sum(x))

# Population characteristics
NAA <- 10:1
MAA <- rep(0.2, nages)
Frate <- 0.3

# Calculate total mortality at age
sel_at_age <- t(growth_matrix) %*% sel_at_length # Convert selectivity from length-based to age-based
FAA = sel_at_age * Frate
ZAA = MAA + sel_at_age * Frate

# Calculate catch-at-age/length
for(l in 1:nlengths) {
  for(a in 1:nages) {
    # approach 1
    pred_CAAL1[l,a] = sel_at_length[l] * growth_matrix[l,a] * Frate *  NAA[a] / ZAA[a] * (1-exp(-ZAA[a]))

    # approach 2
    ztemp <- MAA[a] + sel_at_length[l] * growth_matrix[l,a] * Frate
    pred_CAAL2[l,a] =  sel_at_length[l] * Frate * growth_matrix[l,a] * NAA[a]  /  ztemp * (1-exp(-ztemp))

    # approach 3
    pred_CAAL3[l,a] = FAA[a] * growth_matrix[l,a] *  NAA[a] / ZAA[a] * (1-exp(-ZAA[a]))
  }
}

sum(NAA / ZAA * (1-exp(-ZAA)) * FAA)
sum(pred_CAAL)
sum(pred_CAAL2)
sum(pred_CAAL3)
