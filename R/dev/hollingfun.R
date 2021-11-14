# Single-predator
# i = prey
# u = diet composition
# w = weight
# Bo = Biomass of other prey
# H = Holling form

msvpa_fun <- function(wt = 10, u = 0.15, B0 = 500, AvgN = 1000, AvgNPred = 100, ration = 5, h = 1){

  # Step 2 - Get the percent other food in stomachs: (1 - sum(U)) / biomass of other food
  of_stomkir <- (1-sum(u)) / Bo

  # Stomach proportion: U/Biomass and sum(U/Biomass)
  stom_div_bio2 <- u / ( AvgN^h * wt )
  suma_suit <- stom_div_bio2


  # Step 3 - Suitability: (U/Biomass) / (sum(U/Biomass) + (1 - sum(U)) / biomass of other food)
  suit_main <- stom_div_bio2 / (suma_suit + of_stomkir)
  suit_other <- 1 - suit_main

  # Step 4 - Available food: (sum(suitability * biomass) + biomass of other prey * (1 - sum(suitability))) / nyrs
  avail_food = suit_main * AvgN * wt  + Bo*(1 - suit_main)


  # Step 5 -  Calculate M2: N * suitability * ration / available food
  B_eaten = suit_main * AvgN^(h) * wt / (suit_main * AvgN * wt  + Bo*(1 - suit_main)) * AvgNPred * ration
  # B_eaten = AvgNPred * ration * suit_main / avail_food * AvgN^h

  return(B_eaten)
}

curve(msvpa_fun(AvgN = x), from = 0, to = 100)
curve(msvpa_fun(AvgN = x, h = 2), from = 0, to = 100)
curve(msvpa_fun(AvgN = x, h = 2, AvgNPred  = 1 ), from = 0, to = 100, add = TRUE)

AvgN = seq(from = 0, to =3, length.out =  1000)
Holling2 <- msvpa_fun(h = 2, AvgNPred = 5, AvgN = AvgN)
plot(y = Holling2, x = AvgN)
