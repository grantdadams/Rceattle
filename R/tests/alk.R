library(Rceattle)
data("BS2017MS")
alk <- BS2017MS$age_trans_matrix[1:21,,3]
numbers <- c(1:21)

lengths <- numbers%*%alk

lengths2 <- rep(0, 25)

for(len in 1:25){
    for(age in 1:21){
    lengths2[len] = lengths2[len] + numbers[age] * alk[age, len]
  }
}
