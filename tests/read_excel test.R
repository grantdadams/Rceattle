data("BS2017MS")

# Write data to excel
Rceattle::write_excel(data_list = BS2017MS, file = "BS2017MS.xlsx")

# Change the data how you want in excel
# Read the data back in
mydata2 <- Rceattle::read_excel( file = "BS2017MS.xlsx")


data_names <- names(mydata2)[names(mydata2) %in% names(BS2017MS)]

equal_check <- c()
for(i in 1:length(data_names)){
  equal_check[i] <- sum(mydata2[[data_names[i]]] != BS2017MS[[data_names[i]]], na.rm = TRUE)
}
