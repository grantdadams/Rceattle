library(Rceattle)
data(BS2017SS)
ss_run <- Rceattle(data_list = BS2017SS,
                   inits = NULL, # Initial parameters = 0
                   file_name = NULL, # Don't save
                   debug = 0, # Estimate
                   random_rec = FALSE, # No random recruitment
                   msmMode = 0, # Single species mode
                   avgnMode = 0,
                   silent = TRUE)

data(BS2017MS) # ?BS2017MS for more information on the data

ms_run <- Rceattle(
  data_list = BS2017MS,
  inits = ss_run$estimated_params, # Initial parameters from ss run
  file_name = NULL, # Don't save
  debug = 0, # Estimate
  random_rec = FALSE, # No random recruitment
  niter = 1, # Number of iterations around predation/pop dy functions
  msmMode = 1, # Multi-species holsman mode
  avgnMode = 0 # Use average N
)

save(ms_run, file = "MS_run_classic_1iter.Rdata")
