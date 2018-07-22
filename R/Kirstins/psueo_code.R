psuedo code for CEATTLEv2.0

- include some method for "peel back"
- track year indices of input data seperately from simulation years


0 estimate weight at age using Wt_mode (1= mean obs, 2 = vonB, 3= vonB+cov)
1 run single species mode
2 run esimates of M based on biomass
3 iterate x 3
4 estimate M parameters (maybe using seperate code snippet
5 estimate recuitment
6 project with set estimated values for parameters
	 a. MSE mode 	--> run forward each year
	 				--> simulate catch and survey
	 				--> update input data, and project again
	 b. not MSE 	--> project forward nyrs_fut under F rate/ control rule


TMB Modules:
(0) vonB
(1) Catage
(2) MbyY
(3) CE_Rec

Switches:
1. multispp mode
2. fsh sel
3. srvy sel
4. rec_mode
5. estimate M or calc





