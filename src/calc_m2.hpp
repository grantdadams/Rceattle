
DATA_INTEGER( nTyrs ); // Number of years of hindcast temperature data; n = [1]
DATA_VECTOR( TempC );  // Vector of bottom temperature; n = [ nTyrs ]
DATA_VECTOR( CA ); // Wt specific intercept of Cmax=CA*W^CB; n = [nspp]
DATA_VECTOR( CB ); // Wt specific slope of Cmax=CA*W^CB; n = [nspp]
DATA_VECTOR( C_model ); // Specifies consumption model if == 1, then use Cmax*fT*P ; n = [nspp]
DATA_VECTOR( fdays ); // Number of foraging days for each predator; n = [nspp]

array<Type> ConsumAge( nyrs, max_age, nspp ); // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
array<Type> Consum_livingAge( nyrs, max_age, nspp ); // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
array<Type>  S2Age( nyrs, max_age, nspp ); // Pre-allocate mean stomach weight as a function of sp_age; n = [nyrs, nages, nspp]

matrix<Type>  fT(nspp, nTyrs); // pre-allocate temperature function of consumption - only used if C_model == 1



FUNCTION void CALC_RATION(int fut_pass_number)
  int nyrsUSE;
  if(fut_pass_number==0)
    nyrsUSE=nyrs;
  if(fut_pass_number==1)  
    nyrsUSE=nyrs_fut;



  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){ // Caclulate ration for each species

      if(fut_pass_number==0){
        ConsumAge(y, j, i) = Type(24) * Type(0.0134) * exp(0.0115 * TempC( y )) * Type(91.25) * S2Age(y, j, i) * wt(y, j, i); // Calculate consumption for predator-at-age; units = kg/predator
        Consum_livingAge(y, j, i) = ConsumAge(y, j, i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        
        if(C_model(i)==1){
          ConsumAge(y, j, i)= CA(i) * pow(wt(y, j, i) * 1000, CB( i )) * fT(i, y) * fdays( i ) * wt(y, j, i) * 1000;//g/pred.yr      
          //ConsumAge(predd,yrr)=elem_prod(ConsumAge(predd,yrr),Pvalue(predd)*PAge(predd));
          ConsumAge(predd,yrr)=elem_prod(ConsumAge(predd,yrr),Pvalue(predd)*Pby_yr(predd,yrr));
        }
        ration2Age(predd,yrr)=ConsumAge(predd,yrr)/1000;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }
      if(fut_pass_number==1){  
        
        ConsumAge_fut(predd,yrr)=24*0.0134*mfexp(0.0115*TempC_futUSE(itemp,yrr))*91.25*elem_prod(S2Age(predd,nyrs),wt_fut(predd,yrr));//kg/pred.yr
        Consum_livingAge_fut(predd,yrr)=ConsumAge_fut(predd,yrr); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        if(C_model(predd)==1){
          ConsumAge_fut(predd,yrr)=elem_prod(((CA(predd))*pow(wt_fut(predd,yrr)*1000,CB(predd))*fT_fut(itemp,predd,yrr)*fdays(predd)),wt_fut(predd,yrr)*1000);//g/pred.yr      
          //ConsumAge_fut(predd,yrr)=elem_prod(ConsumAge_fut(predd,yrr),Pvalue(predd)*PAge(predd));
          ConsumAge_fut(predd,yrr)=elem_prod(ConsumAge_fut(predd,yrr),Pvalue(predd)*Pby_yr(predd,nyrs));
        }
        ration2Age_fut(predd,yrr)=ConsumAge_fut(predd,yrr)/1000;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
        
      }  
    }//end predd
  }





FUNCTION void CALC_M2(int fut_pass_number)
  CALC_RATION(fut_pass_number);
  for (int prey=1;prey<=nspp;prey++)  // prey spp loop
  {
    for (int prey_age =1;prey_age<=nages(prey);prey_age++)  // prey sp_age loop
    {
      dvariable Mtmp=0.; 
      dvariable Eatentmp=0.;  
      for (int pred=1;pred<=nspp;pred++)   // pred species loop
      {  
        for (int pred_age=1;pred_age<=nages(pred);pred_age++)  // Pred sp_age loop
        {
          if(fut_pass_number==0)
          {
            if(msmMode==1)
            {
              //Mtmp += (AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,i,pred_age);
            }else
            {
              Mtmp += (AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,i,pred_age);
              Eatentmp += (AvgN(pred,i,pred_age)*overlap(i,pred,prey) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age));
            }
          }
          if(fut_pass_number==1)
          {
            if(msmMode==1)
            {
              //Mtmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration(pred,nyrs,pred_age)*suit_main(pred,pred_age,prey,prey_age))/av_food_fut(pred,i,pred_age);
            }else{
              Mtmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,nyrs,pred_age);
              Eatentmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age));
              // Mtmp += AvgN_fut(pred,i,pred_age)*overlap_fut(i,pred,prey) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/(av_food_fut(pred,i,pred_age));
            }
          }
          /*switch (Kirs)
          {
            case 0:
              if(fut_pass_number==0)
                Mtmp += AvgN(pred,i,pred_age) * ration(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/avail_food(pred,i,pred_age);
              if(fut_pass_number==1)
                Mtmp += AvgN_fut(pred,i,pred_age) * ration(pred,nyrs,pred_age)*suit_main(pred,pred_age,prey,prey_age)/av_food_fut(pred,i,pred_age);  
            break;
            case 1:
            
              if(fut_pass_number==0)
                Mtmp += AvgN(pred,i,pred_age) * ration2Age(pred,i,pred_age)*UobsAge(pred,prey,pred_age,prey_age)/avail_food(pred,i,pred_age);
              if(fut_pass_number==1)
                Mtmp += AvgN_fut(pred,i,pred_age) * ration2Age_fut(pred,i,pred_age)*UobsAge(pred,prey,pred_age,prey_age)/av_food_fut(pred,i,pred_age);
          
            break;
            default:
              if(fut_pass_number==0)
                Mtmp += AvgN(pred,i,pred_age) * ration2Age(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/avail_food(pred,i,pred_age);
              if(fut_pass_number==1)
                Mtmp += AvgN_fut(pred,i,pred_age) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age)/av_food_fut(pred,i,pred_age);
            break;
          }   */     
        }  // end pred sp_age loop
      }  // end pred spp loop
      if(fut_pass_number==0){
        M2(prey,i,prey_age) = Mtmp; //M2_plk(prey,i,prey_age);
        B_eaten(prey,i,prey_age)=Eatentmp;
      } 
      if(fut_pass_number==1){
        M2_fut(prey,i,prey_age) = Mtmp; //M2_plk(prey,i,prey_age);  
        B_eaten_fut(prey,i,prey_age)=Eatentmp;
      }
    }   // end prey sp_age loop  
  }   // end prey spp loop