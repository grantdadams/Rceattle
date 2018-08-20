
DATA_INTEGER( nTyrs ); // Number of years of hindcast temperature data; n = [1]
DATA_VECTOR( TempC );  // Vector of bottom temperature; n = [ nTyrs ]
DATA_VECTOR( CA ); // Wt specific intercept of Cmax=CA*W^CB; n = [nspp]
DATA_VECTOR( CB ); // Wt specific slope of Cmax=CA*W^CB; n = [nspp]
DATA_VECTOR( C_model ); // Specifies consumption model if == 1, then use Cmax*fT*P ; n = [nspp]
DATA_VECTOR( fdays ); // Number of foraging days for each predator; n = [nspp]
DATA_ARRAY( Pby_yr );
DATA_VECTOR( Pvalue ); // This scales the pvalue used if C_model ==1 , proportion of Cmax; Pvalue is P in Cmax*fT*Pvalue*PAge; n = [nspp]
DATA_VECTOR( C_model );          // if ==1, the use Cmax*fT*P ; n = [nspp]

array<Type> ConsumAge( nyrs, max_age, nspp ); // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
array<Type> Consum_livingAge( nyrs, max_age, nspp ); // Pre-allocated indiviudal consumption in grams per predator-age; n = [nyrs, nages, nspp]
array<Type> S2Age( nyrs, max_age, nspp ); // Pre-allocate mean stomach weight as a function of sp_age; n = [nyrs, nages, nspp]
array<Type> ration2Age( nyrs, max_age, nspp ); // Annual age-specific ration kg/yr; n = [nyrs, nages, nspp]
matrix<Type> fT(nspp, nTyrs); // Pre-allocate temperature function of consumption - only used if C_model == 1
array<Type> avail_food( nyrs, max_age, nspp ); // Available food
array<Type> stom_div_bio2( nyrs, nspp, max_age, nspp, max_age ); stom_div_bio2.setZero(); // 


  // ------------------------------------------------------------------------- //
  // 5. PREDATION MORTALITY EQUATIONS                                          //
  // ------------------------------------------------------------------------- //
  // 5.1. Calculate ration

  // Set how many years are looped over.
  int nyrsUSE;
  if(fut_pass_number==0){
      nyrsUSE=nyrs;
  }
  if(fut_pass_number==1){  
      nyrsUSE=nyrs_fut;
  }

  // 5.1.1 Calculate historic ration
  for(i=0; i < nspp; i++){
    for(j=0; j < nages(i); j++){
      for(y=0; y < nyrs; y++){ // Caclulate ration for each species

      if(fut_pass_number==0){
        ConsumAge(y, j, i) = Type(24) * Type(0.0134) * exp( Type(0.0115) * TempC( y )) * Type(91.25) * S2Age(y, j, i) * wt(y, j, i); // Calculate consumption for predator-at-age; units = kg/predator
        Consum_livingAge(y, j, i) = ConsumAge(y, j, i); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
        

        if(C_model( i )==1){
          ConsumAge(y, j, i) = CA(i) * pow(wt(y, j, i) * 1000, CB( i )) * fT(i, y) * fdays( i ) * wt(y, j, i) * 1000;//g/pred.yr      
          //ConsumAge(predd,yrr)=elem_prod(ConsumAge(predd,yrr),Pvalue(predd)*PAge(predd));
          ConsumAge(y, j, i) = ConsumAge(y, j, i) * Pvalue(i) * Pby_yr(y, j, i) // elem_prod(ConsumAge(predd,yrr),Pvalue(predd)*Pby_yr(predd,yrr));
        }

        ration2Age(y, j, i) = ConsumAge(y, j, i) / 1000; //annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
      }


     // if(fut_pass_number==1){  
        
      //  ConsumAge_fut(predd,yrr)=24*0.0134*mfexp(0.0115*TempC_futUSE(itemp,yrr))*91.25*elem_prod(S2Age(predd,nyrs),wt_fut(predd,yrr));//kg/pred.yr
      //  Consum_livingAge_fut(predd,yrr)=ConsumAge_fut(predd,yrr); // vector of specific consumption rates for each size of predator in terms of g/g/yr of all prey consumed per predator per year
       // if(C_model(predd)==1){
       //   ConsumAge_fut(predd,yrr)=elem_prod(((CA(predd))*pow(wt_fut(predd,yrr)*1000,CB(predd))*fT_fut(itemp,predd,yrr)*fdays(predd)),wt_fut(predd,yrr)*1000);//g/pred.yr      
          //ConsumAge_fut(predd,yrr)=elem_prod(ConsumAge_fut(predd,yrr),Pvalue(predd)*PAge(predd));
        //  ConsumAge_fut(predd,yrr)=elem_prod(ConsumAge_fut(predd,yrr),Pvalue(predd)*Pby_yr(predd,nyrs));
       // }
        //ration2Age_fut(predd,yrr)=ConsumAge_fut(predd,yrr)/1000;//annual ration kg/yr //aLW(predd)*pow(lengths(predd,j),bLW(predd));//mnwt_bin(predd,j);
        
      //}  
    }//end predd
  }
}



  // 5.2. Calculate M2: 
  // NOTE: prey "i" of age "j"; predator "p" at age "a"
  // NOTE: Does not include overlap
      Type Mtmp = 0; // Temporary M2 
      Type Eatentmp = 0; // Temporary diet

  for (i = 0; i < nspp; i++){  // prey spp loop
    for (j = 0; j < nages(i); j++){  // prey sp_age loop
      for(y = 0; y < nyrs; y++){
      Mtmp = 0;
      Eatentmp = 0;
      for (p = 0; p < nspp; p++){   // pred species loop  
        for (a = 0; a < nages(p); a ++){  // Pred sp_age loop

          //if(fut_pass_number==0){
              Mtmp += (AvgN(y, a, p) * ration2Age(y, a, p ) * suit_main(p , a, i, j)) / avail_food(y, a, p);
              Eatentmp += (AvgN(y, a, p) * ration2Age(y, a, p) * suit_main(p , a, i, j));
            
          //}
          //if(fut_pass_number==1){
            
              //Mtmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age))/avail_food(pred,nyrs,pred_age);
              //Eatentmp += (AvgN_fut(pred,i,pred_age)*overlap_fut(itemp,pred,prey,i) * ration2Age_fut(pred,i,pred_age)*suit_main(pred,pred_age,prey,prey_age));
            //}
          
          
        }  // end pred sp_age loop
      }  // end pred spp loop
      //if(fut_pass_number == 0){
        M2(y, j, i) = Mtmp; 
        B_eaten(y, j, i) = Eatentmp;
      //} 
      //if(fut_pass_number==1){
      //M2_fut(prey,i,prey_age) = Mtmp; //M2_plk(prey,i,prey_age);  
      //B_eaten_fut(prey,i,prey_age)=Eatentmp;
      }
    }   // end prey sp_age loop  
  }   // end prey spp loop




  // 5.3. Calculate Suitability: 
    Type suit_tmp ;
    
      for (y = 0; y < nyrs; y++){ // year
      for (p = 0; p < nspp; p++){ //predator species
        for (a = 0; a < nages(p); a++){ // pred sp_age 
        for (i = 0; i < nspp; i++){ // prey species 
          for (j = 0; j < nages(i); j++){ // prey sp_age 
           
           suit_tmp = stomKir(year,pred,pred_age,prey,prey_age) / AvgN(y, a, p); 
           if (wt(prey,year,prey_age)!=0.)
             stom_div_bio2(year,pred,pred_age,prey,prey_age) = suit_tmp/wt(prey,year,prey_age);    
  
  //         if(AvgN(prey,year,prey_age)!=0)
  //           stom_div_bio2(year,pred,pred_age,prey,prey_age) = (UobsAge(pred,prey,pred_age,prey_age)*ration2Age(pred,year,pred_age))/(AvgN(prey,year,prey_age)*wt(prey,year,prey_age));    
          }  //end prey sp_age 
        }  // end  prey sp
        }  // end pred sp_age    
      } // end  pred species 
      } // end  year loop  
   }
    //end if kirs_pass_num==0
    // cout<<"************suitabilities calculation*******"<<endl;
    suit_main.initialize(); // sets to zero
    dvariable suma_suit;
    for (int pred=1;pred<=nspp;pred++) //predator species
    {
    for (int pred_age=1;pred_age<=nages(pred);pred_age++) // predator sp 
    {
      for (int prey=1;prey<=nspp;prey++) // prey specie
      {
      for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp
      { 
        switch (msmMode_pass_number)
        {
          case 0:
          break;
          case 1:
            // for (int year=1;year<=nstom;year++) // year
            // {
            //   suma_suit=sum(stom_div_bio(year,pred,pred_age)); // sum of aitemp prey and prey ages in the stom of the pred sp_age j
            //   suit_std(year,pred,pred_age,prey,prey_age) = stom_div_bio(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stom(year,pred,pred_age));  //(sum(stom_div_bio(year,pred,pred_age)+of_stom(year,pred,pred_age)));            
            //   suit_main(pred,pred_age,prey,prey_age)    += suit_std(year,pred,pred_age,prey,prey_age);                                  
            // }// end year
            // suit_main(pred,pred_age,prey,prey_age) /= nstom; 
            cout<<"!!!! THIS MODE NO LONGER USED (MSM MODE = 1) "<<endl;
          break;
          case 2:
            suit_main(pred,pred_age,prey,prey_age)=0;
            for (int year=1;year<=nyrs;year++) // year
            {
              suma_suit=sum(stom_div_bio2(year,pred,pred_age)); // sum of all prey and prey ages in the stom of the pred sp_age j
              suit_main(pred,pred_age,prey,prey_age)    += stom_div_bio2(year,pred,pred_age,prey,prey_age)/(suma_suit+of_stomKir(year,pred,pred_age));                               
            }// end year
            suit_main(pred,pred_age,prey,prey_age) /= nyrs; 
            
          break;
         // REMOVED A3
        }  //end switch              
      }// prey sp_age 
      } // prey
      suit_other(pred,pred_age) = 1. - sum(suit_main(pred,pred_age)); //estimate the other food suitability        
    }// pred sp_age
    }  // pred
  }//if msmMode >0

  FUNCTION void CALC_AVAIL_FOOD(int msmMode_pass_number)
  // cout<<"CALC AVAILABLE FOOD ******************************************************************************"<<endl;
  for (int pred=1; pred<= nspp; pred++)  // predator loop
  {
    for (int pred_age=1;pred_age<=nages(pred);pred_age++)  // predator sp_age loop
    {
      dvariable tmpsuit;
      dvariable tmp_othersuit;
      tmp_othersuit=0.;
      tmpsuit=0.;
      tmpsuit=1-suit_other(pred,pred_age);
      avail_food(pred,i,pred_age)=0.;    
      if(msmMode_pass_number==1){
        for (int prey=1;prey<=nspp;prey++)
        {
          for (int prey_age =1;prey_age<=nages(prey);prey_age++)
            avail_food(pred,i,pred_age) += suit_main(pred,pred_age,prey,prey_age)*overlap(i,pred,prey) *AvgN(prey,i,prey_age) * mn_wt_stom(pred,pred_age,prey,prey_age);   
          tmp_othersuit+=sum(suit_main(pred,pred_age,prey))*overlap(i,pred,prey);
        }
        avail_food(pred,i,pred_age) += other_food(pred)*(1.-(tmp_othersuit));
          //avail_food(pred,i,pred_age) += other_food(pred)* suit_other(pred,pred_age); 
        }else{
          for (int prey=1;prey<=nspp;prey++)
          {
            for (int prey_age =1;prey_age<=nages(prey);prey_age++)
              avail_food(pred,i,pred_age) += suit_main(pred,pred_age,prey,prey_age)*overlap(i,pred,prey) *AvgN(prey,i,prey_age) *wt(prey,i,prey_age); 
            tmp_othersuit+=sum(suit_main(pred,pred_age,prey))*overlap(i,pred,prey);
          }
          avail_food(pred,i,pred_age) += other_food(pred)*(1.-(tmp_othersuit));
          
          //avail_food(pred,i,pred_age) += other_food(pred)* (suit_other(pred,pred_age));     
        }  
      }              // end pred sp_age loop
  }                // end pred loop