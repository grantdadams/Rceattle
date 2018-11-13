for (int yr=1;yr<=nyrs;yr++) // year
for (int pred=1;pred<=nspp;pred++) // pred sp
for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age sp
for (int prey=1;prey<=nspp;prey++) // pred sp_age sp
stomKir(yr,pred,pred_age,prey)=UobsAge(pred,prey,pred_age);

for (int yr=1;yr<=nyrs;yr++) // year
for (int pred=1;pred<=nspp;pred++) // pred sp
for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age sp
of_stomKir(yr,pred,pred_age) =(1-sum(stomKir(yr,pred,pred_age)))/other_food(1);


FUNCTION void CALC_SUIT(int msmMode_pass_number)
if(msmMode_pass_number>0){
  // cout<<" *******FIRST PART  =U/(W*N)*******"<<endl;
  dvariable suit_tmp ;
  stom_div_bio2.initialize();
  if(msmMode_pass_number==1)
  {
    // for (int year=1;year<=2;year++) // year
    // {
      //   for (int pred=1;pred<=nspp;pred++) //predator species
      //   {
        //     for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age
        //     {
          //       for (int prey=1;prey<=nspp;prey++) // prey species
          //       {
            //         for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp_age
            //         {
              //          suit_tmp = stom(year,pred,pred_age,prey,prey_age)/(AvgN(prey,stomyrs(year),prey_age));
              //          if (mn_wt_stom(pred,pred_age,prey,prey_age)!=0.)
                //            stom_div_bio(year,pred,pred_age,prey,prey_age) = suit_tmp/mn_wt_stom(pred,pred_age,prey,prey_age);
                //         }  //end prey sp_age
            //       }  // end  prey sp
          //     }  // end pred sp_age
        //   } // end  pred species
      // } // end  year loop
  }else{
    for (int year=1;year<=nyrs;year++) // year
    {
      for (int pred=1;pred<=nspp;pred++) //predator species
      {
        for (int pred_age=1;pred_age<=nages(pred);pred_age++) // pred sp_age
        {
          for (int prey=1;prey<=nspp;prey++) // prey species
          {
            for (int prey_age=1;prey_age<=nages(prey);prey_age++) // prey sp_age
            {
              suit_tmp = stomKir(year,pred,pred_age,prey,prey_age)/(AvgN(prey,year,prey_age));
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
