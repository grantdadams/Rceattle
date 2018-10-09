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
