//function written by Anders Nielsen to reorder the categories which is needed for OSA residuals of vector-valued observations
// that do not have likelihoods defined in TMB.
template<class Type>
vector<int> order(vector<Type> k){
  int n=k.size();
  vector<int> o(n);
  o.setZero();
  int at=-1;
  for(int i=0; i<n;++i){
    if(k(i)>0.5){o(++at) = i;}
  }
  at=n;
  for(int i=n-1; i>=0;--i){
    if(k(i)<0.5){o(--at) = i;}
  }
  return o;
}



template <class Type>
Type dmultinom(vector<Type> x, vector<Type> p, vector<int> ages, data_indicator<vector<Type>, Type> keep, int give_log, int do_osa)
{
  Type logres = 0;
  vector<Type> p_x(ages.size());
  for(int i = 0; i < ages.size(); i++) p_x(i) = p(ages(i)-1);
  p_x /= sum(p_x);

  if(do_osa){

    // - Reorder
    vector<Type> k=keep;
    vector<int> o=order(k);
    x=x(o);
    p_x=p_x(o);
    k=k(o); //l=l(o); h=h(o);


    Type nUnused=asDouble(sum(x));
    Type pUsed=0;
    //Type cdf;
    for(int i=0; i<x.size(); ++i){
      if(i!=(x.size()-1)){
        vector<Type> x2(2), p2(2);

        x2(0) = x(i);
        x2(1) = nUnused-x(i); // sum(x) - x

        // - [0,1] to (0,1)
        p2(0) = squeeze(p_x(i));
        p2(0) /= squeeze((Type(1)-pUsed));
        p2(0) = squeeze(p2(0));

        p2(1) = 1. - p2(0);
        logres += k(i) * dmultinom(x2,p2,1); //binomial the hard way.

        nUnused -= x(i);
        pUsed += p_x(i);

      }else{ // last index
        logres += k(i)*Type(0);
      }
    }
  } else {
    logres = dmultinom(x,p_x,1);
  }
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}
