#include <TMB.hpp>

// Function written by Anders Nielsen to reorder the categories which is needed for OSA residuals of vector-valued observations
// that do not have likelihoods defined in TMB.
template<class Type>
vector<int> order(vector<Type> k){
  int n=k.size();
  vector<int> ord(n);
  ord.setZero();
  int at=-1;
  for(int i=0; i<n;++i){
    if(k(i)>0.5){ord(++at) = i;}
  }
  at=n;
  for(int i=n-1; i>=0;--i){
    if(k(i)<0.5){ord(--at) = i;}
  }
  return ord;
}

template <class Type>
Type dmultinom(vector<Type> x, vector<Type> p, vector<int> ages, data_indicator<vector<Type>, Type> keep, int give_log, int do_osa)
{
  Type logres = 0;
  vector<Type> p_x(ages.size());
  for(int i = 0; i < ages.size(); i++) p_x(i) = p(ages(i)-1);
  p_x /= sum(p_x);

  if(do_osa){
    vector<Type> k=keep;
    //vector<Type> l=keep.cdf_lower;
    //vector<Type> h=keep.cdf_upper;
    vector<int> o=order(k);
    x=x(o); p_x=p_x(o); k=k(o); //l=l(o); h=h(o);
    Type nUnused=asDouble(sum(x));
    Type pUsed=0;
    //Type cdf;
    for(int i=0; i<x.size(); ++i){
      if(i!=(x.size()-1)){
        vector<Type> x2(2), p2(2);
        //Type p_i = squeeze(p_x(i));
        //Type one_minus_pUsed_i = squeeze(1.0-pUsed);
        x2(0) = x(i);
        x2(1) = nUnused-x(i);
        p2(0) = squeeze(p_x(i));
        p2(0) /= squeeze((Type(1)-pUsed));
        p2(0) = squeeze(p2(0));
        // p2(0) = squeeze(p(i)/(Type(1)-pUsed));//(Type(1)-pUsed_i); //for log of any p = 0
        p2(1) = 1. - p2(0);
        logres += k(i) * dmultinom(x2,p2,1); //binomial the hard way.
        //logres += k(i)*dbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true);
        //cdf = pbinom(x(i),nUnused,p(i)/(Type(1)-pUsed));
        nUnused -= x(i);
        pUsed += p_x(i);
      }else{ // last index
        logres += k(i)*Type(0);
        //cdf = Type(1);
      }
      //cdf = squeeze(cdf);
      //logres += l[i] * log( cdf );       // NaN protected
      //logres += h[i] * log( 1.0 - cdf ); // NaN protected
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
