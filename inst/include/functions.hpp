
// Function for getting max of an IVECTOR
template <class Type>
Type imax(const vector<Type> &x)
{
  int res = x[0];
  for(int i=0; i < x.size(); i++){
    int res2 = x[i];
    if(res < res2){
      res = res2;
    }
  }
  return res;
}

// Function for getting min of an IVECTOR
template <class Type>
Type imin(const vector<Type> &x)
{
  int res = x[0];
  for(int i=0; i < x.size(); i++){
    int res2 = x[i];
    if(res > res2){
      res = res2;
    }
  }
  return res;
}

// Function for mean of a vector
template <class Type>
Type mean(const vector<Type> &x)
{
  Type mean = x.sum()/x.size();
  return mean;
}

// Function to calculate the first difference
template <class Type>
vector<Type> first_difference(const vector<Type> &x)
{
  int length = x.size() - 1;

  if (length <= 0){
    std::cerr << "Error -- vector size too small in first_difference(const dvector&)" << std::endl;
    return(0);
  }

  vector<Type> tmp(length);

  for (int i = 0; i < length; i++){
    tmp(i) = x(i+1) - x(i);
  }

  return tmp;
}

// Function for elementwise division
template <class Type>
matrix<Type> elem_div(matrix<Type> m1, matrix<Type> m2){

  int m1r = m1.rows();
  int m2r = m2.rows();

  int m1c = m1.cols();
  int m2c = m1.cols();

  if (m1r != m2r){
    std::cerr << "Error -- number of rows in matrices does not match" << std::endl;
    return(0);
  }

  if (m1c != m2c){
    std::cerr << "Error -- number of columns in matrices does not match" << std::endl;
    return(0);
  }

  matrix<Type> m3(m1r, m1c);

  // Elementwise division
  for(int r = 0; r < m1r; r++){
    for(int c = 0; c < m1c; c++){
      m3(r, c) = m1(r, c) / m2(r, c);
    }
  }
  return m3;
}

// Function for elementwise matrix exponential functions
template <class Type>
matrix<Type> elem_pow(matrix<Type> m1, Type exponent){

  int m1r = m1.rows();
  int m1c = m1.cols();

  matrix<Type> m2(m1r, m1c);

  // Elementwise division
  for(int r = 0; r < m1r; r++){
    for(int c = 0; c < m1c; c++){
      m2(r, c) = pow( m1(r, c), exponent);
    }
  }
  return m2;
}

// Function for elementwise addition in a matrix
template <class Type>
matrix<Type> elem_add(matrix<Type> m1, Type value){

  int m1r = m1.rows();
  int m1c = m1.cols();

  matrix<Type> m2(m1r, m1c);

  // Elementwise division
  for(int r = 0; r < m1r; r++){
    for(int c = 0; c < m1c; c++){
      m2(r, c) =  m1(r, c) + value;
    }
  }
  return m2;
}

