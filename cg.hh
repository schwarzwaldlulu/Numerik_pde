#ifndef CG_HH
#define CG_HH

#include <cmath>

using namespace std;

// using conjugate gradient method to calculate the system of linear equations Ax = b, where
// A is symmetric and positive definite, usually sparse matrix
// input vector x can be an approximate initial solution or 0
// eps: tolerance, i.e. eps*(b*b) should be larger than the diffenrence between Ax and b, where x is our approximate solution

// our implementation updates the approximate solution x in-place and returns the times of iterations

template< class Matrix, class Vector >
int cg ( const Matrix &A, const Vector &b, Vector &x, double eps, bool detailed = false )
{
  int count = 0;

  double error = eps * sqrt( b * b );

  // first calculate residum r = b - A * x
  Vector r( x );

  A( x, r ); // r = A * x
  r -= b;  // r = r - b = A * x - b
  r *= -1.;  // r = b - A * x

  Vector p( r ), z( r );

  if( detailed )
    cout<<"cghs intial \t"<< r * r << endl;
  
  double residum_norm_squared = r * r;

  while( residum_norm_squared > error )
  {
    A( p, z ); // z = A * p
    
    double alpha = residum_norm_squared / ( p * z ); // alpha = (r*r)/(p*A*p)
    
    x.axpy( alpha, p ); // update x
    
    r.axpy( -alpha, z ); // update residum r
    
    double residum_norm_squared_new = r * r;

    p *= residum_norm_squared_new / residum_norm_squared;
    
    p += r; // update p
    
    residum_norm_squared = residum_norm_squared_new;

    if( detailed )
      cout<<"cghs "<< count << "\t" << residum_norm_squared << endl;
    
    ++count;
  }

  return count;
}


#endif // #ifndef CG_HH
