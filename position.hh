#ifndef POSITION_HH
#define POSITION_HH

#include <algorithm>
#include <array>
#include <iostream>


using namespace std;


// using array<double, 2> to represent points in dimension 2
// define basic operations and overloaded operators
class Position
{
  static const int size = dimension;

public:
  typedef array< double, size >::iterator iterator;
  typedef array< double, size >::const_iterator const_iterator;

  Position ( double s = 0 )
  {
    data_.fill( s );
  }

  Position ( const Position& other )
  {
    copy( other.data_.begin(), other.data_.end(), data_.begin() );
  }
  
  Position& operator= ( const double s )
  {
    data_.fill( s );
    return *this;
  }

  Position& operator= ( const Position& other )
  {
    copy( other.data_.begin(), other.data_.end(), data_.begin() );
    return *this;
  }


  double& operator[] ( const int i ) { return data_[ i ]; }
  const double& operator[] ( const int i ) const { return data_[ i ]; }


  iterator begin () { return data_.begin(); }
  const_iterator begin () const { return data_.begin(); }

  iterator end () { return data_.end(); }
  const_iterator end () const { return data_.end(); }

  // overloading compound assignment operator 
  Position& operator+= ( const Position& other )
  {
    const_iterator oIt = other.begin();
    for( iterator it = begin(); it != end(); ++it, ++oIt )
      *it += *oIt;
    return *this;
  }
  
  Position& operator-= ( const Position& other )
  {
    const_iterator oIt = other.begin();
    for( iterator it = begin(); it != end(); ++it, ++oIt )
      *it -= *oIt;
    return *this;
  }

  Position& operator*= ( const double s )
  {
    for( double &v : data_ )
      v *= s;
    return *this;
  }

  Position& operator/= ( const double s )
  {
    for( double &v : data_ )
      v /= s;
    return *this;
  }

  // overloading arithmetic operators
  Position operator+ ( const Position& other ) const
  {
    Position z( *this );
    return z += other;
  }
  
  Position operator- ( const Position& other ) const
  {
    Position z( *this );
    return z -= other;
  }

  Position operator* ( double s ) const
  {
    Position z( *this );
    return z *= s;
  }
  
  friend Position operator* ( double s, const Position& other )
  {
    Position p( other );
    return p *= s;
  }

  Position operator/ ( double s ) const
  {
    Position z( *this );
    return z /= s;
  }

  // scalar product
  double operator* ( const Position &other ) const
  {
    double value = 0.;
    const_iterator oIt = other.begin();
    for( const_iterator it = begin(); it != end(); ++it, ++oIt )
      value += ( *it ) * ( *oIt );
    return value;
  }

  // add multiple of another position to the receive object
  void axpy ( double s, const Position &other )
  {
    const_iterator oIt = other.begin();
    for( iterator it = begin(); it != end(); ++it, ++oIt )
      *it += s * (*oIt);
  }

protected:
  array< double, size > data_;
};


// read in position
istream& operator>> ( istream& in, Position &v )
{
  Position w;
  for( double &e : w )
    in >> e;
  if( in )
    v = w;
  return in;
}


// print position
ostream& operator<< ( ostream& out, const Position &v )
{
  for( const double &e : v )
    out << e <<" ";
  return out;
}
#endif // #ifndef POSITION_HH
