//  mostly sample code from the tutor

#ifndef VECTOR_HH
#define VECTOR_HH

#include <vector>

class Vector
{
  typedef std::vector< double > Data;

public:

  typedef Data::iterator iterator;
  typedef Data::const_iterator const_iterator;

  Vector ( const int size, double val = 0 )
    : data_( size, val )
  {}

  Vector ( const Vector &other ) = default;

  operator Data & () { return data_; }
  operator const Data & () const { return data_; }

  iterator begin () { return data_.begin(); }
  const_iterator begin () const { return data_.begin(); }

  iterator end () { return data_.end(); }
  const_iterator end () const { return data_.end(); }

  Vector &operator= ( const double s )
  {
    for( double &v : data_ )
      v = s;
    return *this;
  }

  Vector &operator= ( const Vector &other )
  {
    std::copy( other.data_.begin(), other.data_.end(), data_.begin() );
    return *this;
  }

  Vector &operator+= ( const Vector &other )
  {
    const_iterator oIt = other.begin();
    for( iterator it = begin(); it != end(); ++it, ++oIt )
      *it += *oIt;
    return *this;
  }

  Vector &operator-= ( const Vector &other )
  {
    const_iterator oIt = other.begin();
    for( iterator it = begin(); it != end(); ++it, ++oIt )
      *it -= *oIt;
    return *this;
  }

  Vector &operator*= ( const double s )
  {
    for( double &v : data_ )
      v *= s;
    return *this;
  }

  Vector &operator/= ( const double s )
  {
    for( double &v : data_ )
      v /= s;
    return *this;
  }

  double operator* ( const Vector &other ) const
  {
    double value = 0.;
    const_iterator oIt = other.begin();
    for( const_iterator it = begin(); it != end(); ++it, ++oIt )
      value += ( *it ) * ( *oIt );
    return value;
  }

  // add multiple of another vector to receiver object
  void axpy ( double s, const Vector &other )
  {
    const_iterator oIt = other.begin();
    for( iterator it = begin(); it != end(); ++it, ++oIt )
      *it += s * (*oIt);
  }

  int size () const { return data_.size(); }

  double &operator[] ( const int i ) { return data_[ i ]; }
  const double &operator[] ( const int i ) const { return data_[ i ]; }

protected:
  Data data_;
};


// print vector
std::ostream &operator<< ( std::ostream &out, const Vector &v )
{
  for( const auto &e : v )
    out << e << " ";
  return out;
}

#endif // #ifndef VECTOR_HH
