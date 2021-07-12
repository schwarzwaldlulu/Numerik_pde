//  mostly sample code from the tutor


#ifndef MATRIX_HH
#define MATRIX_HH

#include "vector.hh"


class Matrix
{
  typedef std::vector< Vector > Data;

public:
  typedef Data::iterator iterator;
  typedef Data::const_iterator const_iterator;

  Matrix ( int rows, int cols, double value = 0 )
    : data_( rows, Vector( cols, value ) ),
      rows_( rows ),
      columns_( cols )
  {}

  Matrix &operator= ( double s )
  {
    for( Vector &v : data_ )
      v = s;
    return *this;
  }

  iterator begin () { return data_.begin(); }
  const_iterator begin () const { return data_.begin(); }

  iterator end () { return data_.end(); }
  const_iterator end () const { return data_.end(); }

  void operator() ( const Vector &a, Vector &b ) const
  {
    auto k = b.begin();
    for( const_iterator row = begin(); row != end(); ++row, ++k )
      (*k) = ( *row ) * a;
  }

  int rows () const { return rows_; }
  int columns () const { return columns_; }

  Vector &operator[] ( const int i ) { return data_[ i ]; }
  const Vector &operator[] ( const int i ) const { return data_[ i ]; }

protected:
  std::vector< Vector > data_;
  int rows_;
  int columns_;
};


// print position
std::ostream &operator<< ( std::ostream &out, const Matrix &v )
{
  for( const auto &e : v )
    out << e <<"\n";
  return out;
}
#endif // #ifndef MATRIX_HH
