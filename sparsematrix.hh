// sample code from tutor

#ifndef SPMATRIX_HH
#define SPMATRIX_HH

#include <map>
#include <vector>

#include "vector.hh"


class SparseMatrix
{

  typedef std::array< unsigned int, 2 > IndexPair;
  typedef std::map< IndexPair, double > Data;

  typedef typename Data::iterator iterator;
  typedef typename Data::const_iterator const_iterator;

public:

  struct Row
  {
    Row ( SparseMatrix &matrix, unsigned int row ) : matrix_( matrix ), row_( row ) {}

    double& operator[] ( unsigned int col )
    { 
      return data().insert( std::make_pair( IndexPair{{ row_, col }}, 0.0 ) ).first->second;  
    }

    double operator[] ( unsigned int col ) const
    {
      const const_iterator pos = data().find( IndexPair{{ row_, col }} );
      return (pos == data().end() ? 0.0 : pos->second);
    } 

    unsigned int size () const { return matrix_.columns(); }

    void clear () { data().erase( data().lower_bound( IndexPair{{ row_, 0 }} ), data().lower_bound( IndexPair{{ row_, size() }} ) ); }

  protected:
    Data &data () { return matrix_.data_; }
    const Data &data () const { return matrix_.data_; }
  private:
    SparseMatrix &matrix_;
    const unsigned int row_;
  };

  struct ConstRow
  {
    ConstRow ( const SparseMatrix &matrix, unsigned int row ) : matrix_( matrix ), row_( row ) {}

    double operator[] ( unsigned int col ) const
    {
      const const_iterator pos = matrix_.data_.find( IndexPair{{ row_, col }} );
      return (pos == matrix_.data_.end() ? 0.0 : pos->second);
    } 

    unsigned int size () const { return matrix_.columns(); }

  private:
    const SparseMatrix &matrix_;
    const unsigned int row_;
  };


  explicit SparseMatrix ( const int rows, const int cols ) : 
    rows_( rows ),
    cols_( cols )
  {}

  void operator() ( const Vector &a, Vector &b ) const
  {
    b = 0.0;
    for( const std::pair< IndexPair, double > &entry : data_ )
      b[ entry.first[ 0 ] ] += entry.second * a[ entry.first[ 1 ] ];
  }

  Row operator[] ( unsigned int r ) { return Row( *this, r ); }
  ConstRow operator[] ( unsigned int r ) const { return ConstRow( *this, r ); }

  unsigned int rows () const { return rows_; }
  unsigned int columns () const { return cols_; }

  void clear () { data_.clear(); }

private:
  Data data_;
  const unsigned int rows_;
  const unsigned int cols_;
};

#endif
