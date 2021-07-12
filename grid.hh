#ifndef GRID_HH
#define GRID_HH

#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <utility>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>

#include "position.hh"

using namespace std;

struct Vertex
{

  // default constructor, OTHERWISE code like array<Vertex, 3> would be illegal, would get messages like
  // error: use of deleted function 'std::array<Vertex, 3>::array()'

  Vertex () : pos_( 0 ), index_( -1 ) {}

  Vertex ( const Position& pos, int index ) : pos_( pos ), index_( index ) {}

  Vertex ( const Vertex& other ) = default;

  Vertex& operator= ( const Vertex& other ) = default;

  int index () const { return index_; }
  const Position& position () const { return pos_; }

  // overloading == and != for user defined object 'T' is necessary to use the algorithm find(inputIt first, inputIt last, const T& toFind)
  bool operator== ( const Vertex& other ) const { return other.index() == index(); }
  bool operator!= ( const Vertex& other ) const { return other.index() != index(); }

private:
  Position pos_;
  int index_;
};


/*
for the calculation of H_1 error, we need the gradients of lagrange basis phi_i, defined by 
  phi_i(x_i) = 1, phi_i(x_j) = 0 if j != i 
set phi_0 = a*x + b*y + c, then gradient of phi_0 = (a, b).  
set 3 corners of the triangle to pt0 = (x_0, y_0), pt1 = (x_1, y_1), pt2 = (x_2, y_2),
then by definition 
phi_0(x_0, y_0) = 1, phi_0(x_1, y_1) = 0, phi_0(x_2, y_2) = 0, 

since pt1-pt0 and pt2-pt0 are linearly independent,
we can solve it to get 
a = (y_1 - y_2) / ( (y_2 - y_0) * (x_1 - x_0) - (y_1 - y_0) * (x_2 - x_0) )
b = (x_2 - x_1) / ( (y_2 - y_0) * (x_1 - x_0) - (y_1 - y_0) * (x_2 - x_0) ), where
the numerator is perpendicular to pt2-pt1 and the denominator is the det of matrix composed of (pt1-pt0, pt2-pt0)
similarly gradient of phi_1 and phi_2. 
note that det[pt2-pt1,pt0-pt1] = det[pt1-pt0, pt2-pt0] = det[pt0-pt2, pt1-pt2]

geometrically, gradient should be perpendicular to pt2-pt1, as 
the lines parellel to pt2-pt1 are contours of phi_0
*/

struct Triangle
{
    
  Triangle ( const array< Vertex, dimension + 1 >& vertices, int index )
    : vertices_( vertices ), index_( index )
  {
    computeGeometry();
  }

  Triangle ( const vector< Vertex >& vertices, const array< int, dimension +1 >& indices, int index )
    : vertices_ { { vertices[ indices[ 0 ] ], vertices[ indices[ 1 ] ], vertices[ indices[ 2 ] ] } }, index_( index )
  {
    computeGeometry();
  }
  
  Triangle ( const Triangle& other ) = default;
  
  const array< Vertex, dimension + 1 >& vertices () const { return vertices_; }

  const Position& normal ( const Vertex& v ) const
  {
    return normal_[ localIndex( v ) ];
  }

  int localIndex ( const Vertex& v ) const
  {
    //std::find( inputIt first, inputIt last, const T& toFind ) returning iterator to the 1st elem found, defined in <algorithm>
    //std::distance( inputIt first, inputIt last ), defined in <iterator>
    return distance( vertices_.begin(), find( vertices_.begin(), vertices_.end(), v ) );
  }
  
  int index () const { return index_; }

  const Position& baryCenter () const { return baryCenter_; }

  double volume () const { return volume_; }

protected:
  // compute baryCenter, volume and the gradients of the 3 lagrange basis functions
  void computeGeometry ()
  {
    for( const Vertex& v : vertices_ )
      baryCenter_ += 1. / double( dimension + 1 ) * v.position();
      //baryCenter_.axpy( 1. / double( dimension + 1 ), v.position() ) works as well;

    //calculate volume as half of the abs of det of matrix composed of sides of the triangle
    array< Position, 2 > A {{ vertices_[ 1 ].position() - vertices_[ 0 ].position(), vertices_[ 2 ].position() - vertices_[ 0 ].position()}};
    const double detA = A[ 0 ][ 0 ] * A[ 1 ][ 1 ] - A[ 0 ][ 1 ] * A[ 1 ][ 0 ];

    volume_ = abs( detA )* 0.5;

    // compute normals with length of edge
    for( int i = 0; i < dimension +1; ++i )
    {
      Position normal( vertices_[ ( i + 1 ) % ( dimension +1 ) ].position()
                       - vertices_[ ( i + 2 ) % ( dimension +1 ) ].position() );
      normal /= detA;

      normal_[ i ][ 0 ] = -normal[ 1 ];
      normal_[ i ][ 1 ] = normal[ 0 ];
    }
  }

private:
  array< Vertex, dimension + 1 > vertices_;
  array< Position, dimension + 1 > normal_;
  int index_;
  Position baryCenter_;
  double volume_;
};


class Grid
{
public:
  Grid () {}

  Grid ( const string& filename )
  {
    ifstream file( filename.c_str() );

    if( !file.is_open() )
    {
      cerr<< "Can not construct grid from file: " << filename << endl;
      abort();
    }

    cout<< "Constructing grid from file: " << filename << endl;
    
    string line;

    while( getline( file, line ) )
    {
      // skip comments
      line = line.substr( 0, line.find_first_of( "//" ) );

      // Vertices block
      if( line.find( "Vertices" ) != string::npos )
        readVertexBlock( file );

      // Triangle block
      if( line.find( "Triangles" ) != string::npos )
        readTriangleBlock( file );
    }
    
    file.close();

    computeBoundaryVertices();
  }

  Grid ( const vector< Vertex >& vertices, const vector< Triangle >& triangles )
    : vertices_( vertices ), triangles_( triangles )
  {
    computeBoundaryVertices();
  }
  
  Grid ( const Grid &grid ) = default;

  Grid& operator= ( const Grid & ) = default;

  int numVertices () const { return vertices_.size(); }
  int numTriangles () const { return triangles_.size(); }

  const vector< Vertex >& vertices () const { return vertices_; }

  const vector< Triangle >& triangles () const { return triangles_; }

  const vector< Vertex >& boundaryVertices () const { return boundaryVertices_; }
  
  double edgelength( int i , int j ) const
  {
    Position edge = vertices_[i].position() - vertices_[j].position();
    return sqrt( edge * edge );
  }

  // gridwidth is defined as the max length of edges. In the implementation we go through all the edges, calculate its length, fill the map< > edgeLengths_
  // and return the max length
  // the part in the /*..*/ is the sample code from the tutor, which is slow if grid is large because almost every edgelength is calculated twice 
  /*
  double gridWidth ( const Grid &grid )
  {
    double width = numeric_limits< double >::min();
    for( const Triangle &T : grid.triangles() )
      for( const Vertex &v : T.vertices() )
        width = max( width, sqrt( T.normal(v) * T.normal( v ) ) * 2.0 * T.volume() );
    return width;
  }
*/
  double gridWidth( ) const
  {
    double d = numeric_limits< double >::min();

    set< pair< int, int > > set_of_edges;

    for( const Triangle& T: triangles_ )
      for( const Vertex& v: T.vertices() )
        for( const Vertex& w: T.vertices() )
        {
          if( v.index() < w.index() )
          {
            pair< int, int > pr = make_pair( v.index(), w.index() );
            if( set_of_edges.find( pr ) == set_of_edges.end() )
            {
              double length = edgelength( pr.first, pr.second );
              edgeLengths_[ pr ] = length;
              set_of_edges.insert( pr );
              d = max( d, length );
            }
          }
        }

    return d;
  }

protected:
  void readVertexBlock ( ifstream &in )
  {
    string line;
    
    while( getline( in, line ) )
    {
      // return if block is closed
      if( line.find( "#" ) != string::npos )
        break;

      stringstream lineStream;
      lineStream << line;

      Position p;
      lineStream >> p;
      int index = vertices_.size();
      vertices_.emplace_back( p, index );
      lineStream.clear();
    }
  }

  void readTriangleBlock ( ifstream& in )
  {
    string line;
    
    while( getline( in, line ) )
    {
      // return if block is closed
      if( line.find( "#" ) != string::npos )
        break;

      stringstream lineStream;
      lineStream << line;

      array< int, dimension + 1 > c;
      for( int& k : c )
        lineStream >> k;

      int index = triangles_.size();
      triangles_.emplace_back( vertices_, c, index );

      lineStream.clear();
    }
  }

  // for exercise 7.2
  // boundary vertices are vertices on exterior edges, i.e. edges only contained in one triangle
  void computeBoundaryVertices ()
  {
    set< pair< int, int > > boundaryEdges;

    // insert all boundary edges into set
    for( const Triangle &T : triangles_ )
      for( const Vertex &v : T.vertices() )
        for( const Vertex &w : T.vertices() )
        {
          if( v.index() > w.index() )
          {
            auto insert = boundaryEdges.insert( make_pair( w.index(), v.index() ) );
            // if edge already exists we can remove it again, since it is an interior edge
            // note the return value of insert is pair< iterator, bool >
            if( !insert.second )
              boundaryEdges.erase( insert.first );
          }
        }
          

    // the part in /*..*/ is the sample solution from the tutor which is far slow in the case of large grid
    // it's much faster using set.insert() instead of sort on vector, since
    // the complexity of sort is O( n * logn ), while 
    // the complexity of set.insert is O( logn ) in the worst case
    
    /*
    // insert all boundary vertices into vector
    for( auto &e : edges )
    {
      boundaryVertices_.push_back( vertices_[ e.first ] );
      boundaryVertices_.push_back( vertices_[ e.second ] );
    }

    // sort boundary vertices
    std::sort( boundaryVertices_.begin(), boundaryVertices_.end(), 
        [] ( const Vertex &a, const Vertex &b ) { return a.index() < b.index(); } );

    // remove duplicates
    boundaryVertices_.resize( std::distance( boundaryVertices_.begin(), std::unique( boundaryVertices_.begin(), boundaryVertices_.end() ) ) );
    */
    
    set< int > edgeSet;
    for ( auto& e: boundaryEdges )
    {
      edgeSet.insert( e.first );
      edgeSet.insert( e.second );
    }
    
    // using for( int& i: edgeSet ) instead of const int resulting in
    // error: binding reference of type ‘int&’ to ‘const int’ discards qualifiers
    for( const int i: edgeSet )
      boundaryVertices_.push_back( vertices_[ i ] );
  }

private:
  vector< Vertex > vertices_;
  vector< Triangle > triangles_;
  vector< Vertex > boundaryVertices_;
  
  // mutable variables can be modified by const member functions
  mutable map< pair< int, int >, double > edgeLengths_;
};

#endif // #ifndef GRID_HH
