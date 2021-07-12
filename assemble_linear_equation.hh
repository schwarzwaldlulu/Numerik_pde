#include "vector.hh"
#include "grid.hh"
#include "position.hh"
#include "sparsematrix.hh"

// make the mass matrix, stiffness matrix and right-hand side of the linear system using lagrange linear basis
// for solving reactions-diffusion-equation

void mass_matrix( const Grid& grid, SparseMatrix& M )
{
    for ( const Triangle& T: grid.triangles() )
        for ( const Vertex& v: T.vertices() )
            for ( const Vertex& w: T.vertices() )
            {
                if ( v != w )
                    //using edge midpoint quadrature, where phi_i * phi_j != 0 only at the midpoint of the edge connecting vertex i and j
                    M[ v.index() ][ w.index() ] += T.volume() * ( 1. / 3. ) * 0.5 * 0.5;
                else
                    // phi_i = 0.5 at midpoints of the 2 edges going out of vertex i
                    M[ v.index() ][ w.index() ] += T.volume() * ( 1. / 3. ) * 0.5;
            }

}

void stiffness_matrix( const Grid& grid, SparseMatrix& S, double lambda )
{
    for ( const Triangle& T: grid.triangles() )
        for ( const Vertex& v: T.vertices() )
            for ( const Vertex& w: T.vertices() )
            {
                // appromimate integral_T nabla phi_i * nabla phi_j
                S[ v.index() ][ w.index() ] += lambda * T.volume() * ( T.normal( v ) * T.normal( w ) );
            }
}

// calculate the right-hand side of the linear system b_i = /integral f * phi_i
template < typename F >
void rhs_centroid_quadrature( const Grid& g, F f, Vector& b )
{
    if( b.size() != g.numVertices() )
    {
      Vector b_( g.numVertices() );
      b = b_;
    }
    for ( const Triangle& T: g.triangles() )
        for ( const Vertex& v: T.vertices() )
        {
            // phi_i(centroid) = 1/3
            b[ v.index() ] += T.volume() * ( 1. / 3. ) * f( T.baryCenter() );
        }
}

template < typename F >
void rhs_corner_quadrature( const Grid& g, F f, Vector& b )
{
    if( b.size() != g.numVertices() )
    {
      Vector b_( g.numVertices() );
      b = b_;
    }
    for( const Triangle& T: g.triangles() )
        for( const Vertex& v: T.vertices() )
            b[ v.index() ] += T.volume() + ( 1. / 3. ) * f( v.position() );
}

template < typename F >
void rhs_edge_midpoint_quadrature( const Grid& g, F f, Vector& b, double lambda = 0 )
{
    if( b.size() != g.numVertices() )
    {
      Vector b_( g.numVertices() );
      b = b_;
    }

    for ( const Triangle& T: g.triangles() )
        for ( const Vertex& v: T.vertices() )
            for ( const Vertex& w: T.vertices() )
            {
                if ( v != w )
                {
                    Position midpoint = ( v.position() + w.position() ) * 0.5;
                    b[ v.index() ] += ( 1. / 3. ) * T.volume() * f( midpoint, lambda ) * 0.5;
                }
            }
}
