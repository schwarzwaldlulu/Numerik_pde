#include <complex>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif

static const int dimension = 2;

#include "cg.hh"
#include "grid.hh"
#include "sparsematrix.hh"
#include "position.hh"
#include "vector.hh"
#include "vtk.hh"
#include "globalrefine.hh"
#include "assemble_linear_equation.hh"

using namespace std;


// use FEM to calculate the weak solution to a pde. examples are 
// poisson's equation: 
// laplace u = f in G, u = g on boudary G
// for any function phi from sobolev space H^1(omega) the weak solution u_h satisfies:
// integral_omega (gradient u_h * gradient phi) + integral_omega (f*phi) = 0

// reaction-diffusion-equation: 
// u - lambda*laplace u = f in G, u = g on boundary G
// its weak solution u_h satisfies:
// integral_omega (u_h * phi) + lambda * integral_omega (gradient u_h * gradient phi) = integral_omega (f*phi)

// abstract problem definition
struct Problem 
{
  // exact solution
  virtual double u ( const Position& x ) const  = 0;
  
  //boundary condition
  virtual double g ( const Position& x ) const  = 0;
  
  // M = 1 if integral_omega (u_h * phi) is in the weak form, 0 otherwise
  virtual int M () const = 0;
  
  virtual Position gradient ( const Position& x ) const = 0;
  
  virtual double f ( const Position& x, double lambda ) const = 0;
};


struct Sinus
: public Problem
{
  double u ( const Position& x ) const
  {
    return sin( 2.0 * M_PI * x[0] ) * sin( 2.0 * M_PI * x[1] );
  }
  
  double g ( const Position& x ) const
  {
    return u( x );
  }

  int M () const { return 1; }
  
  Position gradient ( const Position& x ) const
  {
    Position grad( 0 );
    grad[ 0 ] = 2.0 * M_PI * cos( 2.0 * M_PI * x[0] ) * sin( 2.0 * M_PI * x[1] );
    grad[ 1 ] = 2.0 * M_PI * cos( 2.0 * M_PI * x[1] ) * sin( 2.0 * M_PI * x[0] );
    return grad;
  }

  // f = u - lambda * laplace u, here laplace u = -8.0 * M_PI * M_PI * u
  double f ( const Position& x, double lambda )  const
  {
    return ( 1.0 + lambda * 8.0 * M_PI * M_PI ) * u( x );
  }
};


struct SaddelPoint
: public Problem
{
  double u ( const Position& x ) const
  {
    return x[ 0 ] * x[ 0 ] - x[ 1 ] * x[ 1 ];
  }

  double g ( const Position& x ) const
  {
    return u( x );
  }

  int M () const { return 1; }
  
  Position gradient ( const Position& x ) const
  {
    Position grad( x );
    grad[ 0 ] *= 2.0;
    grad[ 1 ] *=-2.0;
    return grad;
  }

  // f = u - lambda*laplace u, here laplace u = 0, so f = u
  double f ( const Position& x, double lambda )  const
  {
    return u(x);
  }
};

//returns the phase angle of the complex number x + y*i, return value restricted in [0, 2 * M_PI]
double phaseAngle ( double x, double y )
{
  // return value of arg(complex<..>(x,y)) is atan2(y, x), in [-M_PI, M_PI], in radians
  double phi = arg( complex< double >( x,y ) );
  
  if( y < 0 ) 
      phi += 2. * M_PI;
  
  return phi;
}

// for the question: why are manhole covers round?
// alpha = ( 3 * M_PI ) / 2,  omega = {(x, y) ∈ R ∶ ∣x∣ + ∣y∣ < 1 , arctan(y/x) < alpha }, 
// homogeneous poisson equation:
// -laplace u = 0 in omega, u = g on boudary omega, g(x, y) = r(x, y)^lambda * sin( lambda * arctan(y/x) )
// where lambda = M_PI / alpha = 2./3., r( x, y ) = sqrt( x^2 + y^2 )

struct ReentrantCorner
: public Problem
{
  double lambda = 180. / 270.;// change lambda = 2./3. to 180./270 to falicitate the calculation later.

  double u ( const Position& x ) const
  {
    double phi = phaseAngle( x[ 0 ], x[1] );
    
    double r2 = x * x; // x^2 + y^2
    
    return pow( r2, lambda * 0.5 ) * sin( lambda *phi ); // ||x||^lambda * sin( lambda*phi)
  }

  double g ( const Position& x ) const
  {
    return u( x );
  }

  int M () const { return 0; }

  // partial derivative of u w.r.t. x is:
  // lambda * x * (x^2 + y^2)^(lambda/2. - 1) * sin(lambda * arctan(y/x)) - lambda * y * (x^2 + y^2)^(lambda/2. - 1) * cos(lambda * arctan(y/x))
  // partial derivatice of u w.r.t. y is ( similar to pd w.r.t. x due to symmetry ):
  // lambda * y * (x^2 + y^2)^(lambda/2. - 1) * sin(lambda * arctan(y/x)) + lambda * x * (x^2 + y^2)^(lambda/2. - 1) * cos(lambda * arctan(y/x))
  Position gradient ( const Position& x ) const
  {
    double phi = phaseAngle( x[ 0 ],x[ 1 ]); 
    double r2 = x * x;
    Position grad( x ), grad2( x );

    grad *= lambda * pow( r2, lambda * 0.5 - 1.0 ) * sin( lambda * phi );
    grad2 *= lambda * pow( r2, lambda * 0.5 - 1.0 ) * cos( lambda * phi );

    grad[ 0 ] -= grad2[ 1 ];
    grad[ 1 ] += grad2[ 0 ];
    
    return grad;
  }

  // f = -laplace u = 0
  double f ( const Position& x, double l )  const
  {
    return 0.;
  }
};

// lagrange interpolation for u is Iu = sum_i u_i * phi_i, sum over all basis function phi_i
// and the coeffients u_i = u( pos_i ), i.e. the function value at vertex i
// calculate the coeffients
template< class Vertices, class Function >
void lagrangeInterpolation ( const Vertices &vertices, const Function &f, Vector &x )
{
  for( const Vertex &v : vertices )
    x[ v.index() ] = f( v.position() );
}


// assemble the linear system for solving pde, i.e. build mass matrix, stiffness matrix and right-hand side, where
// mass matrix M_ij = integral_omega phi_i * phi_j
// stiffness matrix S_ij = integral_omega gradient( phi_i ) * gradient( phi_j )
// so the left-hand side matrix is: A_ij = M_ij + lambda * S_ij, or delta_ij, if vertex i is on the boundary
// right-hand side b_i = integral_omega f*phi_i or g_i if vertex i is on the boudary
// bndData = 1 means boundary does not vanish

void assemble( const Grid& grid, const Problem& problem, SparseMatrix& A, Vector& b, Vector& u, double lambda, int bndData = 1 )
{
    A.clear();
    
    stiffness_matrix( grid, A, lambda ); // build lambda * stiffness matrix 

    if( problem.M() == 1 ) // integral_omega (u_h * phi) is in the weak form
      mass_matrix( grid, A ); // add mass matrix to stiffness matrix

    // signature of the function is:
    // template < typename F >
    // void rhs_edge_midpoint_quadrature( const Grid& g, F f, Vector& b, double lambda = 0 )
    rhs_edge_midpoint_quadrature( grid, [ &problem ] ( const Position& pos, double lambda ) -> double { return problem.f( pos, lambda ); }, b, lambda );

    for ( const Vertex& v: grid.boundaryVertices() )
    {
        A[ v.index() ].clear();
        A[ v.index() ][ v.index() ] = 1;
    }

    if ( bndData == 1 ) // boundary does not vanish
    {
        //b_i * = g( p_i ) if p_i is on the boundary, b_i otherwise
        lagrangeInterpolation( grid.boundaryVertices(), [ &problem ] ( const Position& pos) -> double { return problem.g( pos ); }, b );
        // set u_i = g( p_i ) if p_i is on the boudary as A_ii = 1, A_ij = 0 for j != i, then u_i = A_i * b = b_i
        lagrangeInterpolation( grid.boundaryVertices(), [ &problem ] ( const Position& pos ) -> double { return problem.g( pos ); }, u );
    }

}

// || f-g ||_l2 is defined as (integral_omega (f-g)^2 )^(1./2.) 
// use edge middlepoint quadrature to calculate l2 error
// x contains the coeffients of the lagrange interpolation function to u
double l2error ( const Grid& grid, const Problem& problem, const Vector& x )
{
  double error = 0.;

  for( const Triangle &triangle : grid.triangles() )
  {
    for( int i = 0; i < 3; ++i )
    {
      const Vertex &v = triangle.vertices()[ i ];
      const Vertex &w = triangle.vertices()[ (i+1) % 3 ];
      Position midPoint = ( v.position() + w.position() ) * 0.5;
      
      double localError = problem.u( midPoint ) - 0.5* ( x[ v.index() ] + x[ w.index() ]);  
      
      error += 1./ 3. * localError * localError * triangle.volume();
    }
  }

  return sqrt( error );
}


// || f-g ||_h1 is defined as 
//(integral_omega (f-g)^2 + integral_omega ( nabla(f-g) * nabla(f-g) ) )^(1./2.) 
// use edge middlepoint quadrature to calculate h1 error
// x contains the coeffients of the lagrange interpolation function to u
double h1error ( const Grid &grid, const Problem &problem, const Vector &x )
{
  double error = 0.;

  for( const Triangle &triangle : grid.triangles() )
    for( int i = 0; i < 3; ++i )
    {
      const Vertex &v = triangle.vertices()[ i ];
      const Vertex &w = triangle.vertices()[ (i+1) % 3 ];
      Position midPoint = ( v.position() + w.position() ) * 0.5;
      
      double localError = problem.u( midPoint ) - 0.5* ( x[ v.index() ] + x[ w.index() ]);  
      
      error += 1./ 3. * localError * localError * triangle.volume(); // the first part i.e l2 error

      Position gradError = problem.gradient( midPoint );
      
      for( const Vertex &p : triangle.vertices() )
        gradError += x[ p.index() ] * triangle.normal( p );

      error += 1./ 3. * ( gradError * gradError )* triangle.volume();
    }

  return sqrt( error );
}




void print ( const Grid &grid, const Vector &x, int count )
{
  VTKWriter vtk( grid );
  vtk.addVertexData( "solution", x );
  vtk.write( "output-" + to_string( count ) + ".vtu" );
}



int main ( int argc, char **argv )
{
  string filename = "../lmesh3";
  
  // print the message w.r.to. command line arguments, i.e. the following infomation should be given
  // meshfile, problem, lambda, refinements
  // refinements is the number of times the grid will be refined
  if( argc < 4 )
  {
    cout << "Program usage:"<< argv[ 0] <<" <meshfile> <Problem> <lambda> <refinments>" << endl;
    cout << "Using default mesh: "<< filename << endl;
    cout << "Problem '0':  u = sin( 2 * pi *x ) sin( 2 * pi* y) (default)" << endl;
    cout << "Problem '1':  v = x^2 - y^2" << endl;
    cout << "Problem '2':  Reentrant corner" << endl;
    cout << "lambda > 0 ( default = 1 )" << endl;
    cout << "refinments: number of EOC refinements (default 2 )" << endl;
    cout << endl;
  }

  if( argc > 1 )
      filename = argv[ 1 ];

  int probNumber = 0;
  if( argc > 2 )
      probNumber = stoi( argv[ 2 ] );// stoi() function takes a string as an argument and returns its value

  double lambda = 1;
  if( argc > 3 )
      lambda = stoi( argv[ 3 ] );

  int refinements = 2;
  if( argc > 4 )
      refinements = stoi( argv[ 4 ] );
  
  Problem* problem;
  switch( probNumber )
  {
    case 0:
      problem = new Sinus();
      break;
    case 1: 
      problem = new SaddelPoint();
      break;
    case 2:
      problem = new ReentrantCorner();
      lambda = 1.0;
      break;
    default:
      problem = new Sinus();
      break;
  }

  Grid grid( filename );
  vector< double > l2Error( refinements );
  vector< double > h1Error( refinements );
  vector< double > gw( refinements );//gridWidth

  for( int i =0 ; i < refinements; ++i )
  {
    const int numVertices = grid.numVertices();//number of basis functions of the function space
    Vector u( numVertices ), b( numVertices );
    SparseMatrix A( numVertices, numVertices );

    // fill matrix A and right-hand side b for solving pde
    assemble( grid, *problem, A, b, u, lambda );

    const double eps = 1e-12;
    // solve system with cg
    cg( A, b, u, eps, false );

    gw[ i ] = grid.gridWidth();
    l2Error[ i ] = l2error( grid, *problem, u );
    h1Error[ i ] = h1error( grid, *problem, u );

    print( grid, u, i );
    if( i != refinements -1 )
      grid = globalRefine( grid );
  }


  cout << "gridWidth: \t\tL2 error: \t\tEOC:   \t\tH1 error: \t\tEOC: "<< endl;
  cout << gw[ 0 ] << "\t\t" << l2Error[ 0 ] <<"\t\t---   \t\t" << h1Error[ 0 ]<<"\t\t---"<< endl;
  
  for( int i = 1; i< refinements; ++i )
  {
    const double h = log( gw[ i ] / gw[ i-1 ] );
    cout<< gw[ i ] << "\t\t" << l2Error[ i ] <<"\t\t"<< log( l2Error[ i ] / l2Error[ i -1 ] ) / h 
      << "\t\t" << h1Error[ i ] << "\t\t" << log( h1Error[ i ] / h1Error[ i-1 ] ) / h << endl;
  }

  return 0;
}
