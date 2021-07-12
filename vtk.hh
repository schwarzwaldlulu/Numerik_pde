//  mostly sample code from the tutor

#ifndef VTK_HH
#define VTK_HH

#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

#include "position.hh"
#include "vector.hh"
#include "grid.hh"

template< class Data >
struct VTKDataType;

template<>
struct VTKDataType< float >
{
  static std::string type () { return "Float32"; }
  static int components () { return 1; }

  static std::string toString ( const float &value )
  {
    std::ostringstream ret;
    ret << std::scientific << std::setprecision( 8 ) << value;
    return ret.str();
  }
};

template<>
struct VTKDataType< double >
{
  static std::string type () { return "Float64"; }
  static int components () { return 1; }

  static std::string toString ( const double &value )
  {
    std::ostringstream ret;
    ret << std::scientific << std::setprecision( 16 ) << value;
    return ret.str();
  }
};

template<>
struct VTKDataType< std::int32_t >
{
  static std::string type () { return "Int32"; }
  static int components () { return 1; }
  static std::string toString ( const std::int32_t &value ) { return std::to_string( value ); }
};

template<>
struct VTKDataType< std::uint8_t >
{
  static std::string type () { return "UInt8"; }
  static int components () { return 1; }
  static std::string toString ( const std::uint8_t &value ) { return std::to_string( value ); }
};

template<>
struct VTKDataType< Position >
{
  typedef VTKDataType< double > VTKDouble;

  static std::string type () { return VTKDouble::type(); }
  static int components () { return 3; }

  static std::string toString ( const Position &value )
  {
    return VTKDouble::toString( value[ 0 ] ) + " " + VTKDouble::toString( value[ 1 ] ) + " " + VTKDouble::toString( double( 0 ) );
  }
};



// VTKWriter
// ---------

class VTKWriter
{
  typedef std::vector< std::pair< std::string, std::reference_wrapper< const Vector > > > VectorList;

public:
  VTKWriter ( const Grid &grid ) : grid_( grid ) {}

  void addVertexData ( const std::string &name, const Vector &data ) { pointData_.emplace_back( name, data ); }
  void addTriangleData ( const std::string &name, const Vector &data ) { cellData_.emplace_back( name, data ); }

  void write ( const std::string &name ) const
  {
    std::ofstream vtu( name );
    vtu << "<?xml version=\"1.0\"?>" << std::endl;
    vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
        << "byte_order=\"LittleEndian\">" << std::endl;
    vtu << "  <UnstructuredGrid>" << std::endl;
    vtu << "    <Piece NumberOfPoints=\"" << grid_.numVertices() << "\" "
        << "NumberOfCells=\"" << grid_.numTriangles() << "\">" << std::endl;

    writeDataVectors( vtu, "PointData", pointData_ );
    writeDataVectors( vtu, "CellData", cellData_ );

    vtu << "      <Points>" << std::endl;
    writeCoordinates( vtu );
    vtu << "      </Points>" << std::endl;

    vtu << "      <Cells>" << std::endl;
    writeConnectivity( vtu );
    writeOffsets( vtu );
    writeTypes( vtu );
    vtu << "      </Cells>" << std::endl;

    vtu << "    </Piece>" << std::endl;
    vtu << "  </UnstructuredGrid>" << std::endl;
    vtu << "</VTKFile>" << std::endl;
  }

private:
  static void writeDataVectors ( std::ostream &vtu, const std::string &tag, const VectorList &vectors )
  {
    vtu << "      <" << tag;
    if( !vectors.empty() )
    {
      vtu << " Scalars=\"";
      int count( 0 );
      for( const auto &v : vectors )
        vtu << std::string( (count++) > 0 ? "," : "" ) << v.first;
      vtu << "\"";
    }
    vtu << ">" << std::endl;

    for( const auto &v : vectors )
      writeDataArray< double >( vtu, v.first, v.second.get() );

    vtu << "      </" << tag << ">" << std::endl;
  }

  void writeCoordinates ( std::ostream &vtu ) const
  {
    std::vector< Position > points( grid_.numVertices() );
    for( const Vertex &v : grid_.vertices() )
      points[ v.index() ] = v.position();
    writeDataArray( vtu, "Coordinates", points );
  }

  void writeConnectivity ( std::ostream &vtu ) const
  {
    std::vector< int > connectivity( (dimension+1)*grid_.numTriangles() );
    for( const Triangle &t : grid_.triangles() )
    {
      for( int i = 0; i <= dimension; ++i )
        connectivity[ (dimension+1)*t.index() + i ] = t.vertices()[ i ].index();
    }
    writeDataArray( vtu, "connectivity", connectivity );
  }

  void writeOffsets ( std::ostream &vtu ) const
  {
    const int size = grid_.numTriangles();
    std::vector< std::int32_t > offsets( size );
    for( int i = 0; i < size; ++i )
      offsets[ i ] = (dimension+1)*(i+1);
    writeDataArray( vtu, "offsets", offsets );
  }

  void writeTypes ( std::ostream &vtu ) const
  {
    std::vector< std::uint8_t > types( grid_.numTriangles(), 5 );
    writeDataArray( vtu, "types", types );
  }

  template< class Data >
  static void writeDataArray ( std::ostream &vtu, const std::string &name, const std::vector< Data > &data )
  {
    vtu << "        <DataArray type=\"" << VTKDataType< Data >::type() << "\" "
        << "NumberOfComponents=\"" << VTKDataType< Data >::components() << "\" "
        << "Name=\"" << name << "\" format=\"ascii\">" << std::endl;
    int linelength = 0;
    for( const Data &value : data )
    {
      const std::string s = VTKDataType< Data >::toString( value );
      const int sz = s.size();
      if( linelength + sz > 79 )
      {
        vtu << std::endl;
        linelength = 0;
      }
      std::string indent( linelength < 10 ? "          " : " " );
      vtu << indent << s;
      linelength += (sz + indent.size());
    }
    if( linelength > 0 )
      vtu << std::endl;
    vtu << "        </DataArray>" << std::endl;
  }

  const Grid &grid_;
  VectorList pointData_, cellData_;
};

#endif // #ifndef VTK_HH
