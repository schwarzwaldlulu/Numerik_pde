#ifndef GLOBALREFINE_HH
#define GLOBALREFINE_HH

#include <vector>

#include "grid.hh"

using namespace std;

Grid globalRefine ( const Grid &coarseGrid )
{
    // vertices are still vertices in the refined grid. so we copy the old ones and add the new ones
    // triangles are no longer the old ones, so we initialite them as empty vector
    vector< Vertex > vertices( coarseGrid.vertices() );

    vector< Triangle > triangles;
    
    // use map to keep track of the new vertices to avoid repeating calculation 
    // key: index pair of the vertices of the edges,
    // value: the new vertex at the middle point of the edge
    map< pair< int, int >, Vertex > newVertices;
    
    for ( const Triangle& T: coarseGrid.triangles() )
    {
        array< int, 3 > indices{ { T.vertices()[0].index(), T.vertices()[1].index(), T.vertices()[2].index() } };
        
        array< Vertex, dimension + 1 > midpts;
        
        for ( int i = 0; i < 3; ++i )
        {
            // build the opposite edge of the ith vertex
            int v = indices[ ( i + 1 ) % 3 ];
            int w = indices[ ( i + 2 ) % 3 ];
            pair< int, int > pr = make_pair( min( v, w ), max( v, w ) );

            map< pair< int, int >, Vertex >::const_iterator it = newVertices.find( pr );

            if ( it !=  newVertices.end() ) //the vertex has been checked in another triangle earlier
                midpts[ i ] = it -> second; // we just retrieve the vertex
                
            else
            {
                Position pos = ( vertices[ pr.first ].position() + vertices[ pr.second ].position() ) * 0.5;

                Vertex v( pos, vertices.size() );

                vertices.push_back( v );// update vertices

                midpts[ i ] = v; // fill array midpts

                newVertices[ pr ] = v; // update map newVertices to avoid calculating it again

            }
        }
        
        // add the 4 new triangles
        triangles.push_back( Triangle( array< Vertex, dimension + 1 >{{ T.vertices()[ 0 ], midpts[ 1 ], midpts[ 2 ] }}, triangles.size() ) );

        triangles.push_back( Triangle( {{ T.vertices()[ 1 ], midpts[ 0 ], midpts[ 2 ] }}, triangles.size() ) );

        triangles.push_back( Triangle( array< Vertex, dimension + 1 >{{ T.vertices()[ 2 ], midpts[ 0 ], midpts[ 1 ] }}, triangles.size() ) );

        triangles.emplace_back( array< Vertex, dimension + 1 >{{ midpts[ 0 ], midpts[ 1 ], midpts[ 2 ] }}, triangles.size() );

    }       
    
   return Grid( move( vertices ), move( triangles ) );
}


#endif // #ifndef GLOBALREFINE_HH
