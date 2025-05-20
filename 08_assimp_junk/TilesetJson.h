//
// Created on 2024/04/13
//

#ifndef TILESETJSON_H
#define TILESETJSON_H

#include <optional>
#include <list>

#include <json-c/json.h>
#include <json-c/printbuf.h>

#define SIZE(x) (sizeof(x) / sizeof(*x))
#ifndef ASSERT
# define ASSERT(exp) ((void) 0)
#endif

namespace cesiumjs
{
#if 0
}
#endif

class TilesetJson
{
#pragma region types
public:
    struct asset_t
    {
        std::string version;
        /*
          "extras": {
            "ion": {
                "georeferenced": true,
                "movable": true,
                "terrainId": 1
            }
        },
        */
        std::optional<bool> extra_ion_georeferenced;
    };

    class boundingVolume_t
    {
    public:
        virtual ~boundingVolume_t() {}
        //virtual void construct_from_points(
        //    double x0, double y0, double z0,
        //    double x1, double y1, double z1
        //) = 0;
    };

    /// <summary>
    /// Box
    /// The boundingVolume.box property is an array of 12 numbers that define an oriented
    /// bounding box in a right - handed 3 - axis( x, y, z ) Cartesian coordinate system where
    /// the z - axis is up.
    /// The first three elements define the x, y, and z values for the center of the box.
    /// The next three elements( with indices 3, 4, and 5 ) define the x - axis direction and half - length.
    /// The next three elements( indices 6, 7, and 8 ) define the y - axis direction and half - length.
    /// The last three elements( indices 9, 10, and 11 ) define the z - axis direction and half - length.
    /// </summary>
    class boundingVolumeBox_t : public boundingVolume_t
    {
    public:
        double box[12];

        boundingVolumeBox_t() { std::fill( box, box + 12, 0.0 ); }
        // general case
        boundingVolumeBox_t( std::initializer_list<double> const& l )
        {
            ASSERT( l.size() == SIZE( box ) );
            std::fill( box, box + SIZE( box ), 0.0 );
            std::copy( l.begin(), l.end(), box );
        }
        // Assume no rotation. The box is aligned with coordinate axis
        // The P0 and P1 are points on the bounding box diagonal
        boundingVolumeBox_t(
            double x0, double y0, double z0,
            double x1, double y1, double z1
        )
        {
#if 0
            box[0] = ( x0 + x1 ) / 2;
            box[1] = ( y0 + y1 ) / 2;
            box[2] = ( z0 + z1 ) / 2;
            box[3] = ( x1 - x0 ) / 2;
            box[4] = 0;
            box[5] = 0;
            box[6] = 0;
            box[7] = ( y1 - y0 ) / 2;
            box[8] = 0;
            box[9] = 0;
            box[10] = 0;
            box[11] = ( z1 - z0 ) / 2;
#else
            box[0] = ( x1 - x0 ) / 2;
            box[1] = ( y0 - y1 ) / 2;
            box[2] = ( z1 - z0 ) / 2;
            box[3] = ( x1 - x0 ) / 2;
            box[4] = 0;
            box[5] = 0;
            box[6] = 0;
            box[7] = ( y1 - y0 ) / 2;
            box[8] = 0;
            box[9] = 0;
            box[10] = 0;
            box[11] = ( z1 - z0 ) / 2;
#endif
        }
    };

    /// <summary>
    /// Region
    /// The boundingVolume.region property is an array of six numbers
    /// that define the bounding geographic region with latitude, longitude,
    /// and height coordinates with the order
    /// [west, south, east, north, minimum height, maximum height].
    /// Latitudes and longitudes are in the WGS 84 datum as defined in EPSG 4979 and are in radians.
    /// Heights are in meters above (or below) the WGS 84 ellipsoid.
    /// </summary>
    class boundingVolumeRegion_t : public boundingVolume_t
    {
    public:
        double region[6];
        boundingVolumeRegion_t()
        {
            std::fill( region, region + SIZE( region ), 0.0 );
        }

        boundingVolumeRegion_t( std::initializer_list<double> const& l )
        {
            ASSERT( l.size() == SIZE( region ) );
            std::copy( l.begin(), l.end(), region );
        }
    };

    struct transform_t
    {
        double matrix4x4[16];

        transform_t()
        {
            std::fill( matrix4x4, matrix4x4 + SIZE(matrix4x4), 0.0);
        }

        transform_t( std::initializer_list<double> const& t )
        {
            ASSERT( t.size() == SIZE( matrix4x4 ) );
            std::copy( t.begin(), t.end(), matrix4x4 );
        }

    };

    struct content_t
    {
        std::string uri;
    };

    /// <summary>
    /// a tile (?)
    /// </summary>
    class child_t
    {
    public:
        double geometricError;
        std::list<child_t> children;
        content_t content;
        std::unique_ptr<transform_t> transform; // optional?
        std::unique_ptr<boundingVolume_t> boundingVolume; // optional?

        child_t() : geometricError (0) {}

        child_t( double e, content_t const& c )
            : geometricError( e )
            , content( c )
        {
        }

        child_t( std::initializer_list<double> const& t )
            : geometricError( 0 )
            , transform( new transform_t( t ) )
        {
        }
    };

    class root_tile_t : public child_t
    {
    public:
        std::string refine;
    public:
        root_tile_t()
            : child_t{}
            , refine( "ADD" )
        {}

        root_tile_t( std::string const& r, std::initializer_list<double> const &t )
            : child_t( t )
            , refine(r)
        {
        }
    };
#pragma endregion
public:
    /// <summary>
    /// Construct the tileset object from existing json file
    /// </summary>
    /// <param name="filename"></param>
    TilesetJson( std::string const& filename )
    {
        m_root = json_object_from_file( filename.c_str() );
        if ( !m_root )
            throw std::runtime_error( "Document root is not valid" );
    }

    /// <summary>
    /// Construct an empty tileset object
    /// </summary>
    TilesetJson()
    {
        m_root = json_object_new_object();
    }

    virtual ~TilesetJson()
    {
        if ( m_root )
            json_object_put( m_root );
    }

    asset_t asset()
    {
        asset_t r{};
        json_object* asset_obj = 0;
        //    m_root.GetObj( "asset" );
        //r.version = asset_obj.GetInteger( "version", -1 );
        if ( json_object_object_get_ex( m_root, "asset", &asset_obj ) )
        {
            // TODO: get extras
            r.version = get_as_string( asset_obj, "version" );
        }
        return r;
    }

    /*
        "asset": {
            "extras": {
                "ion": {
                    "georeferenced": true,
                    "movable": true,
                    "terrainId": 1
                }
            },
            "version": "1.0"
        },
    */
    void asset( asset_t  const& ass )
    {
        auto o = json_object_new_object();
        json_object_object_add( o, "version", json_object_new_string( ass.version.c_str() ) );

        if ( ass.extra_ion_georeferenced.has_value() )
        {
            auto ion_o = json_object_new_object();
            json_object_object_add( ion_o, "georeferenced", json_object_new_boolean( ass.extra_ion_georeferenced.value() ) );
            json_object_object_add( ion_o, "movable", json_object_new_boolean( true ) );
            auto extras_o = json_object_new_object();
            json_object_object_add( extras_o, "ion", ion_o );
            json_object_object_add( o, "extras", extras_o );
        }

        json_object_object_add( m_root, "asset", o );
    }

    double geometricError()
    {
        return geometricError( m_root );
    }

    void geometricError( double d )
    {
        json_object_object_add( m_root, "geometricError", json_object_new_double_g18( d ) );
    }

    content_t content( json_object* o )
    {
        content_t r;
        //r.uri = o.GetString( "uri" );
        json_object* uri = 0;
        if ( json_object_object_get_ex( o, "uri", &uri ) )
            r.uri = json_object_get_string( uri );
        return r;
    }

    void child_from_obj( json_object* o, child_t& tile )
    {
        tile.geometricError = geometricError( o );
        tile.boundingVolume.reset( boundingVolume( o ) );

        json_object* co = 0;
        if ( json_object_object_get_ex( o, "content", &co ) )
            tile.content = content( co );
        json_object* tr = 0;
        if ( json_object_object_get_ex( o, "transform", &tr ) )
            tile.transform = transform( tr );
        

        json_object* kids_obj;

        if ( json_object_object_get_ex( o, "children", &kids_obj )
             && json_object_get_type( kids_obj ) == json_type_array )
        {
            size_t n = json_object_array_length(kids_obj);
            for ( int i = 0; i < n; ++i )
            {
                child_t ch;
                json_object *o = json_object_array_get_idx(kids_obj, i);
                child_from_obj( o, ch );
                tile.children.push_back( std::move( ch ) );
            }
        }
    }

    root_tile_t root()
    {
        json_object* root_obj;

        root_tile_t tile{};

        if ( json_object_object_get_ex( m_root, "root", &root_obj ) )
        {
            child_from_obj( root_obj, tile );
            json_object* refine_obj;
            if ( json_object_object_get_ex( root_obj, "refine", &refine_obj ) )
                tile.refine = json_object_get_string( refine_obj );
        }
        return tile;
    }

    /// <summary>
    /// Write to json object
    /// </summary>
    /// <param name="tile"></param>
    void root( root_tile_t const& tile )
    {
        auto ro = child_to_obj( tile );
        json_object_object_add( ro, "refine", json_object_new_string( tile.refine.c_str() ) );
        json_object_object_add( m_root, "root", ro );
    }

    json_object* transform_to_obj( transform_t const &tr )
    {
        auto tro = json_object_new_array();
        for ( auto d : tr.matrix4x4 )
            json_object_array_add( tro, json_object_new_double_g18( d ) );
        return tro;
    }

    void children_to_obj( std::list<child_t> const& children, json_object* tile_o )
    {
        auto arr_o = json_object_new_array();
        for ( auto const& cp : children )
            json_object_array_add( arr_o, child_to_obj( cp ) );
        json_object_object_add( tile_o, "children", arr_o );
    }

    json_object* child_to_obj( child_t const& tile )
    {
        auto ro = json_object_new_object();

        json_object_object_add(
            ro, "geometricError", json_object_new_double_g18( tile.geometricError ) );
        {
            if ( auto bounding_box = dynamic_cast<boundingVolumeBox_t*>( tile.boundingVolume.get() ) )
            {
                auto bvo = json_object_new_object();
                auto box_o = json_object_new_array();
                for ( auto d : bounding_box->box )
                    json_object_array_add( box_o, json_object_new_double_g18( d ) );
                json_object_object_add( bvo, "box", box_o );
                json_object_object_add( ro, "boundingVolume", bvo );
            }
            else if ( auto bounding_region = dynamic_cast<boundingVolumeRegion_t*>( tile.boundingVolume.get() ) )
            {
                auto bvo = json_object_new_object();
                auto box_o = json_object_new_array();
                for ( auto d : bounding_region->region )
                    json_object_array_add( box_o, json_object_new_double_g18( d ) );
                    json_object_object_add( bvo, "region", box_o );
                    json_object_object_add( ro, "boundingVolume", bvo );
            }
        }
        if ( !tile.content.uri.empty() )
            json_object_object_add( ro, "content", content_to_obj( tile.content ) );
        if ( auto* tr = tile.transform.get() )
            json_object_object_add( ro, "transform", transform_to_obj( *tr ) );
        if ( tile.children.size() )
            children_to_obj( tile.children, ro );

        return ro;
    }

    json_object* content_to_obj( content_t const& cont )
    {
        auto co = json_object_new_object();
        json_object_object_add( co, "uri", json_object_new_string( cont.uri.c_str() ) );
        return co;
    }

    static aiVector3D make_vector3D( double const* v )
    {
        return aiVector3D( ai_real( v[0] ), ai_real( v[1] ), ai_real( v[2] ) );
    }

    /// <summary>
    /// Construct a bounding volume wire-frame 
    /// </summary>
    /// <param name="vp"></param>
    aiScene* construct_scene( boundingVolume_t const* vp )
    {
        std::unique_ptr<aiScene> model{ new aiScene };

        if ( auto box = dynamic_cast<boundingVolumeBox_t const*>( vp ) )
        {
            auto mesh = bounding_box_mesh( box );

            model->mNumMeshes = 1;
            model->mMeshes = new aiMesh * [1];
            model->mMeshes[0] = mesh;

            aiNode* node_1 = new aiNode( "node[1]" );
            node_1->mMeshes = new unsigned int{};
            node_1->mNumMeshes = 1;

            aiNode* root = new aiNode( "node[0]" );
            node_1->mParent = root;

            root->mChildren = new aiNode * [1] { node_1 };
            root->mNumChildren = 1;
            root->mMeshes = 0;
            root->mTransformation = aiMatrix4x4(
#if 1
                // Y coordinate is elevation
                1, 0, 0, 0,
                0, 0, 1, 0,
                0, -1, 0, 0,
                0, 0, 0, 1
#else
                // Z coordinate is elevation
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1
#endif
            );
            model->mRootNode = root;
        }
        else
        {
            throw std::runtime_error( "Is not supported" );
        }
        return model.release();
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="box"></param>
    /// <returns></returns>
    aiMesh* bounding_box_mesh( boundingVolumeBox_t const* box )
    {
        // the mesh will represent a wire-frame box
        auto origin = make_vector3D( box->box );
        auto ort_x = make_vector3D( box->box + 3 );
        auto ort_y = make_vector3D( box->box + 6 );
        auto ort_z = make_vector3D( box->box + 9 );
#if 0
        CONSOLE( "origin:" << origin );
        CONSOLE( "ort_X:" << ort_x );
        CONSOLE( "ort_Y:" << ort_y );
        CONSOLE( "ort_Z:" << ort_z );
#endif
        aiMesh* mesh = new aiMesh;
        mesh->mPrimitiveTypes = aiPrimitiveType_LINE;

        auto vertices = new aiVector3D[8];
        // a bottom rectangle
        vertices[0] = origin + ort_x + ort_y - ort_z;
        vertices[1] = origin - ort_x + ort_y - ort_z;
        vertices[2] = origin - ort_x - ort_y - ort_z;
        vertices[3] = origin + ort_x - ort_y - ort_z;
        // a top rectangle
        vertices[4] = origin + ort_x + ort_y + ort_z;
        vertices[5] = origin - ort_x + ort_y + ort_z;
        vertices[6] = origin - ort_x - ort_y + ort_z;
        vertices[7] = origin + ort_x - ort_y + ort_z;

        mesh->mVertices = vertices;
        mesh->mNumVertices = 8;

        auto red = new aiColor4D[mesh->mNumVertices];
        for ( int i = 0; i != mesh->mNumVertices; ++i )
            red[i] = aiColor4D( 1, 0, 0, 0 );

        mesh->mColors[0] = red;

        aiFace* faces = new aiFace[12];

        faces[0].mNumIndices = 2;
        faces[0].mIndices = new unsigned int[2] { 0, 1 };

        faces[1].mNumIndices = 2;
        faces[1].mIndices = new unsigned int[2] { 1, 2 };

        faces[2].mNumIndices = 2;
        faces[2].mIndices = new unsigned int[2] { 2, 3 };

        faces[3].mNumIndices = 2;
        faces[3].mIndices = new unsigned int[2] { 3, 0 };

        faces[4].mNumIndices = 2;
        faces[4].mIndices = new unsigned int[2] { 4, 5 };

        faces[5].mNumIndices = 2;
        faces[5].mIndices = new unsigned int[2] { 5, 6 };

        faces[6].mNumIndices = 2;
        faces[6].mIndices = new unsigned int[2] { 6, 7 };

        faces[7].mNumIndices = 2;
        faces[7].mIndices = new unsigned int[2] { 7, 4 };

        faces[8].mNumIndices = 2;
        faces[8].mIndices = new unsigned int[2] { 0, 4 };

        faces[9].mNumIndices = 2;
        faces[9].mIndices = new unsigned int[2] { 1, 5 };

        faces[10].mNumIndices = 2;
        faces[10].mIndices = new unsigned int[2] { 2, 6 };

        faces[11].mNumIndices = 2;
        faces[11].mIndices = new unsigned int[2] { 3, 7 };

        mesh->mFaces = faces;
        mesh->mNumFaces = 12;

        return mesh;
    }

    /// <summary>
    /// Append another bounding volume wire-frame to a scene
    /// (created by construct_mesh(...))
    /// </summary>
    /// <param name="scene"></param>
    /// <param name="vp"></param>
    /// <param name="name"></param>
    void append_node( aiScene* model, boundingVolume_t const* vp, std::string const& name )
    {
        if ( auto box = dynamic_cast<boundingVolumeBox_t const*>( vp ) )
        {
            auto mesh = bounding_box_mesh( box );
            if ( !name.empty() )
                mesh->mName = name.c_str();

            model->mMeshes = (aiMesh**)realloc( model->mMeshes, ( model->mNumMeshes + 1 ) * sizeof( *model->mMeshes ) );
            model->mMeshes[model->mNumMeshes++] = mesh;

            aiNode* node_1 = new aiNode( "node[" + std::to_string( model->mNumMeshes ) + "]" );
            node_1->mMeshes = new unsigned int{ model->mNumMeshes - 1 };
            node_1->mNumMeshes = 1;

            aiNode* root = model->mRootNode;
            node_1->mParent = root;

            root->mChildren = (aiNode**)realloc( root->mChildren, ( root->mNumChildren + 1 ) * sizeof( *root->mChildren ) );
            root->mChildren[root->mNumChildren++] = node_1;
        }
        else
        {
            throw std::runtime_error( "Is not supported" );
        }
    }

    void append_node( aiScene* model, boundingVolume_t const* vp )
    {
        return append_node( model, vp, "" );
    }

    std::string make_mesh_name( child_t const& ch )
    {
        if ( !ch.content.uri.empty() )
            return ch.content.uri;
        return std::to_string( ch.geometricError );
    }

    bool save_as( const char* filename ) const
    {
        auto cp = strdup( filename );
        int flags = 0
            | JSON_C_TO_STRING_NOZERO
            | JSON_C_TO_STRING_PRETTY
            ;
        int rc = json_object_to_file_ext( cp, m_root, flags );
        free( cp );
        return rc == 0;
    }

    bool save_as( std::string const& s ) const
    {
        return save_as( s.c_str() );
    }

    /// @brief Build a scene containing all bounding boxes of the original scene
    /// @return 
    std::unique_ptr<aiScene> build_bounding_boxes_scene()
    {
        auto root_tile = this->root();
        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ this->construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender](
            TilesetJson& actual,
            TilesetJson::child_t const& root_tile,
            aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get(), actual.make_mesh_name( *cp ) );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( *this, root_tile, scene.get() );

        return scene;
    };

protected:
    static int get_as_int( json_object* o, const char* key, int dflt = 0 )
    {
        json_object* ko;
        if ( json_object_object_get_ex( o, key, &ko ) )
            return json_object_get_int( ko );
        return dflt;
    }

    static std::string get_as_string( json_object* o, const char* key, const char* dflt = "" )
    {
        json_object* ko;
        if ( json_object_object_get_ex( o, key, &ko ) )
            return json_object_get_string( ko );
        return dflt;
    }

    static double get_as_double( json_object* o, const char* key, double dflt = 0 )
    {
        json_object* ko;
        if ( json_object_object_get_ex( o, key, &ko ) )
            return json_object_get_double( ko );
        return dflt;
    }

    double geometricError( json_object* o )
    {
        return get_as_double( o, "geometricError", -99999.0 );
    }

    /// <summary>
    /// "boundingVolume": { "box": [ 0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5 ] },
    /// </summary>
    /// <param name="o"></param>
    /// <returns></returns>
    boundingVolume_t* boundingVolume( json_object* o )
    {
        json_object* bv = 0;
        if ( json_object_object_get_ex( o, "boundingVolume", &bv ) )
        {
            json_object* box_obj = 0;
            if ( json_object_object_get_ex( bv, "box", &box_obj ) )
            {
                size_t arr_size = json_object_array_length(box_obj);
                auto p = new boundingVolumeBox_t;
                if (arr_size != SIZE(boundingVolumeBox_t::box))
                    throw std::runtime_error( "Wrong boundingVolume.box array size" );

                for ( int i = 0; i < SIZE( boundingVolumeBox_t::box ); ++i )
                {
                    json_object* oi = json_object_array_get_idx( box_obj, i );
                    p->box[i] = json_object_get_double( oi );
                }
                return p;
            }
            else if ( json_object_object_get_ex( bv, "region", &box_obj ) )
            {
                auto box_size = json_object_array_length( box_obj );
                auto p = new boundingVolumeRegion_t;
                if (box_size != SIZE(boundingVolumeRegion_t::region))
                    throw std::runtime_error( "Wrong boundingVolume.region array size" );

                for ( int i = 0; i < SIZE( boundingVolumeRegion_t::region ); ++i )
                {
                    json_object* oi = json_object_array_get_idx( box_obj, i );
                    p->region[i] = json_object_get_double( oi );
                }
                return p;
            }
            throw std::runtime_error( "Unsupported boundingVolume type" );
        }

        return NULL;
    }

    /// <summary>
    /// "transform": [ 0.08482280432458418, -0.019589455922448373,
    ///     -0.04920716495771534, 0, 0.047019079751234455, 0.07061686167023891,
    ///     0.05293831303690587, 0, 0.024378228111583623, -0.06804051781481994,
    ///     0.06910998429771772, 0, 1557374.754207028, -4346689.358346916, 4385455.219789558, 1 ]
    /// </summary>
    /// <param name="o"></param>
    /// <param name="t"></param>
    std::unique_ptr<transform_t> transform( json_object* t_obj )
    {
        auto t = std::make_unique<transform_t>();

        if ( json_object_get_type( t_obj ) != json_type_array )
            throw std::runtime_error( "transform is not array" );
        size_t t_array_size = json_object_array_length( t_obj );
        if (t_array_size != SIZE(t->matrix4x4))
            throw std::runtime_error( "transform array wrong size" );
        for ( int i = 0; i != t_array_size; ++i )
        {
            json_object* oi = json_object_array_get_idx( t_obj, i );
            t->matrix4x4[i] = json_object_get_double( oi );
        }
        return t;
    }

    std::string top_attribute( std::string const& key )
    {
        std::string s;
#if 0
        json_object* root = m_doco.GetRoot();
        if ( root.IsValid() )
            s = root.GetString( key );
#endif
        return s;
    }

    static int json_double_to_string_g18
        (
        struct json_object* jso,
        struct printbuf* pb,
        int /* level */,
        int /* flags */
        )
    {
        char szBuffer[75] = {};
        snprintf( szBuffer, sizeof( szBuffer ), "%.18g", json_object_get_double( jso ) );
        return printbuf_memappend( pb, szBuffer, static_cast<int>( strlen( szBuffer ) ) );
    }

    static json_object *json_object_new_double_g18( double d )
    {
        json_object*jso = json_object_new_double( d );
        json_object_set_serializer( jso, json_double_to_string_g18, nullptr, nullptr );
        return jso;
    }

    json_object* m_root;
};


} // namespace

#endif // TILESETJSON_H
