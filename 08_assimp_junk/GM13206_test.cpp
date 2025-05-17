//
// Created on 2024/04/08
//

#include "SDKTestF.h"

#include "colorprint.h"

#include <assimp/DefaultLogger.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/aabb.h>
#include <assimp/postprocess.h>
#include <assimp/Exceptional.h>
#include <assimp/SceneCombiner.h>
#include <assimp/cimport.h>
#include <assimp/version.h>

#include "geojson/libjson/json.h"

#include <Eigen/Dense>

#include <filesystem>
#include <map>
#include <algorithm>
#include <limits>


#pragma warning (disable : 4390)

namespace fs = std::filesystem;

std::ostream &operator<<( std::ostream& ss, aiVector3D const& v )
{
    ss << "{" << std::setprecision( 9 ) << v.x
        << "," << std::setprecision( 9 ) << v.y
        << "," << std::setprecision( 9 ) << v.z
        << "}";
    return ss;
};

std::ostream& operator<<( std::ostream& ss, aiAABB const& bb )
{
    ss << "{" << bb.mMin << "," << bb.mMax << "}";
    return ss;
};


#include "../../../FormatSupport/TilesetJson.h"

using namespace cesiumjs;

namespace meshtoolbox
{
#if 0
}
#endif

typedef aiVector3t<double> Vector3D;

/// @brief (origin, dimensions)
typedef std::pair<Vector3D, Vector3D> box_t;

struct cylinder_t
{
    Vector3D origin;
    ai_real height;
    ai_real radius;
};

struct lidarptr_t
{
    ai_real x, y, z;
    uint16_t r, g, b;
};

#ifndef M_PI
/** PI definition */
# define M_PI 3.14159265358979323846
/* 3.1415926535897932384626433832795 */
#endif

/// @brief Helper class to produce meshes and scenes
class Toolbox
{
public:
    /// @brief Create a box mesh
    /// @param origin a corner
    /// @param dimensions another corner
    /// @param name the mesh name. can be null
    /// @return an object pointer. The user owns the object
    aiMesh* mesh_box( aiVector3D const& origin, aiVector3D const& dimensions, const char* name = nullptr)
    {
        auto box = new aiMesh{};
        box->mName = name ? name : "";
        box->mNumVertices = 8;
        box->mVertices = new aiVector3D[box->mNumVertices];

        auto const &v = dimensions;
        box->mVertices[0] = origin;
        box->mVertices[1] = aiVector3D{ origin.x + v.x, origin.y, origin.z };
        box->mVertices[2] = aiVector3D{ origin.x + v.x, origin.y + v.y, origin.z };
        box->mVertices[3] = aiVector3D{ origin.x, origin.y + v.y, origin.z };

        box->mVertices[4] = aiVector3D{ origin.x, origin.y, origin.z + v.z };
        box->mVertices[5] = aiVector3D{ origin.x + v.x, origin.y, origin.z + v.z };
        box->mVertices[6] = origin + dimensions;
        box->mVertices[7] = aiVector3D{ origin.x, origin.y + v.y, origin.z + v.z };

        box->mNumFaces = 12;
        box->mFaces = new aiFace[box->mNumFaces];
        box->mFaces[0].mNumIndices = 3;
        box->mFaces[0].mIndices = new unsigned[3] { 0, 2, 1 };

        box->mFaces[1].mNumIndices = 3;
        box->mFaces[1].mIndices = new unsigned[3] { 3, 2, 0 };

        box->mFaces[2].mNumIndices = 3;
        box->mFaces[2].mIndices = new unsigned[3] { 4, 5, 6 };

        box->mFaces[3].mNumIndices = 3;
        box->mFaces[3].mIndices = new unsigned[3] { 4, 6, 7 };

        box->mFaces[4].mNumIndices = 3;
        box->mFaces[4].mIndices = new unsigned[3] { 0, 4, 3 };

        box->mFaces[5].mNumIndices = 3;
        box->mFaces[5].mIndices = new unsigned[3] { 3, 4, 7 };

        box->mFaces[6].mNumIndices = 3;
        box->mFaces[6].mIndices = new unsigned[3] { 7, 2, 3 };

        box->mFaces[7].mNumIndices = 3;
        box->mFaces[7].mIndices = new unsigned[3] { 6, 2, 7 };

        box->mFaces[8].mNumIndices = 3;
        box->mFaces[8].mIndices = new unsigned[3] { 6, 1, 2 };

        box->mFaces[9].mNumIndices = 3;
        box->mFaces[9].mIndices = new unsigned[3] { 5, 1, 6 };

        box->mFaces[10].mNumIndices = 3;
        box->mFaces[10].mIndices = new unsigned[3] { 0, 1, 4 };

        box->mFaces[11].mNumIndices = 3;
        box->mFaces[11].mIndices = new unsigned[3] { 4, 1, 5 };

        box->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;

        return box;
    };


    /// @brief Create a box point cloud "mesh"
    /// @param origin a corner
    /// @param dimensions another corner
    /// @param name the mesh name. can be null
    /// @return an object pointer. The user owns the object
    aiMesh* mesh_pc_box( aiVector3D const& origin, aiVector3D const& dimensions, const char* name = nullptr )
    {
        auto box = new aiMesh{};
        box->mName = name ? name : "";
        box->mNumVertices = 8;
        box->mVertices = new aiVector3D[box->mNumVertices];

        auto const& v = dimensions;
        box->mVertices[0] = origin;
        box->mVertices[1] = aiVector3D{ origin.x + v.x, origin.y, origin.z };
        box->mVertices[2] = aiVector3D{ origin.x + v.x, origin.y + v.y, origin.z };
        box->mVertices[3] = aiVector3D{ origin.x, origin.y + v.y, origin.z };

        box->mVertices[4] = aiVector3D{ origin.x, origin.y, origin.z + v.z };
        box->mVertices[5] = aiVector3D{ origin.x + v.x, origin.y, origin.z + v.z };
        box->mVertices[6] = origin + dimensions;
        box->mVertices[7] = aiVector3D{ origin.x, origin.y + v.y, origin.z + v.z };

        box->mNumFaces = 8;
        box->mFaces = new aiFace[box->mNumFaces];

        box->mFaces[0].mNumIndices = 1;
        box->mFaces[0].mIndices = new unsigned[1] { 0 };

        box->mFaces[1].mNumIndices = 1;
        box->mFaces[1].mIndices = new unsigned[1] { 1 };

        box->mFaces[2].mNumIndices = 1;
        box->mFaces[2].mIndices = new unsigned[1] { 2 };

        box->mFaces[3].mNumIndices = 1;
        box->mFaces[3].mIndices = new unsigned[1] { 3 };

        box->mFaces[4].mNumIndices = 1;
        box->mFaces[4].mIndices = new unsigned[1] { 4 };

        box->mFaces[5].mNumIndices = 1;
        box->mFaces[5].mIndices = new unsigned[1] { 5 };

        box->mFaces[6].mNumIndices = 1;
        box->mFaces[6].mIndices = new unsigned[1] { 6 };

        box->mFaces[7].mNumIndices = 1;
        box->mFaces[7].mIndices = new unsigned[1] { 7 };

        box->mPrimitiveTypes = aiPrimitiveType_POINT;

        return box;
    };

    /// @brief Create a "mesh" of Lidar points
    /// @param data
    /// @param count
    /// @param name [optional]
    /// @return 
    aiMesh* mesh_lidar_pc( lidarptr_t* data, size_t count, const char* name = nullptr )
    {
        auto box = new aiMesh{};
        box->mName = name ? name : "";
        box->mNumVertices = unsigned(count);
        box->mVertices = new aiVector3D[box->mNumVertices];

        box->mColors[0] = new aiColor4D [box->mNumVertices];

        ai_real d = ai_real( 1.0 / 0xFFFF );
        for ( unsigned i = 0; i < box->mNumVertices; ++i )
        {
            box->mVertices[i] = aiVector3D{ data[i].x, data[i].z, - data[i].y};
            box->mColors[0][i] = aiColor4D{ data[i].r * d, data[i].g * d, data[i].b * d, 1.0 };
        }

        box->mNumFaces = box->mNumVertices;
        box->mFaces = new aiFace[box->mNumFaces];

        for ( unsigned i = 0; i < box->mNumFaces; ++i )
        {
            box->mFaces[i].mNumIndices = 1;
            box->mFaces[i].mIndices = new unsigned[1] { i };
        }

        box->mPrimitiveTypes = aiPrimitiveType_POINT;

        return box;
    };

    /// @brief Create a single box scene
    /// @param origin coordinate of the origin (a corner)
    /// @param dimensions (x,y,z) dimensions of the box
    /// @return a pointer to a single mesh (the box) scene. The user owns the object.
    aiScene* make_box( aiVector3D const& origin, aiVector3D const& dimensions )
    {
        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();

        auto box = mesh_box( origin, dimensions, "box");

        model->mNumMeshes = 1;
        model->mMeshes = new aiMesh * [model->mNumMeshes] { box };

        model->mRootNode = new aiNode( "ROOT" );
        model->mRootNode->mNumChildren = 1;
        model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
        model->mRootNode->mChildren[0] = new aiNode( "box" );
        model->mRootNode->mChildren[0]->mNumMeshes = 1;
        model->mRootNode->mChildren[0]->mMeshes = new unsigned[1] { 0 };

        return model.release();
    }

    /// @brief Create a multi-object scene
    /// @param boxes object coordinates
    /// @return 
    aiScene* make_boxes( std::vector<box_t> const &boxes )
    {
        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();
        model->mRootNode = new aiNode( "ROOT" );
        if ( boxes.size() )
        {
            model->mRootNode->mNumChildren = unsigned(boxes.size());
            model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
            model->mNumMeshes = unsigned(boxes.size());
            model->mMeshes = new aiMesh * [model->mNumMeshes];

            for ( unsigned int i = 0; i != boxes.size(); ++i )
            {
                std::ostringstream box_name; box_name << "box " << boxes[i].second;
                model->mMeshes[i] = mesh_box(boxes[i].first, boxes[i].second, box_name.str().c_str() );
                model->mRootNode->mChildren[i] = new aiNode( "box_" + std::to_string(i) );
                model->mRootNode->mChildren[i]->mNumMeshes = 1;
                model->mRootNode->mChildren[i]->mMeshes = new unsigned[1] { i };
            }
        }
        return model.release();
    }

    /// @brief Create a multi-object scene
    /// @param boxes object coordinates
    /// @return 
    aiScene* make_boxes_with_transform( std::vector<box_t> const& boxes )
    {
        static_assert( sizeof( ai_real ) == 4 );

        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();
        model->mRootNode = new aiNode( "ROOT" );
        if ( boxes.size() )
        {
            Vector3D origin{ 0,0,0 };
            for ( auto& bp : boxes )
            {
                origin += bp.first;
            }

            // compute the origin as a centroid of the boxes
            origin /= static_cast<ai_real>( boxes.size() );

            model->mRootNode->mNumChildren = unsigned( boxes.size() );
            model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
            model->mRootNode->mTransformation = make_offset_transformation(origin);
            model->mNumMeshes = unsigned( boxes.size() );
            model->mMeshes = new aiMesh * [model->mNumMeshes];

            for ( unsigned int i = 0; i != boxes.size(); ++i )
            {
                std::ostringstream box_name; box_name << "box " << boxes[i].second;
                aiVector3D origin_i{ boxes[i].first };
                origin_i -= origin;
                model->mMeshes[i] = mesh_box( origin_i, boxes[i].second, box_name.str().c_str() );
                model->mRootNode->mChildren[i] = new aiNode( "box_" + std::to_string( i ) );
                model->mRootNode->mChildren[i]->mNumMeshes = 1;
                model->mRootNode->mChildren[i]->mMeshes = new unsigned[1] { i };
            }
        }
        return model.release();
    }

    aiMatrix4x4 make_offset_transformation( Vector3D const& origin )
    {
        aiMatrix4x4 r;
        aiMatrix4x4::Translation( origin, r );
        return r;
    }

    /// @brief 
    /// @param v0 
    /// @return 
    Eigen::Vector3f to_eigen( aiVector3D const& v0 )
    {
        return Eigen::Vector3f( v0.x, v0.y, v0.z );
    }

    aiVector3D from_eigen( Eigen::Vector3f const&v )
    {
        return { v.x(), v.y(), v.z() };
    }

    aiVector3D plane_normal( aiVector3D const& v0, aiVector3D const& v1, aiVector3D const& v2 )
    {
        Eigen::Vector3f normal = to_eigen( v1 - v0 ).cross( to_eigen( v2 - v0 ) );
        return { normal.x(), normal.y(), normal.z() };
    }

    aiMesh* mesh_cylinder( aiVector3D const& origin, ai_real H, ai_real R, unsigned N, const char *name )
    {
        auto box = new aiMesh{};
        box->mName = name ? name : "";
        box->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;

        // 
        box->mNumVertices = ( N + 1 ) * 2;
        box->mVertices = new aiVector3D[box->mNumVertices];

        ai_real a = 0, da = ai_real(2 * M_PI / N);
        for ( unsigned i = 0; i < N; ++i )
        {
            ai_real x = R * cos( a + i * da ) + origin.x;
            ai_real z = R * sin( a + i * da ) + origin.z;

            box->mVertices[2*i] = { x, 0, z };
            box->mVertices[2*i + 1] = { x, H, z };
        }
        box->mVertices[N * 2] = origin;
        box->mVertices[N * 2 + 1] = { origin.x, origin.y + H, origin.z };

        box->mNumFaces = 4 * N;
        box->mFaces = new aiFace[box->mNumFaces];

        // the cylinder sides
        for ( unsigned i = 0; i < N; ++i )
        {
            auto i_next = (i + 1) % N;

            box->mFaces[2 * i].mNumIndices = 3;
            box->mFaces[2 * i].mIndices = new unsigned[3] { 2 * i, 2 * i + 1, 2 * i_next };

            box->mFaces[2 * i + 1].mNumIndices = 3;
            box->mFaces[2 * i + 1].mIndices = new unsigned[3] { 2 * i_next, 2 * i + 1, 2 * i_next + 1 };
        }

        {
            // The cylinder bottom
            int pos = 2 * N;
            for ( unsigned i = 0; i < N; ++i )
            {
                auto i_next = ( i + 1 ) % N;
                box->mFaces[pos].mNumIndices = 3;
                box->mFaces[pos].mIndices = new unsigned[3] { 2 * N, 2 * i, 2 * i_next };
                pos += 1;
            }
            // The cylinder top
            for ( unsigned i = 0; i < N; ++i )
            {
                auto i_next = ( i + 1 ) % N;
                box->mFaces[pos].mNumIndices = 3;
                box->mFaces[pos].mIndices = new unsigned[3] { 2 * i + 1, 2 * N + 1, 2 * i_next + 1};
                pos += 1;
            }
        }

        return box;
    }

    /// @brief Create a single cylinder scene
    /// @param origin coordinate of the origin (a center of the cylinder bottom)
    aiScene* make_cylinder( aiVector3D const& origin, ai_real H, ai_real R, unsigned N )
    {
        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();

        auto box = mesh_cylinder( origin, H, R, N, "cylinder" );

        model->mNumMeshes = 1;
        model->mMeshes = new aiMesh * [model->mNumMeshes] { box };

        model->mRootNode = new aiNode( "ROOT" );
        model->mRootNode->mNumChildren = 1;
        model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
        model->mRootNode->mChildren[0] = new aiNode( "cylinder" );
        model->mRootNode->mChildren[0]->mNumMeshes = 1;
        model->mRootNode->mChildren[0]->mMeshes = new unsigned[1] { 0 };

        return model.release();
    }

    /// @brief Create a multi-object scene
    /// @param cylinders object parameters
    /// @param N
    /// @return 
    aiScene* make_cylinders( std::vector<cylinder_t> const& cylinders, unsigned N )
    {
        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();
        model->mRootNode = new aiNode( "ROOT" );
        if ( cylinders.size() )
        {
            model->mRootNode->mNumChildren = unsigned( cylinders.size() );
            model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
            model->mNumMeshes = unsigned( cylinders.size() );
            model->mMeshes = new aiMesh * [model->mNumMeshes];

            for ( unsigned int i = 0; i != cylinders.size(); ++i )
            {
                std::ostringstream box_name; box_name << "cylinder_" << i;
                model->mMeshes[i] = mesh_cylinder(
                    cylinders[i].origin, cylinders[i].height, cylinders[i].radius, N,
                    box_name.str().c_str() );
                model->mRootNode->mChildren[i] = new aiNode( "cyl_" + std::to_string( i ) );
                model->mRootNode->mChildren[i]->mNumMeshes = 1;
                model->mRootNode->mChildren[i]->mMeshes = new unsigned[1] { i };
            }
        }
        return model.release();
    }

    /// @brief Create a single PC-box scene
    /// @param origin coordinate of the origin (a corner)
    /// @param dimensions (x,y,z) dimensions of the box
    /// @return a pointer to a single mesh (the box) scene. The user owns the object.
    aiScene* make_pc_box( aiVector3D const& origin, aiVector3D const& dimensions )
    {
        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();

        auto box = mesh_pc_box( origin, dimensions, "box" );

        model->mNumMeshes = 1;
        model->mMeshes = new aiMesh * [model->mNumMeshes] { box };

        model->mRootNode = new aiNode( "ROOT" );
        model->mRootNode->mNumChildren = 1;
        model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
        model->mRootNode->mChildren[0] = new aiNode( "pc-box" );
        model->mRootNode->mChildren[0]->mNumMeshes = 1;
        model->mRootNode->mChildren[0]->mMeshes = new unsigned[1] { 0 };

        return model.release();
    }

    aiScene* make_lidar_pc( lidarptr_t* data, size_t count )
    {
        std::unique_ptr<aiScene> model = std::make_unique<aiScene>();

        auto box = mesh_lidar_pc( data, count, "box" );

        model->mNumMeshes = 1;
        model->mMeshes = new aiMesh * [model->mNumMeshes] { box };

        model->mRootNode = new aiNode( "ROOT" );
        model->mRootNode->mNumChildren = 1;
        model->mRootNode->mChildren = new aiNode * [model->mRootNode->mNumChildren];
        model->mRootNode->mChildren[0] = new aiNode( "pc-box" );
        model->mRootNode->mChildren[0]->mNumMeshes = 1;
        model->mRootNode->mChildren[0]->mMeshes = new unsigned[1] { 0 };

        return model.release();
    }
};

};

class GM13206_test : public GlobalMapperSDKTestF
{
protected:
    fs::path create_workspace()
    {
        auto ws = fs::path( "out" ) / test_case_name() / test_name();
        if ( fs::is_directory( ws ) )
            fs::remove_all( ws );
        std::error_code err;
        fs::create_directories( ws, err );

        CONSOLE( "ws = " << fs::absolute(ws).string() );

        return ws;
    }

};

TEST_F( GM13206_test, meshtoolbox_0 )
{
    meshtoolbox::Toolbox tb;

    {
        auto actual = std::unique_ptr<aiMesh>(tb.mesh_box( { 0,0,0 }, { 1,1,1}, "TheBox"));
        ASSERT_TRUE( actual );

        EXPECT_STREQ( "TheBox", actual->mName.C_Str() );
        EXPECT_TRUE( actual->HasFaces() );
        EXPECT_TRUE( actual->HasPositions() );
        EXPECT_FALSE( actual->HasNormals() );
        EXPECT_FALSE( actual->HasBones() );

        EXPECT_EQ( 8, actual->mNumVertices );
        EXPECT_EQ( 12, actual->mNumFaces );

        auto get_face_positions = []( aiMesh const* mp, unsigned k )->std::tuple<aiVector3D, aiVector3D, aiVector3D>
            {
                if ( !( k >= 0 && k < mp->mNumFaces ) )
                    throw std::runtime_error( "Failed: k >= 0 && k < mp->mNumFaces" );
                if ( !( mp->mFaces[k].mNumIndices == 3 ) )
                    throw std::runtime_error( "Failed: mp->mFaces[k].mNumIndices == 3" );
                return { mp->mVertices[mp->mFaces[k].mIndices[0]], mp->mVertices[mp->mFaces[k].mIndices[1]],mp->mVertices[mp->mFaces[k].mIndices[2]] };
            };

        /* Compute normals to the faces */
        for ( unsigned i = 0; i < actual->mNumFaces; ++i )
        {
            auto [v0, v1, v2] = get_face_positions( actual.get(), i );

            auto normal = tb.from_eigen( tb.to_eigen( v1 - v0 ).cross( tb.to_eigen( v2 - v0 ) ) );
            CONSOLE_T( i << ": " << "v0: " << v0 << ", v1: " << v1 << ", v2: " << v2
                     << ", normal: " << normal << ", plane: " << tb.plane_normal( v0, v1, v2 ) );
        }
    }
    {
        auto actual = std::unique_ptr<aiMesh>( tb.mesh_box( { 0,0,0 }, { 1,1,1 } ) );
        ASSERT_TRUE( actual );
        EXPECT_STREQ( "", actual->mName.C_Str() );
        EXPECT_EQ( 8, actual->mNumVertices );
    }
}


/// @brief 
/// @param  --gtest_filter=GM13206_test.t_0
/// @param  
TEST_F( GM13206_test, t_0 )
{
    CONSOLE_EVAL( aiGetLegalString() );
}

/// @brief Dump a glTF model as XML file
/// @param --gtest_filter=GM13206_test.BoxTextured_assxml
/// @param  
TEST_F( GM13206_test, BoxTextured_assxml )
{
    // OneDrive=C:\Users\yurikuznetsov\OneDrive
    auto filename = R"(C:\Users\yurikuznetsov\OneDrive\V\Asset Importer Test Models\models\glTF2\BoxTextured-glTF\BoxTextured.gltf)";
    if ( !fs::is_regular_file( filename ) )
        GTEST_SKIP() << filename;

    auto ws = create_workspace();
    Assimp::Importer impo;
    auto model = impo.ReadFile( filename, aiProcess_ValidateDataStructure );
    ASSERT_TRUE( model != nullptr );
    Assimp::Exporter expo;
    ASSERT_EQ( 0, expo.Export( model, "assxml", ( ws / "BoxTextured.xml" ).string() ) );
    CONSOLE( "Saved as: " << fs::absolute( ws / "BoxTextured.xml" ).string() );
}

/// @brief Compose a single box scene and save it as a .GLB file
/// @param --gtest_filter=GM13206_test.tb_0
TEST_F( GM13206_test, tb_0 )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;

    auto actual = std::unique_ptr<aiScene>( tb.make_box( { 0, 0, 0 }, { 100, 50, 10 } ) );
    ASSERT_TRUE( (bool)actual );

    {
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
        exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
    }
}

/// @brief Compose a single cylinder scene and save it as a .GLB file
/// @param --gtest_filter=GM13206_test.tc_0
TEST_F( GM13206_test, tc_0 )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;

    {
        auto actual = std::unique_ptr<aiScene>( tb.make_cylinder( { 1, 0, 2 }, 2.0, 3.0, 6 ) );
        ASSERT_TRUE( actual );
    
        {
            Assimp::Exporter exp;
            exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
            exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
        }
    }
    {
        auto actual = std::unique_ptr<aiScene>( tb.make_cylinder( { 0, 0, 0 }, 5.0, 1.0, 40 ) );
        ASSERT_TRUE( actual );

        {
            Assimp::Exporter exp;
            exp.Export( actual.get(), "glb2", ( ws / "actual2.glb" ).string() );
            exp.Export( actual.get(), "assxml", ( ws / "actual2.xml" ).string() );
        }
    }
    {
        auto actual = std::unique_ptr<aiScene>( tb.make_cylinder( { 0, 0, 0 }, 10.0, 1.0, 400 ) );
        ASSERT_TRUE( actual );

        {
            Assimp::Exporter exp;
            exp.Export( actual.get(), "glb2", ( ws / "actual3.glb" ).string() );
            exp.Export( actual.get(), "assxml", ( ws / "actual3.xml" ).string() );
        }
    }
}

/// @brief Build a multi-object scene
/// @param --gtest_filter=GM13206_test.tb_1
TEST_F( GM13206_test, tb_1 )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;

    auto boxes = std::vector<meshtoolbox::box_t>{
        { { 0, 0, 0 }, { 100, 10, 50 } },
        { {500, 0, 200}, { 50, 80, 30 } },
    };

    auto actual = std::unique_ptr<aiScene>( tb.make_boxes( boxes ) );
    ASSERT_TRUE( (bool)actual );

    {
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
        exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
    }
    {
        // aiProcess_GenNormals
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual2.glb" ).string(), aiProcess_GenNormals );
        exp.Export( actual.get(), "assxml", ( ws / "actual2.xml" ).string(), aiProcess_GenNormals );
    }
}

/// @brief Build a multi-object scene. Use non-zero offsets
/// @param --gtest_filter=GM13206_test.tb_1_offset
TEST_F( GM13206_test, tb_1_offset )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;

    auto boxes = std::vector<meshtoolbox::box_t>{
        { { 0, 0, 0 }, { 100, 10, 50 } },
        { {500, 0, 200}, { 50, 80, 30 } },
    };

    auto actual = std::unique_ptr<aiScene>( tb.make_boxes_with_transform( boxes ) );
    ASSERT_TRUE( (bool)actual );

    {
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
        exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
    }
    {
        // aiProcess_GenNormals
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual2.glb" ).string(), aiProcess_GenNormals );
        exp.Export( actual.get(), "assxml", ( ws / "actual2.xml" ).string(), aiProcess_GenNormals );
    }
}

/// @brief Build a multi-object scene
///
/// Generate a lot of small boxes scattered over a large terrain
/// @param --gtest_filter=GM13206_test.tb_2
TEST_F( GM13206_test, tb_2 )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;


    // A "graveyard" 
    ai_real max_X = 2000;
    ai_real max_Z = 1000;
    ai_real step = 50; // 50m interval between "obelisks" 

    std::vector<meshtoolbox::box_t> boxes;
    boxes.reserve( size_t( max_X / step + 1 ) * size_t( max_Z / step + 1 ) );
    for ( int i = 0; i < max_X / step; ++i )
    {
        for ( int j = 0; j < max_Z / step; ++j )
        {
            boxes.emplace_back( meshtoolbox::box_t{ {i * step, 0, j * step}, { step / 5, step, 4} } );
        }
    }

    auto actual = std::unique_ptr<aiScene>( tb.make_boxes( boxes ) );
    ASSERT_TRUE( (bool)actual );

    {
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
        //exp.Export( actual.get(), "obj", ( ws / "actual.obj" ).string() );
        //exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
    }
}

/// @brief Build a multi-object scene
///
/// Generate a lot of small boxes scattered over a large terrain
/// with huge boxes at the corners
/// @param --gtest_filter=GM13206_test.tb_3
TEST_F( GM13206_test, tb_3 )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;


    // A "graveyard" 
    ai_real max_X = 2000;
    ai_real max_Z = 1000;
    ai_real step = 50; // 50m interval between "obelisks" 

    std::vector<meshtoolbox::box_t> boxes;
    boxes.reserve( size_t( max_X / step + 1 ) * size_t( max_Z / step + 1 ) );
    for ( int i = 0; i < max_X / step; ++i )
    {
        for ( int j = 0; j < max_Z / step; ++j )
        {
            boxes.emplace_back( meshtoolbox::box_t{ {i * step, 0, j * step}, { step / 5, step, 4} } );
        }
    }

    boxes.emplace_back( meshtoolbox::box_t{ {max_X, 0, max_Z}, { step, step * 10, step} } );
    boxes.emplace_back( meshtoolbox::box_t{ {max_X, 0, 0}, { step, step * 10, step} } );
    boxes.emplace_back( meshtoolbox::box_t{ {0, 0, max_Z}, { step, step * 10, step} } );

    auto actual = std::unique_ptr<aiScene>( tb.make_boxes( boxes ) );
    ASSERT_TRUE( (bool)actual );

    {
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
        //exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
    }
}

/// @brief Build a multi-object scene
///
/// Generate a lot of small cylinders scattered over a large terrain
/// with huge cylinders at the corners
/// @param --gtest_filter=GM13206_test.tc_3
TEST_F( GM13206_test, tc_3 )
{
    auto ws = create_workspace();

    meshtoolbox::Toolbox tb;

    // A "graveyard" 
    ai_real max_X = 2000;
    ai_real max_Z = 1000;
    ai_real step = 50; // 50m interval between objects

    std::vector<meshtoolbox::cylinder_t> params;
    params.reserve( size_t( max_X / step + 1 ) * size_t( max_Z / step + 1 ) );
    for ( int i = 0; i < max_X / step; ++i )
    {
        for ( int j = 0; j < max_Z / step; ++j )
        {
            params.emplace_back(
                meshtoolbox::cylinder_t{
                    {i * step, 0, j * step}, // origin
                    step / 3, // height
                    step / 4  // radius
                } );
        }
    }

    params.emplace_back( meshtoolbox::cylinder_t{ {max_X, 0, max_Z}, step, step} );
    params.emplace_back( meshtoolbox::cylinder_t{ {max_X, 0, 0}, step * 10, step} );
    params.emplace_back( meshtoolbox::cylinder_t{ {0, 0, max_Z}, step * 10, step} );

    auto actual = std::unique_ptr<aiScene>( tb.make_cylinders( params, 40 ) );
    ASSERT_TRUE( (bool)actual );

    {
        Assimp::Exporter exp;
        exp.Export( actual.get(), "glb2", ( ws / "actual.glb" ).string() );
        //exp.Export( actual.get(), "assxml", ( ws / "actual.xml" ).string() );
    }
}

#if 0
#   include "../FormatSupport/Cesium3dTiles_Transforms.h"
#else
namespace cesiumjs
{
#if 0
}
#endif

struct Cartesian
{
    double x, y, z;
};

struct Cartographic
{
    double longitude, latitude, height;

    static Cartographic fromDegrees( double lon, double lat, double elev )
    {
        return { lon / 180 * M_PI, lat / 180 * M_PI, elev };
    };

    static Cartesian toCartesian( Cartographic  const& c )
    {
        return c.toCartesian();
    }

    Cartesian toCartesian() const
    {
        GM_Projection_t geo_proj;
        if ( GM_LoadProjectionFromEPSGCode( 4326, &geo_proj ) )
            throw std::runtime_error( "Cannot load EPSG:4326 projection" );
        geo_proj.mUnit = GM_PRJ_UNIT_RADIANS;

        double x, y, z;
        if (GM_ProjectPointToECEF(longitude, latitude, height, &x, &y, &z, GM_DATUM_WGS_84 , &geo_proj) )
            throw std::runtime_error( "Cannot project to ECEF" );
        return {x, y, z};
    }
};

struct Matrix4
{
    double d[16];
};

class Transforms
{
public:
    /// @brief Based on the js original implementation
    ///
    /// See also:
    /// https://github.com/CesiumGS/cesium/blob/main/packages/engine/Source/Core/Transforms.js
    /// @param cartesian 
    /// @return 
    static Matrix4 eastNorthUpToFixedFrame( Cartesian const& cartesian, Cartographic  const& cartographic )
    {
        auto const& lambda = cartographic.longitude;
        auto const &phi = cartographic.latitude;
        Matrix4 res = {
                         - sin(lambda),            cos(lambda),           0, 0,
            - cos( lambda ) * sin(phi), - sin(lambda)*sin(phi),    cos(phi), 0,
                cos(lambda) * cos(phi), sin(lambda) * cos(phi),    sin(phi), 0,
                           cartesian.x,            cartesian.y, cartesian.z, 1
        };

        return res;
    }
};

std::ostream& operator<<( std::ostream& ss, Cartesian const& v )
{
    ss << "{" << std::setprecision( 18 ) << v.x
        << "," << std::setprecision( 18 ) << v.y
        << "," << std::setprecision( 18 ) << v.z
        << "}";
    return ss;
};

std::ostream& operator<<( std::ostream& ss, Matrix4 const& m )
{
    ss << "{\n";
    auto* p = m.d;
    for ( int i = 0; i != 4; ++i )
    {
        for ( int j = 0; j != 4; ++j )
            ss << " " << std::setprecision( 18 ) << *p++;
        ss << "\n";
    }
    ss << "}";
    return ss;
}

}; // namespace cesiumjs
#endif

/// @brief Build a cesium tile transform
///
/// Based on Javascript code from https://github.com/CesiumGS/3d-tiles-tools
/// @param --gtest_filter=GM13206_test.tr_0
/// @param  
TEST_F( GM13206_test, tr_0 )
{
#pragma region transform_json
    auto transform_json = R"(
{
    "geo_location_degrees": [ -105.0375565181, 39.937157238, 1624.6099999999  ],
    "cartographic": {
    "longitude": -1.8332511994904759,
    "latitude": 0.6970348876897846,
    "height": 1624.6099999999
    },
    "cartesian": { "x": -1270909.5708341626, "y": -4730693.317416712, "z": 4073680.810265148 },
    "enuMatrix": {
    "0": 0.965755966816189, "1": -0.25945213924523314,    "2": 0,                   "3": 0,
    "4": 0.1665546943044161, "5": 0.6199647853884153,     "6": 0.7667484585465162,  "7": 0,
    "8": -0.19893452783287865, "9": -0.7404918988884135, "10": 0.6419476624433968, "11": 0,
    "12": -1270909.5708341626, "13": -4730693.317416712, "14": 4073680.810265148,  "15": 1
    },
    "transform": [
    0.965755966816189,
    -0.25945213924523314,
    0,
    0,
    0.1665546943044161,
    0.6199647853884153,
    0.7667484585465162,
    0,
    -0.19893452783287865,
    -0.7404918988884135,
    0.6419476624433968,
    0,
    -1270909.5708341626,
    -4730693.317416712,
    4073680.810265148,
    1
    ]
})";
#pragma endregion
    //
    // a javascript sample
    // Mike_Niwot_mesh_NORMALIZE_3D_tiles:
    // -105.0375565181 39.937157238 1624.6099999999
    //
#if 0
    import { Cartographic } from "cesium";
    import { Matrix4 } from "cesium";
    import { Transforms } from "cesium";

    const lonDegrees = -105.0375565181;
    const latDegrees = 39.937157238;
    const height = 1624.6099999999;
    const cartographic = Cartographic.fromDegrees(
        lonDegrees,
        latDegrees,
        height
    );

    const cartesian = Cartographic.toCartesian( cartographic );
    const enuMatrix = Transforms.eastNorthUpToFixedFrame( cartesian );
    const transform = Matrix4.toArray( enuMatrix );
#endif

    using namespace cesiumjs;

    const auto lonDegrees = -105.0375565181;
    const auto latDegrees = 39.937157238;
    const auto height = 1624.6099999999;
    auto cartographic = Cartographic::fromDegrees( lonDegrees, latDegrees, height );
    EXPECT_DOUBLE_EQ( -1.8332511994904759, cartographic.longitude );
    EXPECT_DOUBLE_EQ( 0.6970348876897846, cartographic.latitude );
    EXPECT_DOUBLE_EQ( 1624.6099999999, cartographic.height );

    auto cartesian = Cartographic::toCartesian( cartographic );
    // Expected:
    // "cartesian": { "x": -1270909.5708341626, "y": -4730693.317416712, "z": 4073680.810265148 },
    CONSOLE_TE( cartesian );
    EXPECT_NEAR( -1270909.5708341626, cartesian.x, 1E-8 );
    EXPECT_NEAR( -4730693.317416712, cartesian.y, 1E-8 );
    EXPECT_NEAR( 4073680.810265148, cartesian.z, 1E-8 );

    auto enuMatrix = Transforms::eastNorthUpToFixedFrame( cartesian, cartographic );
    CONSOLE_TE( enuMatrix );
    // Expected:
    // "enuMatrix": {
    //"0": 0.965755966816189, "1": -0.25945213924523314,    "2": 0,                   "3": 0,
    //"4": 0.1665546943044161, "5": 0.6199647853884153,     "6": 0.7667484585465162,  "7": 0,
    //"8": -0.19893452783287865, "9": -0.7404918988884135, "10": 0.6419476624433968, "11": 0,
    //"12": -1270909.5708341626, "13": -4730693.317416712, "14": 4073680.810265148,  "15": 1
    //},
    EXPECT_NEAR( 0.965755966816189, enuMatrix.d[0], 1E-8 );
    EXPECT_NEAR( -0.25945213924523314, enuMatrix.d[1], 1E-8 );
    EXPECT_EQ( 0, enuMatrix.d[2] );
    EXPECT_EQ( 0, enuMatrix.d[3] );
    EXPECT_NEAR( 0.1665546943044161, enuMatrix.d[4], 1E-6 );
    EXPECT_NEAR( 0.6199647853884153, enuMatrix.d[5], 1E-6 );
    EXPECT_NEAR( 0.7667484585465162, enuMatrix.d[6], 1E-6 );
    EXPECT_EQ( 0, enuMatrix.d[7]);
    EXPECT_NEAR( -0.1989345278328786, enuMatrix.d[8], 1E-6 );
    EXPECT_NEAR( -0.7404918988884135, enuMatrix.d[9], 1E-6 );
    EXPECT_NEAR( 0.6419476624433968, enuMatrix.d[10], 1E-6 );
    EXPECT_EQ( 0, enuMatrix.d[11] );
    EXPECT_NEAR( -1270909.5708341626, enuMatrix.d[12], 1E-8 );
    EXPECT_NEAR( -4730693.317416712, enuMatrix.d[13], 1E-8 );
    EXPECT_NEAR( 4073680.810265148, enuMatrix.d[14], 1E-8 );
    EXPECT_EQ( 1.0, enuMatrix.d[15] );
}

/// @brief Compute transform for the 3D Tileset: "GM13206_test.tc_3"
/// @param --gtest_filter=GM13206_test.tr_1
/// @param  
TEST_F( GM13206_test, tr_1 )
{
    using namespace cesiumjs;

    const auto lonDegrees = -69.9382028426;
    const auto latDegrees = 43.8979944507;
    const auto height = -5.491136942;

    auto cartographic = Cartographic::fromDegrees( lonDegrees, latDegrees, height );
    auto cartesian = Cartographic::toCartesian( cartographic );
    auto enuMatrix = Transforms::eastNorthUpToFixedFrame( cartesian, cartographic );
    CONSOLE_TE( enuMatrix );

    {
        double expected_transform[16] = {
            0.9393231837175706,0.34303346269815504,1.3877787807814454e-17,0,
            -0.23785137742361814,0.6513047191834912,0.7205753846940653,0,
            0.24718146934666044,-0.6768531644593428,0.6933766039988649,0,
            1579099.101018939,-4324022.453400844,4399927.829667903,1
        };
        for ( int i = 0; i != 16; ++i )
            EXPECT_NEAR( expected_transform[i], enuMatrix.d[i], 1E-6 ) << i;
    }
}

/* Moved from Runtime tests */

#include <fstream>

/// <summary>
/// See also https://github.com/CesiumGS/3d-tiles/blob/main/specification/TileFormats/Batched3DModel/README.adoc#tileformats-batched3dmodel-batched-3d-model
/// </summary>
class Batched3DModel
{
public:
    Batched3DModel()
    {}

    virtual ~Batched3DModel() {}

    bool read( std::string const& filename );

    bool false_because( std::string s )
    {
        CONSOLE( "Error: " << s );
        return false;
    }

    std::string featureTableJSON;
    std::vector<char> binary_glTF;
public:
    struct header_t
    {
        char magic[4];
        uint32_t version;
        uint32_t byteLength;
        uint32_t featureTableJSONByteLength;
        uint32_t featureTableBinaryByteLength;
        uint32_t batchTableJSONByteLength;
        uint32_t batchTableBinaryByteLength;
    };
};

bool
Batched3DModel::read( std::string const& filename )
{
    static_assert( 28 == sizeof( header_t ), "Batched3DModel::header_t must be 28 bytes long" );

    CONSOLE_EVAL( filename );

    std::ifstream b3dm( filename, std::ios::binary );

    if ( b3dm.bad() )
        return false_because( "Cannot open " + filename );

    header_t hdr{};
    b3dm.read( (char*)&hdr, sizeof( hdr ) );
    if ( !( hdr.magic[0] == 'b' && hdr.magic[1] == '3' && hdr.magic[2] == 'd' && hdr.magic[3] == 'm' ) )
        return false_because( "Wrong " + filename );
    CONSOLE_EVAL( hdr.byteLength );
    CONSOLE_EVAL( hdr.version );
    CONSOLE_EVAL( hdr.featureTableJSONByteLength );
    CONSOLE_EVAL( hdr.featureTableBinaryByteLength );
    CONSOLE_EVAL( hdr.batchTableJSONByteLength );
    CONSOLE_EVAL( hdr.batchTableBinaryByteLength );

    if ( hdr.version != 1 )
        return false_because( "Only Batched 3D Model version 1 is supported.  Version "
                              + std::to_string( hdr.version )
                              + " is not." );

    // Read the JSON
    if ( auto jsob_buffer_size = hdr.featureTableJSONByteLength )
    {
        char* jsob_buffer = new char[jsob_buffer_size];
        b3dm.read( jsob_buffer, jsob_buffer_size );
        featureTableJSON = std::string( jsob_buffer, jsob_buffer_size );

        CONSOLE_EVAL( featureTableJSON );
        delete[] jsob_buffer;
    }

    if ( auto length = hdr.featureTableBinaryByteLength )
    {
        char* buffer = new char[length];
        b3dm.read( buffer, length );
        delete[] buffer;
    }

    if ( auto length = hdr.batchTableJSONByteLength )
    {
        char* buffer = new char[length];
        b3dm.read( buffer, length );
        delete[] buffer;
    }

    if ( auto length = hdr.batchTableBinaryByteLength )
    {
        char* buffer = new char[length];
        b3dm.read( buffer, length );
        delete[] buffer;
    }

    auto pos = b3dm.tellg();
    CONSOLE_EVAL( pos );

    if ( pos >= hdr.byteLength )
        return false_because( "glTF byte length must be greater than 0." );

    {
        size_t length = hdr.byteLength - pos;
        binary_glTF.resize( length );
        b3dm.read( binary_glTF.data(), length );
    }

    return true;
}


/// @brief Duplicated from Runtime_Tests/Externals_AssImp_Test.cpp
///
class AssImpF : public GlobalMapperSDKTestF
{
protected:
    /// @brief Search for the assimp test model directory
    /// @return 
    fs::path modelsDirectory() const
    {
        const char* d[] = {
            R"(C:\Users\yyk\Downloads\Assimp\models)",
            R"(V:\Asset Importer Test Models\models)",
        };

        auto pos = std::find_if( d, d + sizeof( d ) / sizeof( *d ), []( auto s ) { return fs::is_directory( s ); } );

        return fs::path( pos ? *pos : "<not-found>" );
    }

    fs::path cesiumModelsDirectory() const
    {
        // TODO: check for existence of a local data
        //       to speedup tests
        return fs::path( R"(V:\cesium_3d_tiles)" );
    }

    fs::path workspace_location() const
    {
        return fs::path( "out" ) / test_case_name() / test_name();
    }

    fs::path create_workspace()
    {
        auto ws = workspace_location();
        if ( fs::is_directory( ws ) )
            fs::remove_all( ws );
        std::error_code err;
        fs::create_directories( ws, err );

        CONSOLE_T( "ws = " << ws );

        return ws;
    }

    bool save_assets_as( aiScene const* scene, std::string base_name )
    {
        CONSOLE_TE( base_name );
        if ( !scene )
            return false;

        auto ws = workspace_location();
        if (0) { // verbose?
            for ( int i = 0; i != scene->mNumMeshes; ++i )
            {
                CONSOLE_T( i << " : " << scene->mMeshes[i]->mName.C_Str() );
            }
        }

        {
            auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
            std::string f = ( ws / ( base_name + ".xml" ) ).string();
            Assimp::Exporter exp;
            if ( !aiReturn_SUCCESS, exp.Export( scene, "assxml", f, flags ) )
                return false;
        }
        {
            auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
            std::string f = ( ws / ( base_name + ".glb" ) ).string();
            Assimp::Exporter exp;
            if ( !aiReturn_SUCCESS, exp.Export( scene, "glb2", f, flags ) )
                return false;
        }
        return true;
    }

    bool save_assets_as_glb2( aiScene const* scene, std::string base_name )
    {
        CONSOLE_TE( base_name );
        if ( !scene )
            return false;

        auto ws = workspace_location();
        {
            auto flags = 0
                | aiProcess_GenNormals
                | aiProcess_ValidateDataStructure
                ;
            std::string f = ( ws / ( base_name + ".glb" ) ).string();
            Assimp::Exporter exp;
            if ( !aiReturn_SUCCESS, exp.Export( scene, "glb2", f, flags ) )
                return false;
        }
        return true;
    }

    /// <summary>
    /// Find first existing file in the list
    /// </summary>
    /// <param name="files"></param>
    /// <returns>filename or empty string</returns>
    std::string find_file( std::vector<std::string> const& files )
    {
        for ( auto f : files ) if ( fs::is_regular_file( f ) )
            return f;
        return std::string{};
    }
};

/// Related to GM-13206 (Add Export Support for Cesium 3D Tiles)
/// Currently fails.
/// </summary>
/// <param name="">--gtest_also_run_disabled_tests --gtest_filter=AssImpF.DISABLED_GM_13206_t1</param>
/// <param name=""></param>
TEST_F( AssImpF, DISABLED_GM_13206_t1 )
{
    const char* filename = R"(Z:\GM-13206\CF_BUILDING\model.glb)";

    if ( !fs::is_regular_file( filename ) )
        GTEST_SKIP() << "file not found " << filename;

    {
        // create a logger from both CPP
        Assimp::DefaultLogger::create( "out/AssimpLog_Cpp.txt", Assimp::Logger::VERBOSE,
                aiDefaultLogStream_STDOUT | aiDefaultLogStream_DEBUGGER | aiDefaultLogStream_FILE );

        // .. and C. They should smoothly work together
        aiEnableVerboseLogging( AI_TRUE );
        aiLogStream logstream = aiGetPredefinedLogStream( aiDefaultLogStream_FILE, "out/AssimpLog_C.txt" );
        aiAttachLogStream( &logstream );


        Assimp::Importer importer;
        auto actual = importer.ReadFile( filename, 0 );
        EXPECT_FALSE( actual != nullptr ) << filename; // TODO: FIX!

        if ( actual == nullptr )
            CONSOLE_RED( "The test expected to FAIL. Please FIX" );

        aiDetachAllLogStreams();
    }
}

/// <summary>
/// Related to GM-13206 (Add Export Support for Cesium 3D Tiles)
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t2_1</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t2_1 )
{
    const char* filename = R"(Z:\GM-13206\CF_BUILDING\CF_BUILDING.glb)";

    if ( !fs::is_regular_file( filename ) )
        GTEST_SKIP() << "file not found " << filename;

    {
        // create a logger from both CPP
        Assimp::DefaultLogger::create( "out/AssimpLog_Cpp.txt", Assimp::Logger::VERBOSE,
                aiDefaultLogStream_STDOUT | aiDefaultLogStream_DEBUGGER | aiDefaultLogStream_FILE );

        // .. and C. They should smoothly work together
        aiEnableVerboseLogging( AI_TRUE );
        aiLogStream logstream = aiGetPredefinedLogStream( aiDefaultLogStream_FILE, "out/AssimpLog_C.txt" );
        aiAttachLogStream( &logstream );


        Assimp::Importer importer;
        auto actual = importer.ReadFile( filename, 0 );
        EXPECT_TRUE( actual != nullptr ) << filename;

        aiDetachAllLogStreams();
    }
}


/// <summary>
/// Import some tiles and dump them as .glb files
/// Related to: Add Export Support for Cesium 3D Tiles
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t0</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t0 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    {
        auto model_filename = ( model_dir / "BoxTextured_3d_tiles/0composite0.b3dm" ).string();
        Batched3DModel actual;

        ASSERT_TRUE( actual.read( model_filename ) )
            << "Failed: " << model_filename;

        Assimp::Importer imp;
        auto model = imp.ReadFileFromMemory( actual.binary_glTF.data(), actual.binary_glTF.size(), 0 );
        ASSERT_TRUE( model != nullptr );
        Assimp::Exporter exp;
        auto flags =
            aiProcess_GenNormals |
            aiProcess_ValidateDataStructure |
            0;
        {
            std::string filename_glb = ( ws / "BoxTextured_3d_tiles.glb" ).string();
            EXPECT_EQ( aiReturn_SUCCESS, exp.Export( model, "glb2", filename_glb, flags ) )
                << "Failed export to " << filename_glb;
        }
        {
            std::string filename_glb = ( ws / "BoxTextured_3d_tiles.assxml" ).string();
            EXPECT_EQ( aiReturn_SUCCESS, exp.Export( model, "assxml", filename_glb, flags ) )
                << "Failed export to " << filename_glb;
        }
    }
    {
        auto model_filename = ( model_dir / "CF_BUILDING_3D_Tiles/0material19_03.b3dm" ).string();
        Batched3DModel actual;

        ASSERT_TRUE( actual.read( model_filename ) )
            << "Failed: " << model_filename;
        Assimp::Importer imp;
        auto model = imp.ReadFileFromMemory( actual.binary_glTF.data(), actual.binary_glTF.size(), 0 );
        ASSERT_TRUE( model != nullptr );
        Assimp::Exporter exp;
        auto flags =
            aiProcess_GenNormals |
            aiProcess_ValidateDataStructure |
            0;
        std::string filename_glb = ( ws / "CF_BUILDING_3D_Tiles.glb" ).string();
        EXPECT_EQ( aiReturn_SUCCESS, exp.Export( model, "glb2", filename_glb, flags ) )
            << "Failed export to " << filename_glb;
    }
}

/// <summary>
/// Import CF_BUILDING_3D_Tiles tiles and dump them as .glb files
/// Related to: Add Export Support for Cesium 3D Tiles
/// </summary>
/// <param name="">--gtest_also_run_disabled_tests --gtest_filter=AssImpF.GM13206_t1</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t1 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    const char* tiles[] = {
        "0material0_0.b3dm",
        "0material0_00.b3dm",
        "0material0_001.b3dm",
        "0material14_0.b3dm",
        "0material14_00.b3dm",
        "0material14_000.b3dm",
        "0material14_001.b3dm",
        "0material14_0010.b3dm",
        "0material14_0011.b3dm",
        "0material14_01.b3dm",
        "0material14_010.b3dm",
        "0material14_0100.b3dm",
        "0material14_0101.b3dm",
        "0material14_011.b3dm",
        "0material14_02.b3dm",
        "0material16_0.b3dm",
        "0material16_00.b3dm",
        "0material16_000.b3dm",
        "0material16_001.b3dm",
        "0material16_002.b3dm",
        "0material16_003.b3dm",
        "0material16_01.b3dm",
        "0material19_0.b3dm",
        "0material19_00.b3dm",
        "0material19_01.b3dm",
        "0material19_02.b3dm",
        "0material19_03.b3dm",
        "0material6_0.b3dm",
        "0material6_00.b3dm",
        "0material6_01.b3dm",
    };
    auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure;

    for ( auto tile : tiles )
    {
        CONSOLE_EVAL( tile );
        auto model_filename = ( model_dir / "CF_BUILDING_3D_Tiles" / tile ).string();
        Batched3DModel actual;
        ASSERT_TRUE( actual.read( model_filename ) ) << "Failed: " << model_filename;
        Assimp::Importer imp;
        auto model = imp.ReadFileFromMemory( actual.binary_glTF.data(), actual.binary_glTF.size(), 0 );
        ASSERT_TRUE( model != nullptr );
        Assimp::Exporter exp;
        std::string filename_glb = ( ws / ( std::string( tile ) + ".glb" ) ).string();
        EXPECT_EQ( aiReturn_SUCCESS, exp.Export( model, "glb2", filename_glb, flags ) )
            << "Failed export to " << filename_glb;
    }
}



/// @brief Load and dump BoxTextured_3d_tiles/tileset.json
/// @param --gtest_filter=AssImpF.GM13206_t2
TEST_F( AssImpF, GM13206_t2 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    {
        auto tileset_filename = ( model_dir / "BoxTextured_3d_tiles/tileset.json" ).string();
        TilesetJson actual( tileset_filename );

        {
            CONSOLE_EVAL( actual.asset().version );
            EXPECT_STREQ( "1.0", actual.asset().version.c_str() );

            CONSOLE_EVAL( actual.geometricError() );
            EXPECT_DOUBLE_EQ( 1.7320508075688772, actual.geometricError() );
        }
        {
            auto root_tile = actual.root();
            EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
            EXPECT_DOUBLE_EQ( 1.7320508075688772, root_tile.geometricError );
            EXPECT_TRUE( (bool)root_tile.boundingVolume );
            EXPECT_DOUBLE_EQ( 1.4120970059478493, root_tile.transform->matrix4x4[0] );
            EXPECT_DOUBLE_EQ( 1.0, root_tile.transform->matrix4x4[15] );

            CONSOLE_EVAL( root_tile.children.size() );
            ASSERT_EQ( 1, root_tile.children.size() );
            ASSERT_EQ( 1, root_tile.children.front().children.size() );
            CONSOLE_EVAL( root_tile.children.front().children.front().content.uri );
            ASSERT_EQ( 0, root_tile.children.front().children.front().children.size() );
        }
    }
}

/// <summary>
/// Import CF_BUILDING_3D_Tiles tileset.json file and dump it
/// Related to: Add Export Support for Cesium 3D Tiles
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t3</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t3 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    {
        auto tileset_filename = ( model_dir / "CF_BUILDING_3D_Tiles/tileset.json" ).string();
        TilesetJson actual( tileset_filename );

        {
            CONSOLE_EVAL( actual.asset().version );
            EXPECT_STREQ( "1.0", actual.asset().version.c_str() );

            CONSOLE_EVAL( actual.geometricError() );
            EXPECT_DOUBLE_EQ( 204.46812474527445, actual.geometricError() );
        }
        {
            auto root_tile = actual.root();
            EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
            EXPECT_DOUBLE_EQ( 204.46812474527445, root_tile.geometricError );
            EXPECT_TRUE( (bool)root_tile.boundingVolume );

            auto actual_box = dynamic_cast<TilesetJson::boundingVolumeBox_t*>( root_tile.boundingVolume.get() );
            EXPECT_TRUE( (bool)actual_box );
            if ( actual_box )
            {
                std::unique_ptr<aiScene> model( actual.construct_scene( root_tile.boundingVolume.get() ) );

                {
                    auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
                    std::string f = ( ws / "construct_mesh.assxml" ).string();
                    Assimp::Exporter exp;
                    EXPECT_EQ( aiReturn_SUCCESS, exp.Export( model.get(), "assxml", f, flags ) );
                }
                {
                    auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
                    std::string f = ( ws / "construct_mesh.glb" ).string();
                    Assimp::Exporter exp;
                    EXPECT_EQ( aiReturn_SUCCESS, exp.Export( model.get(), "glb2", f, flags ) );
                    // verify the just exported model
                    {
                        Assimp::Importer imp;
                        auto actual = imp.ReadFile( f, aiProcess_ValidateDataStructure );
                        ASSERT_TRUE( (bool)actual );
                        EXPECT_EQ( 1, actual->mNumMeshes );
                        EXPECT_EQ( 8, actual->mMeshes[0]->mNumVertices );
                    }
                }
            }

            EXPECT_DOUBLE_EQ( 0.08482280432458418, root_tile.transform->matrix4x4[0] );
            EXPECT_DOUBLE_EQ( 1.0, root_tile.transform->matrix4x4[15] );

            CONSOLE_EVAL( root_tile.children.size() );
            ASSERT_EQ( 6, root_tile.children.size() );
            ASSERT_EQ( 5, root_tile.children.front().children.size() );
            ASSERT_EQ( 1, root_tile.children.back().children.size() );
        }
    }
}

/// <summary>
/// Import CF_BUILDING_3D_Tiles tileset.json file and dump it
/// Related to: Add Export Support for Cesium 3D Tiles
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t4</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t4 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    {
        auto tileset_filename = ( model_dir / "CF_BUILDING_3D_Tiles/tileset.json" ).string();
        TilesetJson actual( tileset_filename );


        {
            auto root_tile = actual.root();
            EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
            EXPECT_DOUBLE_EQ( 204.46812474527445, root_tile.geometricError );
            EXPECT_TRUE( (bool)root_tile.boundingVolume );

            auto root_volume = root_tile.boundingVolume.get();
            std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

            actual.append_node( scene.get(), root_tile.children.front().boundingVolume.get() );
            EXPECT_EQ( 2, scene->mNumMeshes );
            EXPECT_EQ( 2, scene->mRootNode->mNumChildren );

            EXPECT_TRUE( save_assets_as( scene.get(), "append_node" ) );
        }
    }
}

/// <summary>
/// Import CF_BUILDING_3D_Tiles tileset.json file and dump it
/// Related to: Add Export Support for Cesium 3D Tiles
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t5</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t5 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    {
        auto tileset_filename = ( model_dir / "CF_BUILDING_3D_Tiles/tileset.json" ).string();
        TilesetJson actual( tileset_filename );

        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 204.46812474527445, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get() );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "bounding_boxes" ) );
    }
}

/// <summary>
/// Import 2CylinderEngine_3D_tiles tileset.json file and dump it
/// Related to: Add Export Support for Cesium 3D Tiles
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t6</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t6 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    {
        auto tileset_filename = ( model_dir / R"(2CylinderEngine_3D_tiles\tileset.json)" ).string();
        if ( !fs::is_regular_file( tileset_filename ) )
            GTEST_SKIP() << "File not found: " << tileset_filename;

        TilesetJson actual( tileset_filename );

        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 404.740140685585, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get() );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( actual, root_tile, scene.get() );
        EXPECT_TRUE( save_assets_as( scene.get(), "bounding_boxes" ) );
    }
}

class MeshSanitizer
{
public:
    MeshSanitizer( aiScene const* scene )
        : m_scene( scene )
        , m_verbose( true )
    {
        if ( m_scene == nullptr )
            throw std::runtime_error( "scene is nullptr" );
    }

    void verbose( bool v )
    {
        m_verbose = v;
    }

    static bool is_close_enough( double a, double b )
    {
        return std::abs( a - b ) < 1E-15;
    }

    static bool is_good_triangle( aiVector3D const& p0, aiVector3D const& p1, aiVector3D const& p2 )
    {
        bool ok = true;

        double d[] = {
            ( p1 - p0 ).Length(),
            ( p2 - p0 ).Length(),
            ( p2 - p1 ).Length()
        };
        std::sort( d, d + 3 );
        if ( is_close_enough( d[0] + d[1], d[2] ) )
            ok = false;

        return ok;
    }

    static bool is_good_triangle( unsigned int i, unsigned int j, unsigned int k, aiVector3D const* v )
    {
        if ( i == j || j == k || i == k )
        {
            //CONSOLE( "BAD INDICES: i,j,k: " << i << ", " << j << ", " << k );
            return false;
        }
        return is_good_triangle( v[i], v[j], v[k] );
    }

    unsigned int sanitize_mesh( aiMesh const* mesh )
    {
        unsigned int cnt = 0;
        if ( mesh->mPrimitiveTypes & aiPrimitiveType_TRIANGLE )
        {
            for ( int i = 0; i != mesh->mNumFaces; ++i )
            {
                aiFace const& face = mesh->mFaces[i];
                if ( face.mNumIndices == 3 )
                {
                    bool ok = is_good_triangle(
                        face.mIndices[0],
                        face.mIndices[1],
                        face.mIndices[2],
                        mesh->mVertices );
                    if ( !ok )
                    {
                        //CONSOLE( "DEGENERATE triangle: " << mesh->mName.C_Str() << ", i=" << i );
                        ++cnt;
                    }
                }
                else
                {
                    CONSOLE( "face.mNumIndices == " << face.mNumIndices );
                }
            }
        }
        else
        {
            CONSOLE( mesh->mName.C_Str() << " has no aiPrimitiveType_TRIANGLE" );
        }
        return cnt;
    }

    void navigate_through( aiNode const* node )
    {
        if ( m_verbose ) CONSOLE( node->mName.C_Str() );
        for ( int i = 0; i != node->mNumMeshes; ++i )
        {
            aiMesh const* mesh = m_scene->mMeshes[node->mMeshes[i]];
            if ( m_verbose )
                CONSOLE( "mesh->mName: " << mesh->mName.C_Str() );
            auto cnt = sanitize_mesh( mesh );
            if ( cnt )
            {
                CONSOLE( mesh->mName.C_Str() << " has " << cnt << " degenerate triangle(s)" );
            }
        }
        for ( int i = 0; i != node->mNumChildren; ++i )
        {
            navigate_through( node->mChildren[i] );
        }
    }

    void summary()
    {
        if ( !m_scene->HasMeshes() )
        {
            std::cout << "Scene does not have meshes\n";
            return;
        }
        auto root = m_scene->mRootNode;
        navigate_through( root );
    }

protected:
    aiScene const* m_scene;
    bool m_verbose;
};

/// <summary>
/// Example of .glb file with many degenerate faces
/// </summary>
/// <param name="">--gtest_filter=AssImpF.2CylinderEngine_SANITIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, 2CylinderEngine_SANITIZE )
{
    std::string model_filename = find_file( {
        R"(C:\Users\yyk\Downloads\2CylinderEngine.glb)",
        R"(c:\Users\yurikuznetsov\Downloads\2CylinderEngine.glb)"
        } );
    if ( model_filename.empty() )
        GTEST_SKIP() << "model file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    MeshSanitizer msan( actual );

    msan.verbose( false );
    msan.summary();

    EXPECT_TRUE( save_assets_as( actual, "2CylinderEngine-out" ) );
}

/// <summary>
/// The model with corrected faces
/// </summary>
/// <param name="">--gtest_filter=AssImpF.2CylinderEngineGOOD_SANITIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, 2CylinderEngineGOOD_SANITIZE )
{
    std::string model_filename = find_file( {
        R"(C:\Users\yurikuznetsov\Downloads\2CylinderEngine (good).glb)",
        R"(Z:\GM-13206\2CylinderEngine (good).glb)"
        } );
    if ( model_filename.empty() )
        GTEST_SKIP() << "model file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    MeshSanitizer msan( actual );

    msan.verbose( false );
    msan.summary();

    EXPECT_TRUE( save_assets_as( actual, "2CylinderEngine-out" ) );
}

/// <summary>
/// 
/// </summary>
/// <param name="">--gtest_filter=AssImpF.BoxTextured</param>
/// <param name=""></param>
TEST_F( AssImpF, BoxTextured )
{
    auto model_dir = modelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "models directory not available: " << model_dir;

    auto ws = create_workspace();
    //Assimp::DefaultLogger::create( "out/AssimpLog_Cpp.txt", Assimp::Logger::VERBOSE,
    //        aiDefaultLogStream_STDOUT | aiDefaultLogStream_DEBUGGER | aiDefaultLogStream_FILE );
    //aiEnableVerboseLogging( AI_TRUE );
    //aiLogStream logstream = aiGetPredefinedLogStream( aiDefaultLogStream_FILE, "out/AssimpLog_C.txt" );
    //aiAttachLogStream( &logstream );

    auto f = R"(glTF2\BoxTextured-glTF\BoxTextured.gltf)";
    auto file = ( model_dir / f ).string();
    Assimp::Importer importer;
    auto actual = importer.ReadFile( file, 0 );
    EXPECT_TRUE( actual != nullptr ) << file;

    {
        auto flags = /*aiProcess_GenNormals | */aiProcess_ValidateDataStructure | 0;
        std::string f = ( ws / "BoxTextured.assxml" ).string();
        Assimp::Exporter exp;
        EXPECT_EQ( aiReturn_SUCCESS, exp.Export( actual, "assxml", f, flags ) );
    }
    {
        auto flags = /*aiProcess_GenNormals | */aiProcess_ValidateDataStructure | 0;
        std::string f = ( ws / "BoxTextured.glb" ).string();
        Assimp::Exporter exp;
        EXPECT_EQ( aiReturn_SUCCESS, exp.Export( actual, "glb2", f, flags ) );
    }
}

/// <summary>
/// 
/// </summary>
/// <param name="">--gtest_filter=AssImpF.Duck_SANITIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, Duck_SANITIZE )
{
    const char* model_filenames[] = {
        R"(C:\Users\yurikuznetsov\Downloads\Duck.glb)",
        R"(C:\Users\yyk\Downloads\Duck.glb)",
    };
    const char* model_filename = 0;
    for ( auto p : model_filenames ) if ( fs::is_regular_file( p ) )
    {
        model_filename = p;
        break;
    }
    if ( !model_filename )
        GTEST_SKIP() << "models file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    MeshSanitizer msan( actual );

    msan.summary();

    EXPECT_TRUE( save_assets_as( actual, "Duck-out" ) );
}

/// <summary>
/// c/Users/yyk/Downloads/CF_BUILDING.glb
/// </summary>
/// <param name="">--gtest_filter=AssImpF.CF_BUILDING_SANITIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, CF_BUILDING_SANITIZE )
{
    std::string model_filename = find_file( {
        R"(C:\Users\yurikuznetsov\Downloads\CF_BUILDING.glb)",
        R"(C:\Users\yyk\Downloads\CF_BUILDING.glb)",
    } );
    if ( model_filename.empty() )
        GTEST_SKIP() << "model file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    MeshSanitizer msan( actual );
    msan.verbose( false );
    msan.summary();

    EXPECT_TRUE( save_assets_as( actual, "CF_BUILDING-out" ) );
}

TEST_F( AssImpF, self_SANITIZE )
{
    std::cout << "sizeof(std::unique_ptr<int>) == " << sizeof( std::unique_ptr<int> ) << std::endl;
    std::cout << "sizeof(int *) == " << sizeof( int* ) << std::endl;
    {
        std::unique_ptr<int> foo{ new int {123} };
        printf( "%d\n", *foo );
    }
    {
        int* foo{ new int {123} };
        printf( "%d\n", *foo );
        delete foo;
    }
}


/// <summary>
/// Load P2P generated mesh (Z-coordinate is UP)
/// </summary>
/// <param name="">--gtest_filter=AssImpF.Mike_Niwot_mesh_SANITIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, Mike_Niwot_mesh_SANITIZE )
{
    std::string model_filename = find_file( {
        R"(V:\cesium_3d_tiles\Mike_Niwot_quick_subset_Zup\Mike_Niwot_quick_subset_Zup.glb)",
        } );
    if ( model_filename.empty() )
        GTEST_SKIP() << "model file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    EXPECT_TRUE( actual->HasMaterials() );
    EXPECT_FALSE( actual->HasTextures() );
    EXPECT_TRUE( actual->HasMeshes() );

    MeshSanitizer msan( actual );

    msan.verbose( false );
    msan.summary();

    EXPECT_TRUE( save_assets_as( actual, "actual" ) );
}

/// @brief MeshNormalizer is a class to "normalize" mesh
///
/// The normalization includes a few mesh transformations:
///     - re-projection
///     ....
class MeshNormalizer
{
public:
    aiScene const* normalize( aiScene const* actual)
    {
        aiScene* norm = 0;
        Assimp::SceneCombiner::CopyScene( &norm, actual );
        if ( is_projected() )
        {
            aiAABB bb = scene_bounding_box( norm );
            for ( unsigned int i = 0; i < norm->mNumMeshes; ++i )
            {
                aiMesh* mesh = norm->mMeshes[i];
                if ( nullptr == mesh )
                {
                    continue;
                }

                for ( unsigned int i = 0; i < mesh->mNumVertices; ++i )
                {
                    mesh->mVertices[i] -= bb.mMin;
                }
            }
        }

        return norm;
    }

    /// @brief Compute the bounding box for the entire scene
    /// @param pScene 
    /// @return the bounding box
    aiAABB scene_bounding_box( aiScene const* pScene )
    {
        ai_real constexpr M = std::numeric_limits<ai_real>::max();
        aiVector3D min( M, M, M ), max( -M, -M, -M );

        for ( unsigned int i = 0; i < pScene->mNumMeshes; ++i )
        {
            aiMesh* mesh = pScene->mMeshes[i];
            if ( nullptr == mesh )
            {
                continue;
            }

            checkMesh( mesh, min, max );
        }
        return { min, max };
    }
protected:
    bool is_projected()
    {
        return m_projection
            // TODO: && is_projected_projection( m_projection )
            && 1;
    }
    void checkMesh( aiMesh* mesh, aiVector3D& min, aiVector3D& max )
    {
        ai_assert( nullptr != mesh );

        if ( 0 == mesh->mNumVertices )
        {
            return;
        }

        for ( unsigned int i = 0; i < mesh->mNumVertices; ++i )
        {
            const aiVector3D& pos = mesh->mVertices[i];
            if ( pos.x < min.x )
            {
                min.x = pos.x;
            }
            if ( pos.y < min.y )
            {
                min.y = pos.y;
            }
            if ( pos.z < min.z )
            {
                min.z = pos.z;
            }

            if ( pos.x > max.x )
            {
                max.x = pos.x;
            }
            if ( pos.y > max.y )
            {
                max.y = pos.y;
            }
            if ( pos.z > max.z )
            {
                max.z = pos.z;
            }
        }
    }

public:
    GM_Projection_t const* m_projection{};
};

/// <summary>
/// NORMALIZE P2P generated mesh (Z-coordinate is UP)
/// </summary>
/// <param name="">--gtest_filter=AssImpF.Mike_Niwot_mesh_NORMALIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, Mike_Niwot_mesh_NORMALIZE )
{
    std::string model_filename = find_file( {
        R"(V:\cesium_3d_tiles\Mike_Niwot_quick_subset_Zup\Mike_Niwot_quick_subset_Zup.glb)",
        } );
    if ( model_filename.empty() )
        GTEST_SKIP() << "model file not available";
    std::string prj_filename = find_file( {
        R"(V:\cesium_3d_tiles\Mike_Niwot_quick_subset_Zup\Mike_Niwot_quick_subset_Zup.prj)",
        } );
    if ( prj_filename.empty() )
        GTEST_SKIP() << "project file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    EXPECT_TRUE( actual->HasMaterials() );
    EXPECT_FALSE( actual->HasTextures() );
    EXPECT_TRUE( actual->HasMeshes() );
    EXPECT_TRUE( save_assets_as( actual, "actual" ) );

    GM_Projection_t projection;
    ASSERT_EQ( 0, GM_LoadProjectionFile( prj_filename.c_str(), &projection ) );
    MeshNormalizer mnorm{ &projection };

    {
        GM_Projection_t geom_proj;
        ASSERT_EQ( 0, GM_LoadProjectionFromEPSGCode( 4326, &geom_proj ) );
        auto actual_bb = mnorm.scene_bounding_box( actual );
        // actual_bb is a mesh point, where Y is UP direction.
        //  x = x, y = -z, z = y;

        CONSOLE_TE( actual_bb );
        // actual_bb: { { 496791.25, 1624.61084, -4420783 }, { 496904.688,1637.79309,-4420651 } }

        GM_Point_t p0, p1;
        GM_ProjectPoint( actual_bb.mMin.x, -actual_bb.mMin.z, &p0.mX, &p0.mY,
                         &projection, &geom_proj );
        GM_ProjectPoint( actual_bb.mMax.x, -actual_bb.mMax.z, &p1.mX, &p1.mY,
                         &projection, &geom_proj );

        CONSOLE_T( "min coordinate: " << p0 );
        CONSOLE_T( "max coordinate: " << p1 );
        // min coordinate : {-105.037556518092288, 39.937157237963973}
        // max coordinate : {-105.036228172496763, 39.935968353625384}

        // From V:\cesium_3d_tiles\Mike_Niwot_mesh_NORMALIZE_3D_tiles\tileset.json
        // The transform from the currently "normalized" mesh to its original
        // geo location at {-105.037556518092288, 39.937157237963973}
        // in geocentric coordinates (EPSG:4978)
        // These numbers correspond to the geo reference:
        // (-105.037556518092288, 39.937157237963973, 1624.6099999999)
        double transform[] = {
             0.965755966816224, -0.259452139245103, -2.77555756156289E-17, 0,
             0.166554694304208,  0.619964785387972,  0.76674845854692,     0,
            -0.198934527832884, -0.74049189888883,   0.641947662442915,    0,
            -1270909.57083419,  -4730693.31741936,   4073680.81026208,     1
        };

        double  ecef_X, ecef_Y, ecef_Z;
        {
            auto rc = GM_ProjectPointToECEF(
                -105.037556518092288, 39.937157237963973, 1624.6099999999,
                &ecef_X, &ecef_Y, &ecef_Z,
                GM_DATUM_WGS_84, &geom_proj );
            EXPECT_EQ( 0, rc ) << "Error: " << GM_Error( rc );
        }
        CONSOLE_T( "ecef_X, ecef_Y, ecef_Z: " << ecef_X << ", " << ecef_Y << ", " << ecef_Z );
    }

    auto actual_norm = mnorm.normalize(actual);
    EXPECT_TRUE( save_assets_as( actual_norm, "actual_norm" ) );
}

// "V:\cesium_3d_tiles\Mike_Niwot_quick_subset_Zup\Mike_Niwot_quick_subset_Zup.glb"

/// <summary>
/// Load P2P generated mesh (Z-coordinate is UP,
/// the origin moved to left upper corner)
/// </summary>
/// <param name="">--gtest_filter=AssImpF.Mike_Niwot_mesh_Zero_SANITIZE</param>
/// <param name=""></param>
TEST_F( AssImpF, Mike_Niwot_mesh_Zero_SANITIZE )
{
    std::string model_filename = find_file( {
        R"(Z:\GM-13206\Mike_Niwot_3D\Mike_Niwot_mesh_Zero.glb)",
        } );
    if ( model_filename.empty() )
        GTEST_SKIP() << "model file not available";

    auto ws = create_workspace();

    Assimp::Importer importer;
    CONSOLE_EVAL( model_filename );
    auto actual = importer.ReadFile( model_filename, 0 );
    ASSERT_TRUE( actual != nullptr ) << model_filename;

    MeshSanitizer msan( actual );

    msan.verbose( false );
    msan.summary();

    EXPECT_TRUE( save_assets_as( actual, "actual" ) );
}

/// <summary>
/// Dump bounding volumes (boxes)
/// The json file was generated by Cesium/ION
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t7</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t7 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    // Generated by Cesium/ION
    // "V:\cesium_3d_tiles\Mike_Niwot_mesh\tileset.json"

    auto tileset_filename = ( model_dir / "Mike_Niwot_mesh/tileset.json" ).string();
    TilesetJson actual( tileset_filename );

    {
        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 31.823295461171298, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get(), actual.make_mesh_name(*cp) );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "bounding_boxes" ) );

        // Step 2. Dump the "tiles" as .GLB files
        const char* tile_filenames[] = {
            "Mike_Niwot_mesh/0material0_0.b3dm",
            "Mike_Niwot_mesh/0material0_00.b3dm",
            "Mike_Niwot_mesh/0material0_01.b3dm",
            "Mike_Niwot_mesh/0material0_02.b3dm",
            "Mike_Niwot_mesh/0material0_03.b3dm",
        };

        for ( auto tile : tile_filenames )
        {
            auto tile_filename = ( model_dir / tile );

            Batched3DModel d3dm;
            ASSERT_TRUE( d3dm.read( tile_filename.string() ) ) << "File: " << tile_filename;

            Assimp::Importer importer;
            auto scene = importer.ReadFileFromMemory( d3dm.binary_glTF.data(), d3dm.binary_glTF.size(), 0 );
            ASSERT_TRUE( scene != nullptr );
            EXPECT_TRUE( save_assets_as( scene, tile_filename.stem().string() ) );
        }
    }
}

/// <summary>
/// Dump bounding volumes (boxes) from CV_Example_3D_tiles
/// The json file was generated by Cesium/ION
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t8</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t8 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    // Generated by Cesium/ION
    // "V:\cesium_3d_tiles\CV_Example_3D_tiles\tileset.json"

    auto tileset_filename = ( model_dir / "CV_Example_3D_tiles/tileset.json" ).string();
    TilesetJson actual( tileset_filename );

    {
        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 2436.9730410411607, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get(), actual.make_mesh_name( *cp ) );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "bounding_boxes" ) );

        // Step 2. Dump the "tiles" as .GLB files
        const char* tile_filenames[] = {
//            "CV_Example_3D_tiles/0composite0.cmpt",
            "CV_Example_3D_tiles/0material1_0.b3dm",
            "CV_Example_3D_tiles/0material1_00.b3dm",
            "CV_Example_3D_tiles/0material1_000.b3dm",
            "CV_Example_3D_tiles/0material1_0000.b3dm",
            "CV_Example_3D_tiles/0material1_0001.b3dm",
            "CV_Example_3D_tiles/0material1_001.b3dm",
            "CV_Example_3D_tiles/0material1_0010.b3dm",
            "CV_Example_3D_tiles/0material1_0011.b3dm",
            "CV_Example_3D_tiles/0material1_01.b3dm",
            "CV_Example_3D_tiles/0material1_010.b3dm",
            "CV_Example_3D_tiles/0material1_0100.b3dm",
            "CV_Example_3D_tiles/0material1_0101.b3dm",
            "CV_Example_3D_tiles/0material1_012.b3dm",
            "CV_Example_3D_tiles/0material1_013.b3dm",
            "CV_Example_3D_tiles/0material1_0130.b3dm",
            "CV_Example_3D_tiles/0material1_0131.b3dm",
            "CV_Example_3D_tiles/0material1_01310.b3dm",
            "CV_Example_3D_tiles/0material1_01311.b3dm",
            "CV_Example_3D_tiles/0material1_01312.b3dm",
//            "CV_Example_3D_tiles/tileset.json",
        };

        for ( auto tile : tile_filenames )
        {
            auto tile_filename = ( model_dir / tile );

            Batched3DModel d3dm;
            ASSERT_TRUE( d3dm.read( tile_filename.string() ) ) << "File: " << tile_filename;

            Assimp::Importer importer;
            auto scene = importer.ReadFileFromMemory( d3dm.binary_glTF.data(), d3dm.binary_glTF.size(), 0 );
            ASSERT_TRUE( scene != nullptr );
            EXPECT_TRUE( save_assets_as( scene, tile_filename.stem().string() ) );
        }
    }
}

/// <summary>
/// Dump bounding volumes (boxes)
/// The json file was generated by Cesium/ION (GM13206_test.tc_3)
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t9</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t9 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    // Generated by Cesium/ION
    // "V:\cesium_3d_tiles\GM13206_test.tc_3\tileset.json"

    auto tileset_filename = ( model_dir / "GM13206_test.tc_3/tileset.json" ).string();
    TilesetJson actual( tileset_filename );

    {
        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 519.6152422706632, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get(), actual.make_mesh_name( *cp ) );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "bounding_boxes" ) );

        // Step 2. Dump the "tiles" as .GLB files
        const char* tile_filenames[] = {
            "GM13206_test.tc_3/0material0_0.b3dm",
            "GM13206_test.tc_3/0material0_01.b3dm",
//            "GM13206_test.tc_3/tileset.json",
        };

        for ( auto tile : tile_filenames )
        {
            auto tile_filename = ( model_dir / tile );

            Batched3DModel d3dm;
            ASSERT_TRUE( d3dm.read( tile_filename.string() ) ) << "File: " << tile_filename;

            Assimp::Importer importer;
            auto scene = importer.ReadFileFromMemory( d3dm.binary_glTF.data(), d3dm.binary_glTF.size(), 0 );
            ASSERT_TRUE( scene != nullptr );
            EXPECT_TRUE( save_assets_as( scene, tile_filename.stem().string() ) );
        }
    }
}

/// <summary>
/// Dump bounding volumes (boxes)
/// The json file was generated by 3d-tiles-tools
/// (see also files produced by AssImpF.GM13206_t9)
/// </summary>
/// <param name="">--gtest_filter=AssImpF.GM13206_t10</param>
/// <param name=""></param>
TEST_F( AssImpF, GM13206_t10 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;

    auto ws = create_workspace();

    // Generated by Cesium/ION
    // "V:\cesium_3d_tiles\GM13206_t9\tileset.json"

    auto tileset_filename = ( model_dir / "GM13206_t9/tileset.json" ).string();
    TilesetJson actual( tileset_filename );

    {
        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 1024, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
        node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
            {
                for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
                {
                    actual.append_node( scene, cp->boundingVolume.get(), actual.make_mesh_name( *cp ) );
                    node_appender( actual, *cp, scene );
                }
            };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "bounding_boxes" ) );
    }
}

/// @brief 
/// @param --gtest_filter=AssImpF.GM13206_t11
/// @param  
TEST_F( AssImpF, GM13206_t11 )
{
    auto model_dir = cesiumModelsDirectory();
    if ( !fs::is_directory( model_dir ) )
        GTEST_SKIP() << "Cesium model directory not available: " << model_dir;
    auto ws = create_workspace();

    // Generated by Cesium/ION
    // "V:\cesium_3d_tiles\Mike_Niwot_mesh_NORMALIZE_3D_tiles\tileset.json"
    auto sub_directory = "Mike_Niwot_mesh_NORMALIZE_3D_tiles";
    {
        auto tileset_filename = ( model_dir / sub_directory / "tileset.json" ).string();
        TilesetJson actual( tileset_filename );

        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 44.211640618984198, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto bb_scene = actual.build_bounding_boxes_scene();
        EXPECT_TRUE( save_assets_as( bb_scene.get(), "bounding_boxes" ) );
    }
    {
        // Step 2. Dump the "tiles" as .GLB files
        const char* tile_filenames[] = {
            "0material0_0.b3dm",
            "0material0_00.b3dm",
            "0material0_000.b3dm",
            "0material0_001.b3dm",
            "0material0_01.b3dm",
            "0material0_010.b3dm",
            "0material0_011.b3dm",
            "0material0_012.b3dm",
            "0material0_013.b3dm",
            "0material0_02.b3dm",
            "0material0_020.b3dm",
            "0material0_021.b3dm",
            "0material0_022.b3dm",
            "0material0_023.b3dm",
            "0material0_03.b3dm",
            "0material0_030.b3dm",
            "0material0_031.b3dm",
            "0material0_032.b3dm",
            "0material0_033.b3dm",
        };

        for ( auto tile : tile_filenames )
        {
            auto tile_filename = ( model_dir / sub_directory / tile );

            Batched3DModel d3dm;
            ASSERT_TRUE( d3dm.read( tile_filename.string() ) ) << "File: " << tile_filename;

            Assimp::Importer importer;
            auto scene = importer.ReadFileFromMemory( d3dm.binary_glTF.data(), d3dm.binary_glTF.size(), 0 );
            ASSERT_TRUE( scene != nullptr );
            EXPECT_TRUE( save_assets_as( scene, tile_filename.stem().string() ) );
        }
    }
}


/// @brief Load IFC file
/// 
/// Related to GM-17192: Unable to load IFC file V:\IFC\AC14-FZK-Haus.ifc
/// @param --gtest_filter=AssImpF.ifc_t0
/// @param  
TEST_F( AssImpF, ifc_t0 )
{
    auto ws = create_workspace();

    auto filename = R"(V:\IFC\AC14-FZK-Haus.ifc)";

    Assimp::Importer impo;
    aiScene const *actual = impo.ReadFile( filename, 0 );
    EXPECT_NE( nullptr, actual );

    save_assets_as_glb2( actual, "AC14-FZK-Haus.ifc" );
}

#pragma region vtk

/// @brief 
/// @param --gtest_filter=AssImpF.ifc_t1
/// @param  
TEST_F( AssImpF, ifc_t1 )
{
    auto ws = create_workspace();

    const char* files[] = {
        R"(V:\IFC\CV_Example.ifc)",
        //R"(V:\IFC\Wellness center Sama.ifc)", // bad?
        //R"(V:\IFC\Example_01_.ifc)", // bad?
    };

    for ( auto filename : files )
    {
        CONSOLE_TE( filename );

        Assimp::Importer impo;
        aiScene const* actual = impo.ReadFile( filename, 0 );
        EXPECT_NE( nullptr, actual );

        save_assets_as_glb2( actual, fs::path(filename).stem().string() );
    }
}

#include <vtkCesium3DTilesWriter.h>
#include <vtkNew.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkOBJReader.h>
#include <vtkLogger.h>
#include <vtkStringArray.h>
#include <vtkFieldData.h>

//------------------------------------------------------------------------------
std::array<double, 3> ReadOBJOffset( const char* comment )
{
    std::array<double, 3> translation = { 0, 0, 0 };
    if ( comment )
    {
        std::istringstream istr( comment );
        std::array<std::string, 3> axesNames = { "x", "y", "z" };
        for ( int i = 0; i < 3; ++i )
        {
            std::string axis;
            std::string s;
            istr >> axis >> s >> translation[i];
            if ( istr.fail() )
            {
                vtkLog( WARNING, "Cannot read axis " << axesNames[i] << " from comment." );
            }
            if ( axis != axesNames[i] )
            {
                vtkLog( WARNING, "Invalid axis " << axesNames[i] << ": " << axis );
            }
        }
    }
    else
    {
        vtkLog( WARNING, "nullptr comment." );
    }
    return translation;
}

vtkSmartPointer<vtkPolyData> ReadOBJMesh( int numberOfBuildings, int vtkNotUsed( lod ),
  const std::vector<std::string>& files, std::array<double, 3>& fileOffset )
{
    vtkNew<vtkAppendPolyData> append;
    for ( size_t i = 0; i < files.size() && i < static_cast<size_t>( numberOfBuildings ); ++i )
    {
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName( files[i].c_str() );
        reader->Update();
        if ( i == 0 )
        {
            fileOffset = ReadOBJOffset( reader->GetComment() );
        }
        auto polyData = reader->GetOutput();
        append->AddInputDataObject( polyData );
    }
    append->Update();
    return append->GetOutput();
}

//------------------------------------------------------------------------------
void SetField( vtkDataObject* obj, const char* name, const char* value )
{
    vtkFieldData* fd = obj->GetFieldData();
    if ( !fd )
    {
        vtkNew<vtkFieldData> newfd;
        obj->SetFieldData( newfd );
        fd = newfd;
    }
    vtkNew<vtkStringArray> sa;
    sa->SetNumberOfTuples( 1 );
    sa->SetValue( 0, value );
    sa->SetName( name );
    fd->AddArray( sa );
}

//------------------------------------------------------------------------------
std::string GetOBJTextureFileName( const std::string& file )
{
    std::string fileNoExt = fs::path( file ).stem().string();
    std::string textureFileName = fileNoExt + ".png";
    return fs::is_regular_file( textureFileName ) ? textureFileName : std::string();
}

vtkSmartPointer<vtkMultiBlockDataSet> ReadOBJBuildings(
    int numberOfBuildings,
    int vtkNotUsed( lod ),
    std::vector<std::string> const & files,
    std::array<double, 3>& fileOffset
)
{
    auto root = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    for ( size_t i = 0; i < files.size() && i < static_cast<size_t>( numberOfBuildings ); ++i )
    {
        vtkNew<vtkOBJReader> reader;
        reader->SetFileName( files[i].c_str() );
        reader->Update();
        if ( i == 0 )
        {
            fileOffset = ReadOBJOffset( reader->GetComment() );
        }
        auto polyData = reader->GetOutput();
        std::string textureFileName = GetOBJTextureFileName( files[i] );
        if ( !textureFileName.empty() )
        {
            SetField( polyData, "texture_uri", textureFileName.c_str() );
        }
        auto building = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        building->SetBlock( 0, polyData );
        root->SetBlock( root->GetNumberOfBlocks(), building );
    }
    return root;
}

TEST_F( GM13206_test, vtk_t00 )
{
    auto ws = create_workspace();

    using namespace vtksys;

    std::array<double, 3> fileOffset;
    vtkSmartPointer<vtkMultiBlockDataSet> actual
        = ReadOBJBuildings( 1, 42, { R"(V:\Asset Importer Test Models\models\OBJ\box.obj)"}, fileOffset);
    ASSERT_TRUE( actual );
    CONSOLE_TE( actual->GetNumberOfBlocks() );
}

/// @brief
///
/// https://gitlab.kitware.com/vtk/vtk/-/blob/v9.3.0/IO/Cesium3DTiles/Testing/Cxx/TestCesium3DTilesWriter.cxx
/// 
/// @param --gtest_filter=GM13206_test.DISABLED_vtk_t0 --gtest_also_run_disabled_tests
/// @param  
TEST_F( GM13206_test, DISABLED_vtk_t0 )
{
    auto ws = create_workspace();

    using namespace vtksys;

    std::string textureBaseDirectory = ws.string(); // TODO: set
    std::array<double, 3> fileOffset;

    const std::vector<std::string> input;
    int inputType = vtkCesium3DTilesWriter::Mesh;
    bool addColor{};
    std::string output = ws.string();
    bool contentGLTF = true;
    int numberOfBuildings{};
    int buildingsPerTile{};
    int lod{};
    std::vector<double> inputOffset;
    bool saveTiles = true;
    bool saveTextures = true;
    std::string crs;
    int utmZone = 19; // Maine
    char utmHemisphere = 'N';

    vtkSmartPointer<vtkMultiBlockDataSet> mbData
        = ReadOBJBuildings( 1, 42, { R"(V:\Asset Importer Test Models\models\OBJ\box.obj)" }, fileOffset );
    ASSERT_TRUE( mbData );
    EXPECT_EQ( 1, mbData->GetNumberOfBlocks() );

    vtkNew<vtkCesium3DTilesWriter> writer;

    writer->SetInputDataObject( mbData );

    writer->SetContentGLTF( contentGLTF );
    writer->ContentGLTFSaveGLBOff();
    writer->SetInputType( inputType );
    writer->SetDirectoryName( output.c_str() );
    writer->SetTextureBaseDirectory( textureBaseDirectory.c_str() );
    writer->SetOffset( fileOffset.data() );
    writer->SetSaveTextures( saveTextures );
    writer->SetNumberOfFeaturesPerTile( buildingsPerTile );
    writer->SetSaveTiles( saveTiles );
    if ( crs.empty() )
    {
        std::ostringstream ostr;
        ostr << "+proj=utm +zone=" << utmZone << ( utmHemisphere == 'S' ? "+south" : "" );
        crs = ostr.str();
    }
    writer->SetCRS( crs.c_str() );
    int rc = writer->Write();

    CONSOLE_T( "rc = " << rc );
}

#pragma endregion vtk

// TODO: ....

/// @brief Example of a single double format
/// @param --gtest_filter=GM13206_test.tileset_t0
/// @param  
TEST_F( GM13206_test, tileset_t00 )
{
    auto ws = create_workspace();

    json_object* actual =
        json_object_new_double( 22.0 / 7.0 );

    char* filename =
        strdup( ( ws / "actual_tileset.json" ).string().c_str() );
    EXPECT_EQ(0, json_object_to_file_ext( filename, actual, JSON_C_TO_STRING_PRETTY| JSON_C_TO_STRING_NOZERO ) );
    free( filename );
};

/// @brief Test TilesetJson in "create" mode
/// @param --gtest_filter=GM13206_test.tileset_t0
/// @param  
TEST_F( GM13206_test, tileset_t0 )
{
    auto ws = create_workspace();
    auto actual_tileset_json = ( ws / "actual_tileset.json" ).string();

    {
        TilesetJson actual;

        actual.asset( TilesetJson::asset_t{ "1.1" } );
        actual.geometricError( 4096 );
        auto root = TilesetJson::root_tile_t(
                         "ADD",
                         { 0.93932803788361, 0.3430201703190752, 0, 0,
                           -0.23782586639925077, 0.6512634642883413, 0.7206210913889448, 0,
                           0.24718756950375376, -0.6768995958319232, 0.6933291012538028, 0,
                           1579137.3428950259, -4324317.0816180855, 4399624.364052873, 1
                         } );

        root.boundingVolume.reset(
            new TilesetJson::boundingVolumeBox_t( {
                999.9999999999997,
                -500.02746582031233,
                250.018310546875,
                0,
                0,
                -250.018310546875,
                1050,
                0,
                0,
                0,
                -550.0274658203125,
                0
        } ) );
        root.children.push_back( TilesetJson::child_t{ 256, { "content.glb" } } );
        root.children.push_back( TilesetJson::child_t{ 256, { "content2.glb" } } );
        root.children.push_back( TilesetJson::child_t{ 256, { "content3.glb" } } );
        actual.root( root );
        EXPECT_TRUE( actual.save_as( actual_tileset_json.c_str() ) );
    }
    // 2. Verify

    {
        TilesetJson actual( actual_tileset_json.c_str() );

        EXPECT_DOUBLE_EQ( 4096, actual.geometricError() );

        auto actual_root = actual.root();
        EXPECT_STREQ( "ADD", actual_root.refine.c_str() );
        EXPECT_EQ( 3, actual_root.children.size() );
       
    }
}


#include "TestDataMgr.h"

/// @brief GM-17462 Export Support for Point Cloud to Cesium 3D Tiles
class GM17462_test : public GlobalMapperSDKTestF
{
protected:
    void SetUp() override
    {
        SetConsoleCP( CP_UTF8 );
        SetConsoleOutputCP( CP_UTF8 );

        close_all();
    }

    /// @brief Check verbosity settings
    /// @return true if env variable "SDK_TESTS_VERBOSE_ON" is defined 
    bool is_verbose() const
    {
        if ( getenv( "SDK_TESTS_VERBOSE_ON" ) )
            return true;
        // TODO: more conditions to check
        return false;
    }

protected:
    fs::path workspace_location() const
    {
        return fs::path( "out" ) / test_case_name() / test_name();
    }

    fs::path create_workspace()
    {
        auto ws = workspace_location();
        if ( fs::is_directory( ws ) )
            fs::remove_all( ws );
        std::error_code err;
        fs::create_directories( ws, err );

        CONSOLE( "ws = " << fs::absolute( ws ).string() );

        return ws;
    }

    auto metadata_list( GM_LayerHandle_t32 lp )
    {
        std::unordered_map<std::string, std::string> mm;

        auto info = GM_GetLayerInfo(lp);
        for ( int i = 0; i != info->mMetadataListSize; ++i )
            mm.insert( { info->mMetadataList[i].mName, info->mMetadataList[i].mVal } );

        return mm;
    }

    std::string save_as( std::string const& text, fs::path const& filename )
    {
        auto f = fs::absolute( filename );
        fs::create_directories( f.parent_path() );

        std::ofstream ff( filename.c_str() );
        ff << text;

        return f.string();
    }

    bool save_assets_as( aiScene const* scene, std::string base_name )
    {
        CONSOLE_TE( base_name );
        if ( !scene )
            return false;

        auto ws = workspace_location();
        if ( 0 )
        { // verbose?
            for ( int i = 0; i != scene->mNumMeshes; ++i )
            {
                CONSOLE_T( i << " : " << scene->mMeshes[i]->mName.C_Str() );
            }
        }

        {
            auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
            std::string f = ( ws / ( base_name + ".xml" ) ).string();
            Assimp::Exporter exp;
            if ( !aiReturn_SUCCESS, exp.Export( scene, "assxml", f, flags ) )
                return false;
        }
        {
            auto flags = aiProcess_GenNormals | aiProcess_ValidateDataStructure | 0;
            std::string f = ( ws / ( base_name + ".glb" ) ).string();
            Assimp::Exporter exp;
            if ( !aiReturn_SUCCESS, exp.Export( scene, "glb2", f, flags ) )
                return false;
        }
        return true;
    }

    bool save_assets_as_glb2( aiScene const* scene, std::string base_name )
    {
        CONSOLE_TE( base_name );
        if ( !scene )
            return false;

        auto ws = workspace_location();
        {
            auto flags = 0
                | aiProcess_GenNormals
                | aiProcess_ValidateDataStructure
                ;
            std::string f = ( ws / ( base_name + ".glb" ) ).string();
            Assimp::Exporter exp;
            if ( !aiReturn_SUCCESS, exp.Export( scene, "glb2", f, flags ) )
                return false;
        }
        return true;
    }

    bool save_text_as( std::string const& text, fs::path const &f )
    {
        std::filesystem::create_directories( f.parent_path() );

        std::ofstream ff( f );
        ff << text;

        return true;
    }
};

#include "../shared_include/GlobalMapperInterface_Lidar.h"

/// @brief 
/// @param --gtest_filter=GM17462_test.t0
/// @param  
TEST_F( GM17462_test, t0 )
{
    TestDataMgr tdb{
        R"(C:\home\work\GM-17462)",
        R"(Z:\GM-17462)",
    };

    auto ws = create_workspace();

    if (0) {
        auto lake_gmw = tdb.find_datafile( "lake.gmw" );
        if ( !fs::is_regular_file( lake_gmw ) )
            GTEST_SKIP() << tdb.getError();

        auto layers = load_layer_list( lake_gmw );
        ASSERT_EQ( 1, layers.size() );
        if ( is_verbose() )
        {
            for ( auto pp : metadata_list( layers.front() ) )
                CONSOLE( pp.first << " => " << pp.second );
        }
        CONSOLE( "" );
        close_all();
    }

    {
        auto filename = tdb.find_datafile( "mike_niwot_Generated_Point_Cloud.gmw" );
        if ( !fs::is_regular_file( filename ) )
            GTEST_SKIP() << tdb.getError();
        auto layers = load_layer_list( filename );
        ASSERT_EQ( 1, layers.size() );
        if ( is_verbose() )
        {
            for ( auto pp : metadata_list( layers[0] ) )
                CONSOLE( pp.first << " => " << pp.second );
            CONSOLE( "" );
        }

        {
            CONSOLE("*** STATISTICS ***");

            GM_LidarStats_t lidarStats;
            auto rc = GM_GetLayerLidarStats( layers[0], &lidarStats, 0 );
            ASSERT_EQ( 0, rc );
            CONSOLE_EVAL( lidarStats.mAllPointsStats.mCount );
            EXPECT_EQ( 1001824, lidarStats.mAllPointsStats.mCount );


            std::vector<meshtoolbox::lidarptr_t> points( lidarStats.mAllPointsStats.mCount );

            for ( uint64 i = 0; i != lidarStats.mAllPointsStats.mCount; ++i )
            {
                GM_LidarPoint_t lidarPoint;
                GM_GetFeatureFlags_t32 flags = 0;
                auto rc = GM_GetLidarPoint( layers[0], i, &lidarPoint, flags, 0 );
                ASSERT_EQ( 0, rc );
                points[i] = {
                    ai_real( lidarPoint.mPos.mX ), ai_real( lidarPoint.mPos.mY ), ai_real( lidarPoint.mElevMeters ),
                    lidarPoint.mRed, lidarPoint.mGreen, lidarPoint.mBlue
                };
            }
            meshtoolbox::Toolbox tb;

            auto actual = std::unique_ptr<aiScene>( tb.make_lidar_pc( points.data(), points.size() ) );
            Assimp::Exporter expo;
            {
                aiReturn r = expo.Export( actual.get(), "gltf2", ( ws / "export-pc-lidar.gltf" ).string() );
                EXPECT_EQ( 0, r );
            }
            {
                aiReturn r = expo.Export( actual.get(), "glb2", ( ws / "export-pc-lidar.glb" ).string() );
                EXPECT_EQ( 0, r );
            }
            {
                aiReturn r = expo.Export( actual.get(), "assxml", ( ws / "export-pc-lidar.xml" ).string() );
                EXPECT_EQ( 0, r );
            }
        }
        close_all();
    }
}

/// @brief Export the lake point cloud using the SDK functionality 
/// @param --gtest_filter=GM17462_test.t0_lake_as_pdf
/// @param  
TEST_F( GM17462_test, t0_lake_as_pdf )
{
    TestDataMgr tdb{
        R"(C:\home\work\GM-17462)",
        R"(Z:\GM-17462)",
    };

    auto ws = create_workspace();

    {
        auto lake_gmw = tdb.find_datafile( "lake.gmw" );
        if ( !fs::is_regular_file( lake_gmw ) )
            GTEST_SKIP() << tdb.getError();

        auto layers = load_layer_list( lake_gmw );
        ASSERT_EQ( 1, layers.size() );
        if ( is_verbose() )
        {
            for ( auto pp : metadata_list( layers.front() ) )
                CONSOLE( pp.first << " => " << pp.second );
            CONSOLE( "" );
        }

        {
            const char export_cmd[] = R"(// A comment 
            SET_LOG_FILE FILENAME="lake.pdf.log" APPEND_TO_FILE=NO
            EXPORT_PDF FILENAME="lake.pdf" PDF_PAGE_SIZE="Letter" PDF_PAGE_ORIENTATION="PORTRAIT" \
                 PDF_MARGINS="0.750000,0.750000,0.750000,0.750000" PDF_HEADER="" PDF_FOOTER="" PDF_COMBINE_RASTERS="YES" \
                 PDF_COMBINE_ELEV_GRIDS="YES" PDF_CREATE_IMAGE_DRAPED="YES" PDF_LIDAR_AS_LINES="YES" \
                 PDF_USE_FILENAME_FOR_NODES="NO" PDF_USE_LDP="YES" PDF_USE_ADOBE_EXT="NO"
            )";

            auto export_cmd_file = ( ws / "export_cmd_file.gms" ).string();

            (void)save_as( export_cmd, export_cmd_file );

            GM_LayerHandle_t32* ignore_list = 0;
            uint32 ignore_list_size = 0;
            GM_LoadFlags_t32 flags = GM_LoadFlags_HideAllPrompts;

            EXPECT_EQ( 0, GM_RunScript( export_cmd_file.c_str(), &ignore_list, &ignore_list_size, flags, 0 ) );
            EXPECT_TRUE( fs::is_regular_file( ws / "lake.pdf.log" ) );
            EXPECT_TRUE( fs::is_regular_file( ws / "lake.pdf" ) );
        }
        close_all();
    }
}

/// @brief Export the lake point cloud using the SDK functionality 
/// @param --gtest_filter=GM17462_test.t0_lake_as_glb
/// @param  
TEST_F( GM17462_test, t0_lake_as_glb )
{
    TestDataMgr tdb{
        R"(C:\home\work\GM-17462)",
        R"(Z:\GM-17462)",
    };

    auto ws = create_workspace();

    {
        auto lake_gmw = tdb.find_datafile( "lake.gmw" );
        if ( !fs::is_regular_file( lake_gmw ) )
            GTEST_SKIP() << tdb.getError();

        auto layers = load_layer_list( lake_gmw );
        ASSERT_EQ( 1, layers.size() );
        if ( is_verbose() )
        {
            for ( auto pp : metadata_list( layers.front() ) )
                CONSOLE( pp.first << " => " << pp.second );
            CONSOLE( "" );
        }

        {
            const char export_cmd[] = R"(// A comment 
            SET_LOG_FILE FILENAME="lake.glb.log" APPEND_TO_FILE=NO
            EXPORT_VECTOR EXPORT_LAYER="lake.laz" FILENAME="lake.glb" TYPE="GLB" \
                Y_UP = "NO" CREATE_BINARY = "NO" GEN_PRJ_FILE = "NO"
            )";

            auto export_cmd_file = ( ws / "export_cmd_file.gms" ).string();

            (void)save_as( export_cmd, export_cmd_file );

            GM_LayerHandle_t32* ignore_list = 0;
            uint32 ignore_list_size = 0;
            GM_LoadFlags_t32 flags = GM_LoadFlags_HideAllPrompts;

            EXPECT_EQ( 0, GM_RunScript( export_cmd_file.c_str(), &ignore_list, &ignore_list_size, flags, 0 ) );
            EXPECT_TRUE( fs::is_regular_file( ws / "lake.glb.log" ) );
            EXPECT_TRUE( fs::is_regular_file( ws / "lake.glb" ) );
        }
        close_all();
    }
}

/// @brief Export a workspace with two point cloud layers using the SDK functionality 
/// @param --gtest_filter=GM17462_test.t0_lake_as_glb
/// @param  
TEST_F( GM17462_test, t0_two_lakes_as_glb )
{
    TestDataMgr tdb{
        R"(C:\home\work\GM-17462)",
        R"(Z:\GM-17462)",
    };

    auto ws = create_workspace();

    {
        auto lake_gmw = tdb.find_datafile( "two_lakes.gmw" );
        if ( !fs::is_regular_file( lake_gmw ) )
            GTEST_SKIP() << tdb.getError();

        auto layers = load_layer_list( lake_gmw );
        ASSERT_EQ( 2, layers.size() );
        if ( is_verbose() )
        {
            for ( auto pp : metadata_list( layers.front() ) )
                CONSOLE( pp.first << " => " << pp.second );
            CONSOLE( "" );
        }

        {
            const char export_cmd[] = R"(// A comment - will be ignored
SET_LOG_FILE FILENAME="two_lakes.glb.log" APPEND_TO_FILE=NO
EXPORT_VECTOR EXPORT_LAYER="lake.laz" EXPORT_LAYER="lake.laz [COPY]" FILENAME="two_lakes.glb" \
    TYPE="GLB" Y_UP="NO" CREATE_BINARY="NO" GEN_PRJ_FILE="NO")";

            auto export_cmd_file = ( ws / "export_cmd_file.gms" ).string();

            (void)save_as( export_cmd, export_cmd_file );

            GM_LayerHandle_t32* ignore_list = 0;
            uint32 ignore_list_size = 0;
            GM_LoadFlags_t32 flags = GM_LoadFlags_HideAllPrompts;

            EXPECT_EQ( 0, GM_RunScript( export_cmd_file.c_str(), &ignore_list, &ignore_list_size, flags, 0 ) );
            EXPECT_TRUE( fs::is_regular_file( ws / "two_lakes.glb.log" ) );
            EXPECT_TRUE( fs::is_regular_file( ws / "two_lakes.glb" ) );
        }
        close_all();
    }
}


/// @brief Create PC glTF output
/// @param --gtest_filter=GM17462_test.ai_max
TEST_F( GM17462_test, ai_max )
{
    CONSOLE( "AI_MAX_FACE_INDICES: " << AI_MAX_FACE_INDICES );
    CONSOLE( "AI_MAX_BONE_WEIGHTS: " << AI_MAX_BONE_WEIGHTS );
    CONSOLE( "AI_MAX_VERTICES: " << AI_MAX_VERTICES );
    CONSOLE( "AI_MAX_FACES: " << AI_MAX_FACES );
    CONSOLE( "AI_MAX_NUMBER_OF_COLOR_SETS: " << AI_MAX_NUMBER_OF_COLOR_SETS );
}

/// @brief Create PC glTF output ("box")
/// @param --gtest_filter=GM17462_test.meshtoolbox_0
TEST_F( GM17462_test, meshtoolbox_0 )
{
    meshtoolbox::Toolbox tb;

    {
        auto actual = std::unique_ptr<aiMesh>( tb.mesh_pc_box( { 0,0,0 }, { 1,1,1 }, "TheBox" ) );
        ASSERT_TRUE( actual );

        EXPECT_STREQ( "TheBox", actual->mName.C_Str() );
        EXPECT_TRUE( actual->HasFaces() );
        EXPECT_TRUE( actual->HasPositions() );
        EXPECT_FALSE( actual->HasNormals() );
        EXPECT_FALSE( actual->HasBones() );

        EXPECT_EQ( 8, actual->mNumVertices );
        EXPECT_EQ( 8, actual->mNumFaces );
#if 0
        auto get_face_positions = []( aiMesh const* mp, unsigned k )->std::tuple<aiVector3D, aiVector3D, aiVector3D>
            {
                if ( !( k >= 0 && k < mp->mNumFaces ) )
                    throw std::runtime_error( "Failed: k >= 0 && k < mp->mNumFaces" );
                if ( !( mp->mFaces[k].mNumIndices == 3 ) )
                    throw std::runtime_error( "Failed: mp->mFaces[k].mNumIndices == 3" );
                return { mp->mVertices[mp->mFaces[k].mIndices[0]], mp->mVertices[mp->mFaces[k].mIndices[1]],mp->mVertices[mp->mFaces[k].mIndices[2]] };
            };

        /* Compute normals to the faces */
        for ( unsigned i = 0; i < actual->mNumFaces; ++i )
        {
            auto [v0, v1, v2] = get_face_positions( actual.get(), i );

            auto normal = tb.from_eigen( tb.to_eigen( v1 - v0 ).cross( tb.to_eigen( v2 - v0 ) ) );
            CONSOLE_T( i << ": " << "v0: " << v0 << ", v1: " << v1 << ", v2: " << v2
                     << ", normal: " << normal << ", plane: " << tb.plane_normal( v0, v1, v2 ) );
        }
#endif
    }
    //{
    //    auto actual = std::unique_ptr<aiMesh>( tb.mesh_box( { 0,0,0 }, { 1,1,1 } ) );
    //    ASSERT_TRUE( actual );
    //    EXPECT_STREQ( "", actual->mName.C_Str() );
    //    EXPECT_EQ( 8, actual->mNumVertices );
    //}

    {
        auto ws = create_workspace();

        auto actual = std::unique_ptr<aiScene>( tb.make_pc_box( { 0,0,0 }, { 1,1,1 } ) );
        ASSERT_TRUE( actual );

        EXPECT_TRUE( actual->HasMeshes() );
        {
            auto mp = actual->mMeshes[0];
            ASSERT_TRUE( mp );
        }

        Assimp::Exporter expo;
        {
            aiReturn r = expo.Export( actual.get(), "gltf2", ( ws / "export-pc-box.gltf" ).string() );
            EXPECT_EQ( 0, r );
        }
        {
            aiReturn r = expo.Export( actual.get(), "assxml", ( ws / "export-pc-box.xml" ).string() );
            EXPECT_EQ( 0, r );
        }
    }
}

/// @brief Generate bounding boxes
/// @param --gtest_filter=GM17462_test.bounding_boxes
/// @param  
TEST_F( GM17462_test, bounding_boxes )
{
    if ( !is_explicitly_called() )
        GTEST_SKIP() << "is_explicitly_called";

    auto ws = create_workspace();

    TestDataMgr tdb{
        R"(C:\home\work\GM-17462)",
        R"(Z:\GM-17462)",
    };

    std::function<void( TilesetJson&, TilesetJson::child_t const&, aiScene* )> node_appender;
    node_appender = [&node_appender]( TilesetJson& actual, TilesetJson::child_t const& root_tile, aiScene* scene )
        {
            for ( auto cp = root_tile.children.begin(); cp != root_tile.children.end(); ++cp )
            {
                actual.append_node( scene, cp->boundingVolume.get() );
                node_appender( actual, *cp, scene );
            }
        };

    {
        auto tileset_filename = tdb.find_datafile( "mike_niwot_Generated_Point_Cloud/tileset.json" );
        if ( tileset_filename.empty() )
            GTEST_SKIP() << tdb.getError();

        TilesetJson actual( tileset_filename );

        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 1.12152826801556, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "mike_niwot_Generated_Point_Cloud" ) );
    }
    {
        // "C:\home\work\GM-17462\lake_no_draco\tileset.json"
        auto tileset_filename = tdb.find_datafile( R"(lake_no_draco\tileset.json)" );
        if ( tileset_filename.empty() )
            GTEST_SKIP() << tdb.getError();

        TilesetJson actual( tileset_filename );

        auto root_tile = actual.root();
        EXPECT_STREQ( "ADD", root_tile.refine.c_str() );
        EXPECT_DOUBLE_EQ( 4.9856285997004699, root_tile.geometricError );
        EXPECT_TRUE( (bool)root_tile.boundingVolume );

        auto root_volume = root_tile.boundingVolume.get();
        std::unique_ptr<aiScene> scene{ actual.construct_scene( root_volume ) };

        node_appender( actual, root_tile, scene.get() );

        EXPECT_TRUE( save_assets_as( scene.get(), "lake_no_draco" ) );
    }
}

/// @brief Point Cloud file (.pnts)
/// @see-also https://github.com/CesiumGS/3d-tiles/tree/main/specification/TileFormats/PointCloud
/// 
class PointCloud
{
public:
    PointCloud()
        : feature_table{}
    {};

    bool read( std::string const& filename );

    bool false_because( std::string s )
    {
        CONSOLE( "Error: " << s );
        return false;
    }

    bool parse_feature_table();
    bool is_draco_enabled() const
    {
        return extensions.count( "3DTILES_draco_point_compression" );
    }
public:
    struct header_t
    {
        char magic[4];
        uint32_t version;
        uint32_t byteLength;
        uint32_t featureTableJSONByteLength;
        uint32_t featureTableBinaryByteLength;
        uint32_t batchTableJSONByteLength;
        uint32_t batchTableBinaryByteLength;
    };

    typedef float triplet_t[3];

    struct feature_table_t
    {
        triplet_t *POSITION;
        uint16_t( *POSITION_QUANTIZED )[3];
        uint32_t POINTS_LENGTH;
        float RTC_CENTER[3];
    };

    std::string featureTableJSON;
    std::vector<char> featureTableBinary;
    std::string batchTableJSON;
    std::vector<char> batchTableBinary;
    std::vector<char> binary_glTF;
    feature_table_t feature_table;
    std::unordered_map<std::string, std::string> extensions;
};

bool PointCloud::read( std::string const& filename )
{
    static_assert( 28 == sizeof( header_t ), "PointCloud::header_t must be 28 bytes long" );

    CONSOLE_EVAL( filename );

    std::ifstream pnts( filename, std::ios::binary );

    if ( pnts.bad() )
        return false_because( "Cannot open " + filename );

    header_t hdr{};
    pnts.read( (char*)&hdr, sizeof( hdr ) );
    if ( !( hdr.magic[0] == 'p' && hdr.magic[1] == 'n' && hdr.magic[2] == 't' && hdr.magic[3] == 's' ) )
        return false_because( "Wrong " + filename );
    CONSOLE_EVAL( hdr.byteLength );
    CONSOLE_EVAL( hdr.version );
    CONSOLE_EVAL( hdr.featureTableJSONByteLength );
    CONSOLE_EVAL( hdr.featureTableBinaryByteLength );
    CONSOLE_EVAL( hdr.batchTableJSONByteLength );
    CONSOLE_EVAL( hdr.batchTableBinaryByteLength );

    if ( hdr.version != 1 )
        return false_because( "Only Point Cloud version 1 is supported.  Version "
                              + std::to_string( hdr.version )
                              + " is not." );

    // Read the JSON
    if ( auto jsob_buffer_size = hdr.featureTableJSONByteLength )
    {
        char* jsob_buffer = new char[jsob_buffer_size];
        pnts.read( jsob_buffer, jsob_buffer_size );
        featureTableJSON = std::string( jsob_buffer, jsob_buffer_size );
        delete[] jsob_buffer;

        CONSOLE_EVAL( featureTableJSON );
    }

    if ( auto length = hdr.featureTableBinaryByteLength )
    {
        featureTableBinary.resize(length);
        pnts.read( featureTableBinary.data(), length);
    }

    if ( auto length = hdr.batchTableJSONByteLength )
    {
        char* buffer = new char[length];
        pnts.read( buffer, length );
        batchTableJSON = std::string( buffer, length );
        delete[] buffer;

        CONSOLE_EVAL( batchTableJSON );
    }

    if ( auto length = hdr.batchTableBinaryByteLength )
    {
        batchTableBinary.resize(length);
        pnts.read( batchTableBinary.data(), length);
    }

    auto pos = pnts.tellg();
    CONSOLE_EVAL( pos );

    {
        size_t length = hdr.byteLength - pos;
        if ( length )
        {
            binary_glTF.resize( length );
            pnts.read( binary_glTF.data(), length );
        }
    }

    return true;
}

bool PointCloud::parse_feature_table()
{
    /*
    * No draco
      {
        "POINTS_LENGTH":57363,
        "POSITION":{"byteOffset":0},
        "RTC_CENTER":[477074.955,4366597.995,2747.0150000000003]
      }
    */
    /*
    * With draco
        {
            "POINTS_LENGTH":57363,
            "POSITION":{"byteOffset":0},
            "RTC_CENTER":[477074.955,4366597.995,2747.0150000000003],
            "extensions":{
                "3DTILES_draco_point_compression":{
                    "byteLength":780259,
                    "byteOffset":0,
                    "properties":{"POSITION":0}
                }
            }
        }
    */
    memset( &feature_table, 0, sizeof( feature_table ) );

    auto root = json_tokener_parse( featureTableJSON.c_str() );
    if ( !root )
        return false_because( "Cannot parse featureTableJSON" );

    std::unique_ptr<json_object, std::function<void( json_object* )>>
        on_return( root, []( json_object* r ) { json_object_put( r ); CONSOLE( "Goodbye cruel world..." ); } );

    {
        auto o = json_object_object_get( root, "POINTS_LENGTH" );
        if ( !o )
            return false_because( "POINTS_LENGTH is not defined" );
        feature_table.POINTS_LENGTH = json_object_get_int( o );
    }

    if ( auto eo = json_object_object_get( root, "extensions" ) )
    {
        if ( auto e_draco = json_object_object_get( eo, "3DTILES_draco_point_compression" ) )
        {
            // TODO: parse the draco extensions parameters
            this->extensions.insert( { "3DTILES_draco_point_compression", "1" } );
        }
    }

    if ( auto po = json_object_object_get( root, "POSITION" ) )
    {
        unsigned byteOffset = 0;
        if ( auto bo = json_object_object_get( root, "byteOffset" ) )
            byteOffset = json_object_get_int( bo );
        feature_table.POSITION = (triplet_t*)( featureTableBinary.data() + byteOffset );
    }
    if ( auto o = json_object_object_get( root, "RTC_CENTER" ) )
    {
        if ( json_type_array != json_object_get_type( o ) )
            return false_because( "RTC_CENTER must be array" );
        array_list* arr = json_object_get_array( o );
        if (arr->length != 3 )
            return false_because( "RTC_CENTER must be array of 3" );
        for ( int i = 0; i < arr->length; ++i )
            feature_table.RTC_CENTER[i] = json_object_get_double( (json_object*)array_list_get_idx( arr, i ) );
    }
    if ( feature_table.POSITION == nullptr && feature_table.POSITION_QUANTIZED == nullptr )
        return false_because( "Nether POSITION nor POSITION_QUANTIZED defined" );

    return true;
}

/// @brief Read 0/0.pnts tile
/// @param --gtest_filter=GM17462_test.mike_niwot_Generated_Point_Cloud_0_0
/// @param  
TEST_F( GM17462_test, mike_niwot_Generated_Point_Cloud_0_0 )
{
    if ( !is_explicitly_called_suite() )
        GTEST_SKIP() << "is_explicitly_called";

    auto ws = create_workspace();

    {
        TestDataMgr tdb{
            R"(C:\home\work\GM-17462)",
            R"(Z:\GM-17462)",
        };

        auto model_filename = tdb.find_datafile( R"(mike_niwot_Generated_Point_Cloud\0\0.pnts)" );
        if ( model_filename.empty() )
            GTEST_SKIP() << tdb.getError();

        PointCloud actual;

        ASSERT_TRUE( actual.read( model_filename ) )
            << "Failed: " << model_filename;
        EXPECT_EQ( 0, actual.binary_glTF.size() );
        EXPECT_EQ( 276, actual.featureTableJSON.size() );
        EXPECT_EQ( 381360, actual.featureTableBinary.size() );
        EXPECT_EQ( 368, actual.batchTableJSON.size() );
        EXPECT_EQ( 0, actual.batchTableBinary.size() );

        save_text_as( actual.featureTableJSON, ws / "featureTable.json" );
        save_text_as( actual.batchTableJSON, ws / "batchTable.json" );

        if ( actual.binary_glTF.size() )
        {
            Assimp::Importer imp;
            auto model = imp.ReadFileFromMemory( actual.binary_glTF.data(), actual.binary_glTF.size(), 0 );
            ASSERT_TRUE( model != nullptr );
            save_assets_as( model, "tile_0_0" );
        }
    }
}

/// @brief Read 0/0.pnts tile w/o Draco compression
///
/// E.g. "C:\home\work\GM-17462\lake_no_draco\0\0.pnts"
/// @param --gtest_filter=GM17462_test.lake_no_draco
/// @param  
TEST_F( GM17462_test, lake_no_draco )
{
    if ( !is_explicitly_called_suite() )
        GTEST_SKIP() << "is_explicitly_called";

    auto ws = create_workspace();

    {
        TestDataMgr tdb{
            R"(C:\home\work\GM-17462)",
            R"(Z:\GM-17462)",
        };

        auto model_filename = tdb.find_datafile( R"(lake_no_draco\0\0.pnts)" );
        if ( model_filename.empty() )
            GTEST_SKIP() << tdb.getError();

        PointCloud actual;

        ASSERT_TRUE( actual.read( model_filename ) )
            << "Failed: " << model_filename;
        EXPECT_EQ( 0, actual.binary_glTF.size() );
        EXPECT_EQ( 108, actual.featureTableJSON.size() );
        EXPECT_EQ( 688360, actual.featureTableBinary.size() );
        EXPECT_EQ( 504, actual.batchTableJSON.size() );
        EXPECT_EQ( 860448, actual.batchTableBinary.size() );

        save_text_as( actual.featureTableJSON, ws / "featureTable.json" );
        save_text_as( actual.batchTableJSON, ws / "batchTable.json" );

        EXPECT_TRUE( actual.parse_feature_table() );
        EXPECT_EQ( 57363, actual.feature_table.POINTS_LENGTH );
        EXPECT_FLOAT_EQ( 477074.97, actual.feature_table.RTC_CENTER[0] );
        EXPECT_FLOAT_EQ( 4366598, actual.feature_table.RTC_CENTER[1] );
        EXPECT_FLOAT_EQ( 2747.0149, actual.feature_table.RTC_CENTER[2] );

        {
            std::ofstream pc( ws / "point_cloud.csv" );
            for ( int i = 0; i != actual.feature_table.POINTS_LENGTH; ++i )
            {
                float* f3 = actual.feature_table.POSITION[i];
                pc << f3[0] << "," << f3[1] << "," << f3[2] << std::endl;
            }
            CONSOLE( "Saved to " << ( ws / "point_cloud.csv" ).string() );
        }
    }
}

/// @brief Read 0/0.pnts tile with Draco compression
///
/// E.g. "C:\home\work\GM-17462\lake_no_draco\0\0.pnts"
/// @param --gtest_filter=GM17462_test.lake
/// @param  
TEST_F( GM17462_test, lake )
{
    if ( !is_explicitly_called_suite() )
        GTEST_SKIP() << "is_explicitly_called";

    auto ws = create_workspace();

    {
        TestDataMgr tdb{
            R"(C:\home\work\GM-17462)",
            R"(Z:\GM-17462)",
        };

        auto model_filename = tdb.find_datafile( R"(lake\0\0.pnts)" );
        if ( model_filename.empty() )
            GTEST_SKIP() << tdb.getError();

        PointCloud actual;

        ASSERT_TRUE( actual.read( model_filename ) )
            << "Failed: " << model_filename;
        EXPECT_EQ( 0, actual.binary_glTF.size() );
        EXPECT_EQ( 228, actual.featureTableJSON.size() );
        EXPECT_EQ( 780264, actual.featureTableBinary.size() );
        EXPECT_EQ( 640, actual.batchTableJSON.size() );
        EXPECT_EQ( 0, actual.batchTableBinary.size() );

        save_text_as( actual.featureTableJSON, ws / "featureTable.json" );
        save_text_as( actual.batchTableJSON, ws / "batchTable.json" );

        EXPECT_TRUE( actual.parse_feature_table() );
        EXPECT_EQ( 57363, actual.feature_table.POINTS_LENGTH );
        EXPECT_FLOAT_EQ( 477074.97, actual.feature_table.RTC_CENTER[0] );
        EXPECT_FLOAT_EQ( 4366598, actual.feature_table.RTC_CENTER[1] );
        EXPECT_FLOAT_EQ( 2747.0149, actual.feature_table.RTC_CENTER[2] );
    }
}
