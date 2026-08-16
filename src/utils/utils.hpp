// Utils.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <glm/vec3.hpp>
#include <vector>

// MAJOR TODO: For any instance of a std::vector<Eigen::Vector2d>, you MUST add an Eigen allocator to prevent bad behavior when 
// compiling with c++11 or c++14. See documentation: https://libeigen.gitlab.io/eigen/docs-nightly/group__TopicStlContainers.html
// While default compiler is c++17, this prevents annoyances if other people want to work on the code with a different c++.

namespace Utils {

    // GLM::vec3 to Eigen::Vector3d converter
    Eigen::Vector3d glmToEigen(const glm::vec3 input);

    // Eigen::Vector3d to GLM::vec3 converter
    glm::vec3 eigenToGLM(const Eigen::Vector3d input);

    // Convert an eigen matrix with 3 columns to a std::vector of vector3d's
    void EigM3toStdV(const Eigen::MatrixXd& mat, std::vector<Eigen::Vector3d>& vec);

    // Entire mesh conversion routine Eigen to GLM
    void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM);

    // Entire mesh conversion routine GLM to Eigen
    // Note that GLM is float, while Eigen prefers double
    // Do NOT convert back and forth, you will lose information
    void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM);

    // Copy positions and connectivity into a copied container
    void copyPositions(const std::vector<Eigen::Vector3d>& V_old, std::vector<Eigen::Vector3d>& V_new);

    void copyConnectivity(const std::vector<std::vector<int>>& T_old, std::vector<std::vector<int>>& T_new);

    // SORTING
    // Sort a reference list of doubles while sorting their indices in the same way
    void doubleListIdxSort(std::vector<double>& ref_List, std::vector<int>& idx_List);

    // Insert at index between a pair of indices in a list
    bool insertIdxBetweenPair(std::vector<int>& idxList, int a, int b, int new_idx);

    // Flattens an Eigen::Matrix3d into a 9x1 row vector
    // NOTE: Does so column-wise!
    Eigen::VectorXd flattenMatrix3d(const Eigen::Matrix3d& F);
    // Compresses a 9x1 Eigen::VectorXd into an Eigen::Matrix3d
    // Assumes column-wise storage
    Eigen::Matrix3d compressVector9d(const Eigen::VectorXd& f);

    // MESH INIT HELPERS
    std::pair<int, int> undirectedKey(int a, int b);
    int edgeDirRelativeToKey(int a, int b);
    bool orientFacesConsistently(std::vector<std::vector<int>>& F_List);

    // VECTOR/PROJECTION HELPERS
    // Rotation variant SVD
    void rotationVariantSVD(Eigen::Matrix3d& mat, Eigen::Matrix3d& U, Eigen::Vector3d& Sigma, Eigen::Matrix3d& V);
    void polarDecomposition(Eigen::Matrix3d& mat, Eigen::Matrix3d& R, Eigen::Matrix3d& S);

    // Compute signed angle between two vectors in 3D given the axis
    double signedAngle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& axis, bool positive = false);
    
    // Computes the closest point to a triangle
    Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d> triVerts, const Eigen::Vector3d p);
    // Overload with projection type
    Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d>& triVerts, const Eigen::Vector3d& p, 
                                         int& projType, int& projIdx, double snapTol = 0.0);

    // Computes the closest point to a bilinear patch
    Eigen::Vector3d bilinearPatchClosestPoint(const std::vector<Eigen::Vector3d>& patchVerts, const Eigen::Vector3d& p, double eps = 1e-6, int max_iter = 15);
    int bilinearPatchClosestPoint(const std::vector<Eigen::Vector3d>& patchVerts, const Eigen::Vector3d& p, double& u, double& v, double eps = 1e-6, int max_iter = 15);
    // Overload with projection type
    Eigen::Vector3d bilinearPatchClosestPoint(const std::vector<Eigen::Vector3d>& patchVerts, const Eigen::Vector3d& p,
                                              int& projType, int& projIdx, double snapTol = 0.0, double eps = 1e-8, int max_iter = 15);
    // Evaluate a bilinear patch query point given u and v
    Eigen::Vector3d bilinearPatch(const std::vector<Eigen::Vector3d>& patchVerts, double u, double v);
    
    // Computes closest point to a Newell polygon
    Eigen::Vector3d polygonClosestPointNewell(const std::vector<Eigen::Vector3d>& polyVerts, const Eigen::Vector3d& p,
                                              const Eigen::Vector3d& polyNormal, int& projType, int& projIdx, double snapTol = 0.0);


    // Compute the angle of a vector starting at an origin and ending at a target, when projected onto a basis spanned by t1 and t2.
    bool directionAngleInPlane(
        const Eigen::Vector3d& origin,
        const Eigen::Vector3d& target,
        const Eigen::Vector3d& normal,
        const Eigen::Vector3d& t1,
        const Eigen::Vector3d& t2,
        double& theta
    );

    // Given two unit vectors, compute the rotation from one to the other
    Eigen::Matrix3d computeRotation(const Eigen::Vector3d& u, const Eigen::Vector3d& v);
    // Overload given axis and rotation
    Eigen::Matrix3d computeRotation(const Eigen::Vector3d& axis, const double& theta);


    // Returns true if two angular values are effectively the same direction.
    bool anglesCoincident(double a, double b, double eps = 1e-10);

    // Project a vector onto a tangent plane, given the normal to the plane
    // Returns -1 if degenerate (shouldn't happen but we should handle it)
    double projectVectorOntoTangentPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& vec, Eigen::Vector3d& proj, double scale = 1.0);

    // Projects a point onto the tangent plane of a normal given a center 
    Eigen::Vector3d projectPointOntoPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& center, const Eigen::Vector3d& p);

    // Build basis (t1, t2) for a plane given a normal
    void buildPlaneBasis(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);

    // Given a point on a plane basis and the plane basis, convert to 2D planar point
    Eigen::Vector2d convertTo2D(const Eigen::Vector3d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2);

    // Given a 2D planar point and the plane basis, revert to its 3D counterapart
    Eigen::Vector3d revertTo3D(const Eigen::Vector2d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2);

    // Check if a 2D point is in a 2D polygon
    // To do this, we do raycasting to the segment
    bool pointInPolygon2D(const Eigen::Vector2d& p, const std::vector<Eigen::Vector2d>& poly);

    // Returns the magnitude of the "cross product" of two 2D vectors.
    double cross2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b);
    // Find the intersection between a ray and a segment in 2D, if one exists.
    bool raycastToSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& direc, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, double& t, double& u);

    // Get the closest point on a segment in 2D and 3D, where the endpoints are defined
    // To do this, project onto parameterized segment and snap t to [0, 1]
    // TODO: Can we combine the 2D and 3D cases using VectorXd?
    Eigen::Vector2d closestPointOnSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, bool clip = true);

    Eigen::Vector3d closestPointOnSegment3D(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, bool clip = true);

    // MEAN VALUE COORDINATES
    // Helpers
    double vectorAngle(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3);
    double computeSign(const double& value);

    // Returns the weights only
    // Follows method of Fuda and Hormann [2024]
    void meanValueCoordinates(const Eigen::Vector2d& target, const std::vector<Eigen::Vector2d>& cage, Eigen::VectorXd& weights);

} // namespace Utils