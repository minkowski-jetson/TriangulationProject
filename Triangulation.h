#ifndef Triangulation_H
#define Triangulation_H

#include <vector>
#include <iostream>
#include <stdio.h>
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/QR"
#include "Eigen/SparseCore"
#include "Eigen/SVD"
#include "Eigen/StdVector"


typedef Eigen::Vector2f Vec2f;
typedef Eigen::Vector3f Vec3f;
typedef Eigen::Vector4f Vec4f;
typedef Eigen::VectorXf Vecf;


typedef Eigen::Matrix<float, 3, 4> Mat34f;
typedef Eigen::MatrixXf Matf;
typedef Eigen::Matrix4f Mat4f;
typedef Eigen::Matrix3f Mat3f;

using namespace std;

namespace UTILS
{

    // ** added section ** //
    void test_print();

    void test_func(const Mat3f& k,const Mat3f& r);

    // ** END added section ** //



    template <typename TMat, typename TVec>
    float nullspace(TMat *A, TVec *x);

    void residule(const Vec3f& pt3d, const std::vector<Vec2f>& pts, const std::vector<Mat34f>& Q, Vecf& e, Matf& J);

    class Triangulation
    {

    public:
        Triangulation(/* args */);
        ~Triangulation();

        void Linear(const std::vector<Mat34f>& poses, const std::vector<Vec2f>& pts, Vec3f& pt3d);
        void Sampson(const std::vector<Mat34f>& poses, const std::vector<Vec2f>& pts, Vec3f& pt_triangulated);
    };
    
}

#endif // start of file closing
