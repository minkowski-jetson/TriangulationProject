
#include "Triangulation.h"

using namespace std;


// namespace UTILS{

    // void test_print()
    // {
    //     cout << "test print" << endl;
    // }

    // void test_func(const Mat3f& k,const Mat3f& r)
    // {
    //     std::cout << "K" << std::endl;
    //     std::cout << k.format(HeavyFmt) << std::endl;
    //     std::cout << "R" << std::endl;
    //     std::cout << r.format(HeavyFmt) << std::endl;

    // }
// }

UTILS::Triangulation::Triangulation(/* args */)
{
    std::cout << "triangulation object instantiated" << std::endl;
};

UTILS::Triangulation::~Triangulation()
{
    std::cout << "triangulation destructor called" << std::endl;
};

template <typename TMat, typename TVec>
float UTILS::nullspace(TMat *A, TVec *x)
{
    if ( A->rows() >= A->cols() )
    {
        Eigen::JacobiSVD<TMat> svd( *A, Eigen::ComputeFullV );
        ( *x ) = svd.matrixV().col( A->cols() - 1 );
        return svd.singularValues()( A->cols() - 1 );
    }
    // More columns than rows, extend A with rows of zeros to make it square.
    TMat A_extended( A->cols(), A->cols() );
    A_extended.block( A->rows(), 0, A->cols() - A->rows(), A->cols() ).setZero();
    A_extended.block( 0, 0, A->rows(), A->cols() ) = ( *A );
    return nullspace(&A_extended, x);
};

void UTILS::residule(const Vec3f& pt3d, const std::vector<Vec2f>& pts, const std::vector<Mat34f>& Q, Vecf& e, Matf& J)
{
    int sz = static_cast<int>(pts.size());
    Vec3f x;

    for( int i=0; i<sz; ++i )
    {
        Mat3f q = Q[i].block<3,3>(0,0);
        Vec3f x0 = Q[i].block<3,1>(0,3);
        x = q * pt3d + x0;
        e.block<2,1>(2*i, 0) = x.block<2,1>(0,0)/x(2) - pts[i];
        J.block<1,3>(2*i, 0) = (x(2) * q.row(0) - x(0) * q.row(2))/ powf(x(2), 2.0);
        J.block<1,3>(2*i+1, 0) = (x(2) * q.row(1) - x(1) * q.row(2)) / powf(x(2), 2.0);
    }
};


void UTILS::Triangulation::Linear(const std::vector<Mat34f>& poses, const std::vector<Vec2f>& pts, Vec3f& pt3d)
{
    Matf A( 2 * pts.size(), 4 );

    for( int i=0; i<pts.size(); i++ )
    {
        float x = pts[i](0);
        float y = pts[i](1);
        const Mat34f pose = poses[i];

        Vec4f p1t = pose.row(0);
        Vec4f p2t = pose.row(1);
        Vec4f p3t = pose.row(2);

        Vec4f a = x * p3t - p1t;
        Vec4f b = y * p3t - p2t;

        a.normalize();
        b.normalize();

        A.block<1, 4>( 2*i, 0)   = a.transpose();
        A.block<1, 4>( 2*i+1, 0) = b.transpose();
    }

    Vec4f X;
    nullspace<Matf, Vec4f>(&A, &X);
    // UTILS::nullspace<Matf, Vec4f>(&A, &X);
    pt3d = X.head(3)/ X(3);

};

void UTILS::Triangulation::Sampson(const std::vector<Mat34f>& poses, const std::vector<Vec2f>& pts, Vec3f& pt3d)
{
    int sz = static_cast<int>(pts.size());
    Linear(poses, pts, pt3d);

    Mat4f T, Ttmp;

    Eigen::Matrix<float, 1, 4> pt3d_homo = pt3d.transpose().homogeneous();
    Mat4f A;
    A.setZero();
    A.block<1,4>(0, 0) = pt3d_homo;
    Eigen::JacobiSVD<Mat4f> asvd(A, Eigen::ComputeFullV);
    T = asvd.matrixV();
    Ttmp.block<4, 3>(0, 0) = T.block<4, 3>(0, 1);
    Ttmp.block<4, 1>(0, 3) = T.block<4, 1>(0, 0);
    T = Ttmp;

    std::vector<Mat34f> Q(sz);

    for( int i=0; i<sz; ++i )
    {
        Q[i] = poses[i] * T;
    }

    // Gauss Newton Section
    Vec3f Y(0, 0, 0);
    Vecf err_prev(2 * sz), err(2 * sz);
    Matf J(2 * sz, 3);
    err_prev.setOnes();
    err_prev *= std::numeric_limits<float>::max();
    const int niters = 30;

    for( int i=0; i<niters; ++i )
    {
        residule(Y, pts, Q, err, J);
        if( 1 - ( err.norm()/err_prev.norm() ) < 10e-8 )
            break;

        err_prev = err;
        Y = Y - (J.transpose() * J).inverse() * (J.transpose() * err);
    }

    Vec4f X = T * Y.homogeneous();
    pt3d = X.head<3>() / X(3);

    bool is_debug = false;
    
    if( is_debug )
    {
        float err = 0.0f;
        for( int i=0; i<sz; i++)
        {
            Vec3f x_homog = poses[i] * X;
            Vec2f dx = (x_homog.head<2>() / x_homog(2)) - pts[i];
            err += sqrtf(dx.dot(dx));
        }

        std::cout << "error (after refinement): " << err << std::endl;
    }

};
