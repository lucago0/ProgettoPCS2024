 #include <gtest/gtest.h>
 #include "GeometryLibrary.hpp"

using namespace std;
using namespace fracturesLib;

double tol = pow(10,-10);

 TEST(areCloseTest, areClose) {
     Fractures fractures;
     Matrix<double,3,4> poly1;
     Matrix<double,3,4> poly2;
     poly1 << 0.68,0.21,0.08,0.55,
         0.52,0.99,0.86,0.39,
         0.19,0.19,0.6,0.6;
     poly2 << 1.17,0.19,0.09,1.07,
         0.85,1.05,0.54,0.34,
         0.02,0.02,0.44,0.44;
     fractures.vertices.resize(2);
     fractures.vertices[0] = poly1;
     fractures.vertices[1] = poly2;;
     unsigned int id1 = 0;
     unsigned int id2 = 1;

     EXPECT_TRUE(areClose(fractures,id1, id2,tol));
 }

 TEST(areCloseTest, areClose2) {
     Fractures fractures;
     Matrix<double,3,4> poly1;
     Matrix<double,3,4> poly2;
     poly1 << 5.31,5.31,6.72,6.72,
         0,2,2,0,
         2.48,2.48,1.07,1.07;
     poly2 << 5,5,6.42,6.42,
         0,2,2,0,
         6.31,6.31,7.72,7.72;
     fractures.vertices.resize(2);
     fractures.vertices[0] = poly1;
     fractures.vertices[1] = poly2;;
     unsigned int id1 = 0;
     unsigned int id2 = 1;

     EXPECT_FALSE(areClose(fractures,id1, id2,tol));
 }

 TEST(planeTest, plane){
     Fractures fractures;
     Matrix<double,3,4> poly1;
     poly1 << 0.68,0.21,0.08,0.55,
         0.52,0.99,0.86,0.39,
         0.19,0.19,0.6,0.6;
     fractures.vertices.resize(1);
     fractures.vertices[0] = poly1;

     unsigned int id = 0;
     Vector4d expected = {0.1926999999999, 0.1927, 0.1222, -0.254458};
     EXPECT_TRUE(plane(id, fractures).isApprox(expected, tol));
 }

 TEST(planeTest, plane2){
     Fractures fractures;
     Matrix<double,3,4> poly1;
     poly1 << 5.31,5.31,6.72,6.72,
         0,2,2,0,
         2.48,2.48,1.07,1.07;
     fractures.vertices.resize(1);
     fractures.vertices[0] = poly1;

     unsigned int id = 0;

     Vector4d expected = {-2.82, 0, -2.82, 21.967799999999};
     EXPECT_TRUE(plane(id, fractures).isApprox(expected, tol));
 }

 TEST(planesIntersectionTest, planesIntersection){
     Vector4d plane1 = {0.19, 0.19, 0.12, -0.25};
     Vector4d plane2 = {0.09, 0.42, 0.52, -0.46};

     line r;
     r.point = {0.2807017543859,1.0350877192982,0};
     r.direction = {-0.0484000000000,0.0880000000000,-0.062699999999};

     EXPECT_TRUE(planesIntersection(plane1,plane2,tol).point.isApprox(r.point,tol));
     EXPECT_TRUE(planesIntersection(plane1,plane2,tol).direction.isApprox(r.direction,tol));
 }


 TEST(linesIntersectionTest, linesIntersection){
     line r1;
     r1.point = {1,0,1};
     r1.direction = {-1,1,2}; //Y-Z
     line r2;
     r2.point = {0,-2,0};
     r2.direction = {1,2,1};

     Vector3d M = r1.point+linesIntersection(r1,r2)[3]*r1.direction;
     Vector3d expected = {1,0,1}; //(x,y,z)
     EXPECT_TRUE(M.isApprox(expected,tol));
     EXPECT_NEAR(linesIntersection(r1,r2)[4],1,tol);
 }

  TEST(intervalsIntersectionTest, intervalsIntersection){
      Vector4d Q = {-1, 1, 0, 2};
      Vector4d expected = {0,1,-1,2}; //(x,y,z,t,s)

      EXPECT_TRUE(intervalsIntersection(Q,tol).isApprox(expected,tol));
     }


 int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
 }

