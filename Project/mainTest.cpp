 #include <gtest/gtest.h>
 #include "GeometryLibrary.hpp"

using namespace std;
using namespace FracturesLib;

double tol = pow(10,-10);

 TEST(AreCloseTest, areClose) {
     Fractures frac;
     Matrix<double,3,4> poly1;
     Matrix<double,3,4> poly2;
     poly1 << 0.68,0.21,0.08,0.55,
         0.52,0.99,0.86,0.39,
         0.19,0.19,0.6,0.6;
     poly2 << 1.17,0.19,0.09,1.07,
         0.85,1.05,0.54,0.34,
         0.02,0.02,0.44,0.44;
     frac.Vertices.insert(make_pair(0,poly1));
     frac.Vertices.insert(make_pair(1,poly2));
     unsigned int id1 = 0;
     unsigned int id2 = 1;

     EXPECT_TRUE(areClose(frac,id1, id2,tol));
 }

 TEST(AreCloseTest, areClose2) {
     Fractures frac;
     Matrix<double,3,4> poly1;
     Matrix<double,3,4> poly2;
     poly1 << 5.31,5.31,6.72,6.72,
         0,2,2,0,
         2.48,2.48,1.07,1.07;
     poly2 << 5,5,6.42,6.42,
         0,2,2,0,
         6.31,6.31,7.72,7.72;
     frac.Vertices.insert(make_pair(0,poly1));
     frac.Vertices.insert(make_pair(1,poly2));
     unsigned int id1 = 0;
     unsigned int id2 = 1;

     EXPECT_FALSE(areClose(frac,id1, id2,tol));
 }

 TEST(PianoTest, Piano){
     Fractures frac;
     Matrix<double,3,4> poly1;
     poly1 << 0.68,0.21,0.08,0.55,
         0.52,0.99,0.86,0.39,
         0.19,0.19,0.6,0.6;
     frac.Vertices.insert(make_pair(0,poly1));

     unsigned int id = 0;
     Vector4d expected = {0.1926999999999, 0.1927, 0.1222, -0.254458};
     EXPECT_TRUE(Piano(id, frac).isApprox(expected, tol));
 }

 TEST(PianoTest, Piano2){
     Fractures frac;
     Matrix<double,3,4> poly1;
     poly1 << 5.31,5.31,6.72,6.72,
         0,2,2,0,
         2.48,2.48,1.07,1.07;
     frac.Vertices.insert(make_pair(0,poly1));

     unsigned int id = 0;

     Vector4d expected = {-2.82, 0, -2.82, 21.967799999999};
     EXPECT_TRUE(Piano(id, frac).isApprox(expected, tol));
 }

 TEST(InterTest, Inter){
     Vector4d plane1 = {0.19, 0.19, 0.12, -0.25};
     Vector4d plane2 = {0.09, 0.42, 0.52, -0.46};



     Line r;
     r.point = {0.28070175438,1.03508771929,0};
     r.direction = {-0.04840000000,0.088000000000,-0.062699999999};

     EXPECT_TRUE(Inter(plane1,plane2,tol).point.isApprox(r.point,tol));
     EXPECT_TRUE(Inter(plane1,plane2,tol).direction.isApprox(r.direction,tol));
 }


 TEST(PuntiIntersRettaTest, PuntiIntersRetta){
     Line r1;
     r1.point = {1,0,1};
     r1.direction = {-1,1,2}; //Y-Z
     Line r2;
     r2.point = {0,-2,0};
     r2.direction = {1,2,1};

     Vector3d M = r1.point+PuntiIntersRetta(r1,r2)[3]*r1.direction;
     Vector3d expected = {1,0,1}; //(x,y,z)
     EXPECT_TRUE(M.isApprox(expected,tol));
     EXPECT_NEAR(PuntiIntersRetta(r1,r2)[4],1,tol);
 }

  TEST(IntersectionTest, intersection){
      Vector4d Q = {-1, 1, 0, 2};
      Vector4d expected = {0,1,-1,2}; //(x,y,z,t,s)

      EXPECT_TRUE(intersection(Q,tol).isApprox(expected,tol));
     }


 int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
 }

