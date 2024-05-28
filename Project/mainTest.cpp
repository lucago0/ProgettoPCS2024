 #include <gtest/gtest.h>
 #include "GeometryLibrary.hpp"

 TEST(AreCloseTest, areClose) {
     Fractures frac;
     MatrixXd poly1; MatrixXd poly2;
     poly1 << 0.68,0.21,0.08,0.55,
         0.52,0.99,0.86,0.39,
         0.19,0.19,0.6,0.6;
     poly2 << 1.17,0.19,0.09,1.07;
         0.85,1.05,0.54,0.34;
         0.02,0.02,0.44,0.44;
     frac.Vertices.insert(make_pair(0,poly1)); frac.Vertices.insert(make_pair(1,poly2));

     EXPECT_TRUE(areClose(frac,poly1, poly2));
 }

 TEST(PianoTest, Piano){
     Fractures frac;
     MatrixXd poly1;
     poly1 << 0.68,0.21,0.08,0.55,
         0.52,0.99,0.86,0.39,
         0.19,0.19,0.6,0.6;
     frac.Vertices.insert(make_pair(0,poly1));

     Vector4d expected = {0.19, 0.19, 0.12, -0.25};
     EXPECT_TRUE(Piano(0,frac)==expected);
 }

 TEST(InterTest, Inter){
     Vector4d plane1 = {0.19, 0.19, 0.12, -0.25};
     Vector4d plane2 = {0.09, 0.42, 0.52, -0.46};

     Line r;
     r.point = {0.26,1.07,-0.01};
     r.direction = {0.05,-0.09,0.06};

     EXPECT_TRUE(Inter(plane1,plane2).point==r.point);
     EXPECT_TRUE(Inter(plane1,plane2).direction==r.direction);
 }


 TEST(PuntiIntersRettaTest, PuntiIntersRetta){
     Line r1;
     r1.point = {0.09,0.54,0.44};
     r1.direction = {0.98,-0.2,0}; //Y-Z
     Line r2;
     r2.point = {0.26,1.07,-0.01};
     r2.direction = {0.05,-0.09,0.06};

     Vector3d M = r1.point+PuntiIntersRetta(r1,r2)[3]*r1.direction;
     Vector3d expected = {0.59,0.44,0.44}; //(x,y,z,t,s)
     EXPECT_TRUE(M==expected);
     EXPECT_TRUE(PuntiIntersRetta(r1,r2)[4]==209./391);
 }

 TEST(IntersectionTest, intersection){
     Vector4d Q = {-1, 1, 0, 2};
     VectorXd expected = {0,1,-1,2}; //(x,y,z,t,s)

     EXPECT_TRUE(intersection(Q)==expected);
 }

 TEST(IntersectionTest, intersection){
     Vector4d Q = {-1, 1, 0, 2};
     VectorXd expected = {0,1,-1,2}; //(x,y,z,t,s)

     EXPECT_TRUE(intersection(Q)==expected);
 }

 int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
     return RUN_ALL_TESTS();
 }
