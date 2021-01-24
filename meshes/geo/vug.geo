// Gmsh project created on Tue Oct 20 14:40:31 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.33333, 0.33333, 0, 1.0};
//+
Point(6) = {0.66667, 0.33333, 0, 1.0};
//+
Point(7) = {0.66667, 0.66667, 0, 1.0};
//+
Point(8) = {0.33333, 0.66667, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 2};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(1) = {3, 4, 1, -2};
//+
Curve Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Curve("west") = {4};
//+
Physical Curve("south") = {1};
//+
Physical Curve("east") = {2};
//+
Physical Curve("north") = {3};
//+
Physical Surface("outer") = {1};
//+
Physical Surface("inner") = {2};
