//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {1, 10, 0, 1.0};
//+
Point(5) = {0.5, 10, 0, 1.0};
//+
Point(6) = {0, 10, 0, 1.0};
//+
Line(1) = {5, 6};
//+
Line(2) = {6, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 2};
//+
Curve Loop(1) = {2, 3, -7, 1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 4};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Body") = {1, 2};
//+
Physical Curve("North") = {1, 6};
//+
Physical Curve("South") = {3, 4};
//+
Physical Curve("West") = {2};
//+
Physical Curve("East") = {5};
