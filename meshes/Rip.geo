//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 1, 0, 1.0};
//+
Point(3) = {1, 0, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 0.48, 0, 1.0};
//+
Point(6) = {0, 0.52, 0, 1.0};
//+
Point(7) = {0.2, 0.5, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 4};
//+
Line(4) = {4, 6};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 5};
//+
Line(7) = {5, 1};
//+
Curve Loop(1) = {3, 4, 5, 6, 7, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("NORTH") = {3};
//+
Physical Curve("SOUTH") = {1};
//+
Physical Curve("EAST") = {2};
//+
Physical Curve("WEST") = {4, 7};
//+
Physical Curve("RIP") = {5, 6};
//+
Physical Surface("BODY") = {1};
