// Gmsh project created on Thu Dec 10 10:52:36 2020
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 0.5, 0, 1.0};
//+
Point(4) = {0, 0.5, 0, 1.0};
//+
Point(5) = {-0.6, 0.6, -0, 1.0};
//+
Recursive Delete {
  Point{5}; 
}
//+
Circle(1) = {1, 4, 2};
//+
Point(5) = {-0.5, 0.5, -0, 1.0};
//+
Line(2) = {1, 5};
//+
Line(3) = {5, 4};
//+
Line(4) = {4, 2};
//+
Curve Loop(1) = {4, -1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("Body") = {1};
//+
Physical Curve("Boundary") = {1, 3, 4, 2};
