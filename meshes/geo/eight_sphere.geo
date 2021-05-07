//+
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
//+
Box(2) = {-1, -1, -1, 1, 2, 2};
//+
Box(3) = {0, -1, -1, 1, 2, 1};
//+
Box(4) = {0, -1, 0, 1, 1, 1};
//+
BooleanDifference{ Volume{1}; Delete; }{ Volume{2}; Volume{4}; Volume{3}; Delete; }
//+
Physical Surface("Outer") = {1};
//+
//+
Physical Surface("InnerZ") = {3};
//+
Physical Surface("InnerY") = {4};
//+
Physical Surface("InnerX") = {2};
//+
Physical Volume("Body") = {1};
