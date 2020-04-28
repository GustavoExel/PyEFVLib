//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Curve("NORTH") = {3};
//+
Physical Curve("SOUTH") = {1};
//+
Physical Curve("WEST") = {4};
//+
Physical Curve("EAST") = {2};
//+
Physical Surface("BODY") = {1};
