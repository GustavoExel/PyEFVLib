//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Surface("TOP") = {6};
//+
Physical Surface("BOTTOM") = {5};
//+
Physical Surface("NORTH") = {3};
//+
Physical Surface("SOUTH") = {4};
//+
Physical Surface("EAST") = {1};
//+
Physical Surface("WEST") = {2};
//+
Physical Volume("BODY") = {1};
