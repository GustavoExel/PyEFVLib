//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Line("SOUTH") = {3};
//+
Physical Line("WEST") = {1};
//+
Physical Line("NORTH") = {4};
//+
Physical Line("EAST") = {2};
//+
Physical Surface("BODY") = {1};
