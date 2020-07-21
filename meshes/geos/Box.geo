//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Surface("Top") = {6};
//+
Physical Surface("Bottom") = {5};
//+
Physical Surface("North") = {4};
//+
Physical Surface("South") = {3};
//+
Physical Surface("West") = {1};
//+
Physical Surface("East") = {2};
//+
Physical Volume("Body") = {1};
