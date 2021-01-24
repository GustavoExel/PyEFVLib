//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 0.2, 0};
//+
Physical Line("South") = {3};
//+
Physical Line("West") = {1};
//+
Physical Line("North") = {4};
//+
Physical Line("East") = {2};
//+
Physical Surface("Body") = {1};
