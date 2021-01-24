//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Line("North") = {3};
//+
Physical Line("South") = {1};
//+
Physical Line("West") = {4};
//+
Physical Line("East") = {2};
//+
Physical Surface("Body") = {1};
