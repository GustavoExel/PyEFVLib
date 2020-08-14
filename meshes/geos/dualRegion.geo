//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 0.5, 1, 0};
//+
Rectangle(2) = {0.5, 0, 0, 0.5, 1, 0};
//+
Physical Surface("Body1") = {1};
//+
Physical Surface("Body2") = {2};
//+
Physical Curve("NW") = {3};
//+
Physical Curve("NE") = {7};
//+
Physical Curve("W") = {4};
//+
Physical Curve("E") = {6};
//+
Physical Curve("SW") = {1};
//+
Physical Curve("SE") = {5};
//+
Physical Curve("M") = {2};
