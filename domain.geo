// 2D axisymmetric domain as a rectangular
Point(1) = {0, 0, 0};
Point(2) = {12.5/1000, 0, 0 };
Point(3) = {12.5/1000, 12.5/1000, 0};
Point(4) = {0, 1.25/100., 0};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(5) = {1,2,3,4}; 
Plane Surface(1) = {5};


Physical Curve("TopAnode", 11) = {3}; //Anode as a boundary
Physical Curve("BottomCathode",12) = {1};  //Cathode as a boundary


Physical Surface("domain_rec", 2) = {1};

Field[1] = Box;  //Box field for the streaemr present area
Field[1].VIn = 4e-6;
Field[1].VOut = 500e-5;
Field[1].XMax = 0.88e-3;
Field[1].XMin = 0;
Field[1].YMax = 12.5/1000;
Background Field = 1;

Mesh.Algorithm = 5; //2D Delaunay meshing algorithm



