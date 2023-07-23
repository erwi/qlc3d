// Define structure dimensions
width = 5;      // cell width along x-axis
depth = 0.1;    // cell depth along y-axis
height = 2;     // cell height along z-axis

// additional x-coordinates for start and end of electrod on bottom surface
electrode_width = 1;
electrode_xmin = 2;
electrode_xmax = electrode_xmin + electrode_width;

// Points front surface, y = 0
Point(1) = {0, 0, 0};
Point(2) = {electrode_xmin, 0, 0};
Point(3) = {electrode_xmax, 0, 0};

Point(4) = {width, 0, 0};
Point(5) = {width, 0, height};
Point(6) = {0, 0, height};

// Points back surface, y = depth
Point(7) = {0, depth, 0};
Point(8) = {electrode_xmin, depth, 0};
Point(9) = {electrode_xmax, depth, 0};

Point(10) = {width, depth, 0};
Point(11) = {width, depth, height};
Point(12) = {0, depth, height};

// lines front surface
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// lines back surface
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 7};

// lines left side
Line(13) = {1, 7};
Line(14) = {6, 12};

// lines right side
Line(15) = {4, 10};
Line(16) = {5, 11};

// lines bottom surface
Line(17) = {2, 8};
Line(18) = {3, 9};

FRONT_SURFACE = 1;
Curve Loop(FRONT_SURFACE) = {1, 2, 3, 4, 5, 6}; // front surface
Plane Surface(FRONT_SURFACE) = {1};

BACK_SURFACE = 2;
Curve Loop(BACK_SURFACE) = {7, 8, 9, 10, 11, 12}; // back surface
Plane Surface(BACK_SURFACE) = {2};

LEFT_SURFACE = 3;
Curve Loop(LEFT_SURFACE) = {13, -12, -14, 6}; // left surface
Plane Surface(LEFT_SURFACE) = {LEFT_SURFACE};

RIGHT_SURFACE = 4;
Curve Loop(RIGHT_SURFACE) = {15, 10, -16, -4}; // right surface
Plane Surface(RIGHT_SURFACE) = {RIGHT_SURFACE};

TOP_SURFACE = 5;
Curve Loop(TOP_SURFACE) = {14, -11, -16, 5}; // top surface
Plane Surface(TOP_SURFACE) = {TOP_SURFACE};

BOTTOM_SURFACE1 = 6;
Curve Loop(BOTTOM_SURFACE1) = {13, 7, -17, -1}; // bottom surface 1
Plane Surface(BOTTOM_SURFACE1) = {BOTTOM_SURFACE1};

BOTTOM_SURFACE2 = 7;
Curve Loop(BOTTOM_SURFACE2) = {17, 8, -18, -2}; // bottom surface 2 (electrode)
Plane Surface(BOTTOM_SURFACE2) = {BOTTOM_SURFACE2};

BOTTOM_SURFACE3 = 8;
Curve Loop(BOTTOM_SURFACE3) = {18, 9, -15, -3}; // bottom surface 3 (electrode)
Plane Surface(BOTTOM_SURFACE3) = {BOTTOM_SURFACE3};

// Ensure left/right and front/back surface meshses are periodic
// During mesh generation, the mesh on RIGHT_SURFACE will be crated by copying
// the mesh from the LEFT_SURFACE
Periodic Surface{RIGHT_SURFACE} = {LEFT_SURFACE} Translate {width, 0, 0};
Periodic Surface{BACK_SURFACE} = {FRONT_SURFACE} Translate {0, depth, 0};


// Define single LC volume
LC_VOLUME = 1;
Surface Loop(LC_VOLUME) = {1, 2, 3, 4, 5, 6, 7, 8};
Volume(LC_VOLUME) = {LC_VOLUME};

// Materials names for every surface and volume
Physical Surface("periodic", 1)     = {LEFT_SURFACE, RIGHT_SURFACE, FRONT_SURFACE, BACK_SURFACE};
Physical Surface("FixLC1", 2)       = {TOP_SURFACE};
Physical Surface("Electrode1", 3)   = {TOP_SURFACE};
Physical Surface("FixLC2", 4)       = {BOTTOM_SURFACE1, BOTTOM_SURFACE2, BOTTOM_SURFACE3};
Physical Surface("Electrode2", 5)   = {BOTTOM_SURFACE2};

Physical Volume("domain1", 6)       = {LC_VOLUME};


// Create 3D mesh. Set the element size equal to the mesh depth so that the mesh
// is a single element deep.
Mesh.MeshSizeMin = depth;
Mesh.MeshSizeMax = depth;
Mesh 3;

