// Generate a quarter-circle mesh

// Define parameters
radius = 100e-9;
meshSize = 10e-9;  // Mesh element size

// Points for the quarter-circle
Point(1) = {0, 0, 0, meshSize}; 
Point(2) = {radius, 0, 0, meshSize};
Point(3) = {0, radius, 0, meshSize};

// Create the quarter-circle arc and lines
Circle(1) = {2, 1, 3};
Line(2) = {3, 1};
Line(3) = {1, 2};

// Create the surface
Line Loop(4) = {2, 1, 3};
Plane Surface(5) = {4};

// Physical groups for boundaries
Physical Line("arc") = {1}; // Boundary ID for the arc
Physical Line("vertical") = {2}; // Boundary ID for y-axis
Physical Line("horizontal") = {3}; // Boundary ID for x-axis

// Physical group for the surface
Physical Surface("quarter_circle") = {5};


Mesh 2;
Save "mesh1.msh";