The shape toolbox is a set of routines for manipulating and analyzing shapes, including estimation of the MAP skeleton (see Feldman & Singh, PNAS 2006). 

USING THE SHAPE TOOLBOX

The toolbox can be used in two ways:

1. Run shape_tool at the command line, which starts the shape toolbox GUI. In the GUI, you can draw and modify a shape with the mouse, estimate its MAP skeleton, draw its MAT (medial axis transform, Blum 1973), or draw a skeleton manually and then analyze it with respect to the Feldman & Singh Bayesian framework. The buttons each have tooltips (mouseover to see them) which provide a bit more explanation of what they do. Running shape_tool is the easiest way to use the toolbox, since it doesn't require understanding any of the underlying data structures.

Note: In the GUI, shapes are drawn in black, and manually drawn skeletons are drawn in blue. After you finish drawing a shape, the system automatically switches to drawing a skeleton next. That is, after the shape is complete the next mouse click begins a (blue) skeleton. Of course if you intend to have the program estimate a skeleton (e.g. the MAP skeleton), there's no reason to manually draw a skeleton. Just draw your shape, then click MAP skeleton, and after a short time you the program will draw the estimated skeleton with ribs. By default, MAP skeletons are drawn in "rainbow" style, with a different color for each axis; note that this implies a part decomposition, which is one of the main goals of the theory (see Feldman & Singh 2006 for discussion). Or, if you have a pre-made shape you want to analyze (from a matlab variable or from a file), use load_shape (see below) at the command line to install it into the GUI, then analyze or modify it using the buttons. 

2. Alternatively, one can run several functions from the command line, including mapskeleton, medial_axis, and draw_shape. Each of these takes a shape as an argument, which is an nx2 array of numbers corresponding to the n boundary points of a polygonal approximation to a shape. Use help to see details. 

Quick guide to command-line functions: 

1. draw_shape draws a shape, e.g. draw_shape(myshape);

2. mapskeleton estimates the MAP skeleton, e.g. skeleton = mapskeleton(myshape);

3. medial_axis computes the Medial Axis Transform (MAT), e.g. draw_shape(myshape); medial_axis(myshape,1);

4. draw_skeleton draws a skeleton, e.g. skel = mapskeleton(myshape); draw_skeleton(skel);

5. load_shape takes a shape from a file or from a variable in the matlab workspace and loads it into the shape_tool GUI. sample_shape.txt is an example file in the correct format. Note the oridinary matlab command load is sufficient to load a shape stored in a file into the workspace (see example below); load_shape is specifically for loading one into the shape_tool GUI. 

For example, to load the sample shape, compute its MAP skeleton, and draw the shape and skeleton:

>> s = load('sample_shape.txt'); draw_shape(s); draw_skeleton(mapskeleton(s));

Generally, the functions visible at the top level of the shape_toolbox directory are intended to be run from the command line, while the functions in the subdirectory Functions are not intended to be called directly. However once you understand the data structures they use, some of them can be useful. All have documentation in "help." 

NOTE ON THE MAT:

Note the difference between medial_axis.m and Functions/make_medial_axis.m.  medial_axis.m simply computes the MAT and (optionally) draws it. For example, to draw the MAT of the sample shape, use:

>> s = load('sample_shape.txt'); draw_shape(s); medial_axis(s,1);

make_medial_axis.m is a much more complicated (and slower) procedure to compute the MAT and hierarchically organize it into a skeleton (in the format used for skeletons elsewhere in the toolbox), and thus includes a great deal of computation beyond the simple MAT. For example, the simple MAT doesn't really have distinct "axes" (the MAT itself in Blum's definition only has points), but the output of make_medial_axis.m does. 

You can play around with shapes in the GUI and display their MATs using the MAT button, which uses make_medial_axis.m.

COMPETENCE VS. PERFORMANCE:

Note that for some shapes mapskeleton yields a "bad" skeletal estimate, though in our view this is due to deficiencies in the estimation procedure rather than in the underlying theory. For example, rectangles often have skeletons that "fork" at one or both ends. However, if you manually draw a skeleton consisting of a single line segment along the major axis, then optimize it, you can easily confirm that its DL is better (lower) than the result of the mapskeleton estimation procedure. 

NOTES ON DATA STRUCTURES:

Throughout the toolbox, a shape is an nx2 matrix of numbers, defining the points on the outer boundary of the shape. Generally shapes are assumed to be closed, meaning that shapes are drawn so that the last point connects to the first point.

A skeleton is a vector of axes, each of which is a structure with fields:
 skeleton(i).index
 skeleton(i).parent:  the index of the axis to which this axis is attached
 skeleton(i).contour: the points along this axis (an open curve)

Ribs (sometimes called coribs) are correspondences between a shape and a skeleton. See Functions/compute_coribs for details. 

GLOBAL VARIABLES: Some variables used by the shape_tool GUI are stored as globals and hence are available at the command line. For example, the shape currently displayed is stored in (global) current_shape; the skeleton just estimated (either by the MAPskeleton button or the medial axis button) is in (global) current_skeleton; and the associated ribs are in (global) current_coribs. So for example if you analyze a shape in the GUI and estimate its skeleton, and then want to save its skeleton for future use, you can access the skeleton by setting current_skeleton to global.

KNOWN BUGS/ISSUES:

- Mapskeleton and other functions do not work correctly on shapes with self-intersections. However, the system does not currently implement a check for self-intersections. (Algorithms to check are O(n log n) and would slow the system down noticeably.) Be aware that some subtle and peculiar misbehavior arises due to self-intersections, even very small ones. So when drawing a shape manually, be sure it does not intersect itself. "Fractalizing" a shape in the GUI can accidentally introduce small self-intersections. Random_axial_shape sometimes creates small self-intersections (usually near cusps) and these will lead to aberrant behavior. 

- When you draw a skeleton manually in the GUI (see above), the second and subsequent axis are supposed to automatically connect to the nearest existing axis in order to create a connected, hierarchically organized skeleton. However it takes a short time to find the nearest point on an axis, and if you begin to move the mouse too quickly before it has completed, the axis won't connect correctly. Workaround: when clicking the first point of the second or subsequent axis, pause for a moment (a half second on my system) before beginning to move the mouse. The axis should connect correctly.

- The time mapskeleton takes to run depends on the number of points in the shape, and it can be slow for shapes with many points. Drawing a shape in the GUI with default constants yields a shape with about 150-200 points, and mapskeleton has been tested primarily on shapes of about that size. With shapes of 1000 points for example it is quite a bit slower. But note that "extra" points often do not affect the appearance of the shape, so most shapes can be down-sampled to a smaller number of points (e.g. with Functions/resample_shape) without any perceptual difference whatsoever. If your shape takes a long time (> 1 minute on a recent machine) try resampling it (e.g. mapskeleton(resample_shape(myshape,200)).

GENERAL CAVEAT: The routines in this toolbox are hastily written and thinly documented, and contain numerous examples of poor programming practice. We know. Their main purpose is to demonstrate and facilitate scientific ideas described in the 2006 paper, not to efficiently implement a well-established algorithm. Version 1.0 attempts to clean these problems up enough for public distribution, but many fundamental problems remain. Needless to say, many theoretical ideas developed in subsequent work have not yet been included in this release. 

VERSION HISTORY:
- Version 0.1: Original. Written during 2004-6 by Jacob Feldman and Manish Singh.
- Version 0.2-4. Sept 2006-10. Added several functions. Thanks to John Wilder and Seha Kim for various contributions. 
- Version 1.0: Feb 2016. Many, many bug fixes. Cleaned up code and added this documentation. 
