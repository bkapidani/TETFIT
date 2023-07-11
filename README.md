# TETFIT

<h1>The TetFIT toolbox</h1>
<p>TetFIT a simulation toolbox which was developed during my Ph.D. thesis work at the University of Udine, which now has been superseeded by arbitrary polynomial degree versions available in Netgen/NGSolve (see <a href="https://ngsolve.org">[https://ngsolve.org]</a>). The source-code is written in C++ and works under Unix architectures (tested on Debian 9), MacOs and on Windows, compiled in the CygWin Posix compatibility layer or using WSL.
All the documentation that follows is taken from chapter 7 of my final thesis report available <a href="https://air.uniud.it/bitstream/11390/1142992/2/thesis_kapidani_pdfA.pdf">here</a>.</p>

<h2>Installation</h2>
<p>After cloning the git repository one should use simply <code>cmake</code> and <code>make -j install</code> with the appropriate targets to get the executable compiled.</p>

<h2>The user interface</h2>
<p>The executable is called on the terminal with the instruction:</p>

<pre><code>tetfit simulation.fdtd
</code></pre>
<p>where <code>simulation.fdtd</code> is an input file written in a scripting language developed by the author. An example of an input file is the following:</p>

<pre><code>########## file simulation.fdtd

########## mesh
DEFINE mesh 1
   SET type           tetrahedral
   #SET xgrid         {0,0.01,1}
   #SET ygrid         {0,0.01,1}
   #SET zgrid         {0,0.01,1}
   SET mesher         netgen
   SET file           box.mesh
   SET name	          waveguide
   SET scalefactor    1
END mesh 1

########## Materials definitions
DEFINE material 1
   SET epsilon      1
   SET mu           1
   SET sigma        0
   SET chi          0
END material 1

########## Boundary conditions
DEFINE bc 1
   SET surface 2
   SET surface 3
   SET surface 4
   SET surface 5
   SET surface 6
   SET type pec
END bc 1

########## Sources
DEFINE source 1
   SET surface            1
   SET type               h # magnetic field
   SET profile            wave
   SET mode               { sin , cos , cos }
   SET center             {0, 0, 0} # x-xcenter, y-ycenter, z-zcenter
   SET direction          x
   SET amplitude          1
   SET frequency          2e8
   SET carrier            sin
   SET wavevector         { 0.5, 0, 0 }
END source 1

########## Output setups
DEFINE output 1
   SET name               zsections
   SET mode               probepoint
   SET xgrid              { 0.1, 0.1,0.6}
   SET ygrid              { 0.5, 1, 0.5}
   SET zgrid              { 0, 0.04, 1.05}
   SET grid               on
   SET period             0
END output    1

######## Simulations
DEFINE simulation 1
   SET method          dga
   SET solver          cg
   SET tolerance       1e-8
   SET mesh            1
   SET source          1
   SET duration        4e-8
   SET output          1
   SET courant         0.98
END simulation    1
</code></pre>
<p>First of all the <code>#</code> token is reserved to indicate comments, which can start anywhere on a line (as can be seen in this particular example) but end only with an EOL character. The scripting language should be very straightforward to understand due to its simplicity. There are in fact just three types of instructions, terminated by end-of-line (EOL) characters (or by comments) and with the following syntax:</p>
<ul>
<li><code>DEFINE primitive label</code>: defines an object of type <code>primitive</code> with an integer <code>label</code> attached to it.</li>
<li><code>END primitive label</code>: signals the end of the definition of an object of type <code>primitive</code>. The integer <code>label</code> must match the one of the associated <code>DEFINE</code> instruction.</li>
<li><code>SET property value</code>: sets a property of an object, all <code>SET</code> instructions must be in between a <code>DEFINE</code> instruction and the matching <code>END</code> instruction.</li>
</ul>
<p>There are six (plus one) types of primitives. To illustrate them more clearly e we reintroduce an annotated declaration example for every single one.</p>

<h3>The mesh primitive</h3>
<pre><code>DEFINE mesh 1
   SET type        tetrahedral # alternative: cartesian
   SET xgrid       {0,0.01,1}  # only for cartesian, value={xmin,step,xmax}
   SET ygrid       {0,0.01,1}  # only for cartesian, value={ymin,step,ymax}
   SET zgrid       {0,0.01,1}  # only for cartesian, value={zmin,step,zmax}
   SET mesher      netgen      # alternative: gmsh
   SET file        box.mesh    # input file (relative path)
   SET name	       waveguide   # just a name for the mesh
   SET scalefactor 1           # scales the mesh, default unit [m]
END mesh 1
</code></pre>
<p>The <code>mesher</code> field tells TetFIT which syntax it must expect from the mesh input file, given by the <code>file</code> field. For string values <code>netgen</code> and <code>gmsh</code> the respective 3D neutral mesh format provided by the two meshers (see <a href="http://gmsh.info/">[gmsh]</a> and <a href="http://sourceforge.net/projects/netgen-mesher/">[netgen]</a>) is assumed. If the mesh is <code>cartesian</code>, an in-house developed mesher is used and the value for the <code>mesher</code> property is ignored. In that case the 3-vectors given as <code>value</code> fields of <code>xgrid</code>, <code>ygrid</code>, <code>zgrid</code> define the computational background domain (more on this later).</p>

<p>There is a further hidden primitive type linked with the mesh primitive: the <code>solid</code> primitive, which is used only in the FDTD method. When the user wants to use a Cartesian orthogonal grid, an in-house mesher is called instead of the one set by the <code>mesher</code> field. Furthermore, in this case the <code>file</code> field must point to the name of a further input file, which is written with the same syntax of the input file, but can only contain object definitions of type <code>solid</code>:</p>
<pre><code>#### solids file input
DEFINE solid 1
	SET material 1 # the material of the solid
	SET type box  # alternative: sphere, cylinder
	SET corner { -10,-10,-10} # coordinates of the lower left corner of the box
	SET size 	{ 20, 20, 20} # dimensions
END solid 1
</code></pre>
<pre><code>
DEFINE solid 2 # overrides solid 1 where they intersect
	SET material 2
	SET type sphere
	SET center { 0.025,0.0125,0.05 } #obvious
	SET radius 0.005 #obvious
END solid 2
</code></pre>
<p>This example was used in the conductive sphere test-case of Section 1.3. Available types of solid are <code>box</code>, <code>sphere</code>, <code>cylinder</code>. We remark that, for all solid types, the actual simulated solid is given by its intersection with the computational domain (defined in the main script by <code>xgrid</code>, <code>ygrid</code>, <code>zgrid</code>). On the other hand if no solid definition fills a given region of the computational domain, that region is excluded from the simulation (it is given default material 0).</p>

<h3>The material primitive</h3>
<pre><code>DEFINE material 1
   SET epsilon      1 # relative dielectric permittivity
   SET mu           1 # relative magnetic permeability
   SET sigma        0 # electric conductivity
   SET chi          0 # magnetic conductivity
END material    1
</code></pre>
<p>The material label should match a material label present in the mesh input file (both NETGEN and GMSH allow defining material labels). If this is not the case, the declaration is ignored. The default values for the four properties are given in the above example.</p>

<h3>The boundary condition primitive</h3>
<pre><code>DEFINE bc 1
   SET surface    2   # associated surface in mesh file
   SET surface    3   # multiple surfaces are possible
   SET surface    4
   SET surface    5
   SET surface    6   
  #SET materials  {1,0} # 0 is default background material
   SET type       pec # alternative pmc, free (default)
END bc 1
</code></pre>
<p>The value of the <code>surface</code> field can be set many times to different values without overwriting anything (it is just an <code>std::set</code> insertion, if the reader is familiar with C++ programming). If the value does not match a physical surface label present in the mesh file, it is simply ignored.
When the mesh is Cartesian orthogonal surface labels from 1 to 6 are reserved for the plane limiting faces of the computational domain: surface 1 is the one with constant <code>z=zmin</code>, surface 2 denotes <code>x=xmin</code>, surface 3 denotes <code>y=ymin</code>, surface 4 denotes <code>x=xmax</code>, surface 5 denotes <code>y=ymax</code> and surface 6 denotes <code>z=zmax</code>.</p>

<p>Another way to set boundary conditions is by setting the <code>materials</code> field to a 2-vector containing the label of two materials, which applies the boundary condition at their interface. This is especially useful in FDTD simulations in which the domain we want to simulate is not a cuboid or in simulations in which we want to exclude some conductive volume. Values <code>pec</code> and <code>pmc</code> stand for Perfect Electric Conductor and Perfect Magnetic Conductor, respectively. PML absorbing boundary conditions are under development.</p>

<h3>The source primitive</h3>
<pre><code>DEFINE source 1
   SET type       h # magnetic field, alternative: e
   SET amplitude  1 # V/m if type is e, A/m if type is h
   SET mode       { sin , cos , cos } # x,y,z basis functions
   SET wavevector {0.5, 0, 0} # in 1/meters
   SET center     {0, 0, 0} # Spatial source origin
   SET delay      0 # in seconds, works as a unit step function in time
   SET frequency  2e+8 # in Hz
   SET carrier    sin # alternative gaussian or cos
  #SET width      1e-10 # in seconds, only for gaussian carriers
   SET direction  x # alternative y,z,r,phi
   SET surface    1 # label of the mesh surface on which the wave impinges 
END source 1
</code></pre>
<p>The source primitive is a bit heavier on the eyes. Suffices to say that for this particular example, the source we set amounts to forcing a tangential magnetic field</p>

$$
\vec{h} \left(x,y,z=0,t\right)  = \left( 1 \text{A/m} \right) \times \text{sin} \left(2\pi (x-0)\times 1/2\right) \times \text{cos} \left(2\pi (y-0)\times  0\right) \times \text{cos} \left(2\pi (z-0)\times  0\right) \times \text{sin} \left(2\pi t\times  2\times 10^8\right) \times \Theta(t-0) \times \hat{\mathbf{x}}
$$

or simply:

$$
\vec{h} \left(x,y,z=0,t\right) = \text{sin} (\pi x) \text{sin} (2\pi f t) \hat{\mathbf{x}}
$$

<p>where \(f= 200\) MHz. An important feature is that the vector direction of each source field must be aligned with one of the Cartesian axes, but sources can be linearly combined on the same target surface, to obtain (the tangent field of) an arbitrarily impinging waveguide mode. It is also interesting to note that the <code>wavevector</code> field combines with the trigonometric functions in the <code>mode</code> field.</p>

<h3>The output primitive</h3>
<pre><code>DEFINE output 1
   SET name       zsections # becomes the output file name preamble
   SET mode       probepoint # alternative: silo
   SET xgrid      { 0.1, 0.1,0.6} # as in the cartesian mesh definition
   SET ygrid      { 0.5, 1, 0.5} # as in the cartesian mesh definition
   SET zgrid      { 0, 0.04, 1.05} # as in the cartesian mesh definition
   SET grid       on # switches on the grid of probes
  #SET probe      {0.5,0.5,0.5} # alternative, xyz probe coordinates 
   SET period     0 # in seconds!
   SET axes       cartesian # alternative: cylindrical, spherical
END output    1
</code></pre>
<p>The output primitive defines the type of output we want to obtain. In general <code>probepoint</code> is the most used and works in two flavors: either it defines a three-dimensional structured grid of points on which it interpol

ates and stores the fields at various time-steps (with a user-defined period which is a true time, not a period in number of time-steps), or, if the grid is not <em>switched on</em><sup>[^1^]</sup>, multiple <code>SET probe</code> instructions can be used to measure the fields at arbitrary points. A routine, inspired by this <a href="http://theory.cm.utexas.edu/jacob/cpmd2001/node20.html#SECTION00083000000000000000">book</a>, to directly compute the Discrete Fourier Transform of the field at a particular point is also available as <code>SET fprobe</code>. An alternative output mode is <code>mode silo</code>, which instead yields the fields (and the mesh) in the widely used open-source <em>Silo</em> format<a href="http://wci.llnl.gov/codes/silo/">[silo]</a>, which can be very useful both for animated field time evolutions and for debugging. Finally, via <code>SET axes</code>, it is always possible to obtain field values in polar (spherical or cylindrical) coordinates, even when it does not make any practical sense.</p>

<h3>The simulation primitive</h3>
<pre><code>DEFINE simulation 1
   SET method      dga  # alternative: fem, fdtd, default: dga
   SET solver      cg   # only for fem! alternative: agmg
   SET tolerance   1e-8 # accepted relative residual
   SET mesh        1    # the mesh used by this simulation
   SET source      1    # the source used, possibly more
   SET duration    4e-8 # the simulated time, in seconds!
   SET output      1    # the output setup we want to use
   SET courant     0.98 # courant factor
END simulation    1
</code></pre>
<p>The simulation primitive is a kind of wrapper for the rest of the script. One script can have multiple simulation definitions, which will run sequentially<sup>[^2^]</sup>. The <code>solver</code> (which can be a conjugate gradient or algebraic multigrid<a href="http://wci.llnl.gov/codes/silo/">[AGMG]</a>) and <code>tolerance</code> fields are ignored if the <code>method</code> field is not set to the FEM (the only implicit approach among the three). The <code>courant</code> property sets a coefficient which multiplies the maximum time-step computed with spectral methods, therefore it should always be \(0 <\) <code>courant</code> \(< 1\) (0.98 was used in most of the results obtained in the the thesis).</p>

<p>Incidentally, the script we just described was the one used to obtain the outputs shown in Appendix A. Furthermore, for every primitive, the following general property applies:</p>

<p>Finally, we remark that, for any primitive <code>P</code> any number of objects of type <code>P</code> can be defined. If a <code>P</code> object instantiated with the same label <code>N</code> as a previously defined one (of the same type) overwrites the pre-existing definition.</p>
<p><em>[^1^]: This is the closest thing to a hack I feel I have gotten to in the development of this tool.</em></p>
<p><em>[^2^]: For now: the object of future work is also obviously parallelization of independent tasks.</em></p>
