# DG1D-CPP-development
Frozen state of the discontinuous Galerkin C++ code based on the Fortran code found in the DG1D-Fortran-test repository, which is a code sample given to me by Peter Diener to port to C++ at his request. 
I am not positive whether this is a working state or a non working state because I no longer have access to all the associated software needed to run the code. 
This was the primary component of my masters research at LSU from 2014-2017. The results of this research can be found on academia.edu.

At the time this software was frozen in Spring 2017, I had implemented a generalized object for solving wave equations and tested it with a wave equation in a vaccuum in flat spacetime with gaussian and sinusoidal initial conditions. 

Then I implemented a Schwarzschild metric from general relativity-- that's a spherically symmetric non rotating black hole with no charge. I used a scalar approximation to general relativity. In the absence of an orbiting black hole, I tested a perturbation to the "gravitational field" or scalar field at the location of the black hole. It produced the correct quasinormal mode ringdown that would be seen post merger in the merger of a smaller black hole into a supermassive black hole. 

Then I implemted the effective source model of the self-force. For a stellar mass black hole orbiting a supermassive black hole (an extreme mass ratio inspiral or EMRI), the central black hole curves spacetime, and the smaller black hole's mass interacts with it, curving its path away from a geodesic that a massless particle would follow and causing it to lose energy through gravitational waves and inspiral toward the supermassive black hole. I evolved and plotted the self force as a function of time for a stellar mass black hole held artifically fixed on a circular orbit in the scalar Schwarzchild model with an effective source. In other words, in nature the stellar mass black hole would have spiraled inward, but I computed the force the stellar mass "would have experienced" over one orbit if it had been on a circular orbit instead. 

To run this code, it is necessary to install the TNT and JAMA matrix
libraries, found at http://math.nist.gov/tnt/overview.html.

It is also necessary to install the libconfig library, found at
http://www.hyperrealm.com/libconfig. This can also be installed
through Macports using
sudo port install libconfig-hr
if you are running this on a Mac. 

Set the parameters you desire in params.cfg. If you wish to alter the
equation evolved by the program, edit the A and B matrices in
DiffEq.cpp. du/dt + A du/dx + Bu = 0. 

--------------------------------------------

main.cpp controls program flow. First, some preliminary setup occurs
in thegrid in Grid.cpp. The grid initializes the reference element
(ReferenceElement.cpp) to determine a nodes, a derivative matrix, and
a lift matrix for an element of a specific order with boundaries at -1
and 1. It then maps this to physical elements based on upper and lower
limits for the physical domain. It also calculates matrices associated
with the characteristic form of the differential equation in
preparation for use with the characteristic flux.

Ignoring details about output, the rest of the program works something
like this: main.cpp calls rk4lowstorage in Evolution.cpp with every
time step. That, in turn, evaluates the right hand side of the the
differential equation in DiffEq::characteristicFlux and DiffEq::RHS with
every sub-time step. This is used to advance the time step in the
Runga Kutta routine.

GridFunctions and VectorGridFunctions are storage
classes. GridFunctions store values at each point over the nodes of
the entire grid. VectorGridFunctions store values at each point over
the nodes of the grid for each variable in the differential equation.

The Grid stores information about spatial variables and coordinates
associated with the computational grid.

ConfigParams is involved in reading in parameters to create a global
variable, params, from params.cfg.

CharacteristicFlux stores information about the the differential
equation in characteristic form.

orbit currently works for a circular orbit.

namespaces stores information about a couple of namespaces. Include this in a file and
say "using <whatever the namespace is>" to acces the namespace.

Modes handles setting up and summing the spherical harmonics.

The ReferenceElement calculates the lift matrix (for the flux) and the
spatial derivative matrix needed in the right hand side. It does this
for a reference element, which can then be scaled to a physical
element based simply on the size of the physical element.

source_interface is legacy code. It is a wrapper to interfacing with
the effective source package, slightly modified before given to me.

TNT2, tnt_array1D_extn, and tnt_array2D_extn extend the TNT matrix
manipulation library.

WriteFile outputs files in specific formats.


-----------------------

To run this code it is necessary to install LibConfig (for reading in
parameters), the Template Numerical Toolkit (TNT), and JAMA, which is
associated with it. You will also need libgsl and read access to the 
scalar1deffective source repository on BitBucket.
