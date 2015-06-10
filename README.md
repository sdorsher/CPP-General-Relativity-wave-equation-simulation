# DG1D-CPP-development
Development of the discontinuous Galerkin C++ code based on Peter Diener's Fortran code found in the DG1D-Fortran-test repository.

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
differential equation in Grid::characteristicFlux and Grid::RHS with
every sub-time step. This is used to advance the time step in the
Runga Kutta routine.

GridFunctions and VectorGridFunctions are storage
classes. GridFunctions store values at each point over the nodes of
the entire grid. VectorGridFunctions store values at each point over
the nodes of the grid for each variable in the differential equation.
