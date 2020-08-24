
module global
implicit none

!!!!!!!!!! external libraries !!!!!!!!!!

include 'mpif.h'
include 'fftw_f77.i'

!!!!!!!!!! hardware parameters !!!!!!!!!!

! processor information
integer,parameter :: npr = 12, npc = 6 ! number of processors in x & z direction
integer,parameter :: np = npr * npc ! total number of processors

! grid information
integer,parameter :: ndx = 384, ndy = 864, ndz = 72, ndyr = 864 ! number of points in x,y,z directions, number of ypoints in input file
integer,parameter :: ndxh = ndx / 2, ndzh = ndz / 2 ! half of number of points in x,z directions
integer,parameter :: ndx2 = ndx*3/2, ndz2 = ndz*3/2 ! 3/2 of number of points in x,z directions for fft dealising

integer,parameter :: ndxnpr = ndx/npr, ndynpr = ndy/npr, ndznpr = ndz/npr, ndyrnpr = ndyr/npr ! number of the above parameters on each processor
integer,parameter :: ndxhnpr = ndxnpr/2
integer,parameter :: ndz2npr = ndznpr*3/2

integer,parameter :: ndxnpc = ndx/npc, ndynpc = ndy/npc, ndznpc = ndz/npc, ndyrnpc = ndyr/npc
integer,parameter :: ndxhnpc = ndxnpc/2

! buffer sizes in transposing mpi division
integer,parameter :: x1count = ndxhnpr * (ndynpc+1) * ndznpc
integer,parameter :: y1count = ndxhnpr * (ndynpc+1) * ndznpc
integer,parameter :: x2count = ndxhnpr * (ndynpc+1) * ndz2npr
integer,parameter :: z2count = ndxhnpr * (ndynpc+1) * ndz2npr
integer,parameter :: ilrecxnpy = (ndy+1) * 8, ilrecxnpz = 8

! max number of: scalars to be computed
integer,parameter :: ns = 0

!!!!!!!!!! MPI information !!!!!!!!!!

integer,save :: impi_myid ! id of the current processor
integer,save :: impi_errorinfo
integer,save :: coord(0:1) ! coordinate of the current processor in 2D processor grid
integer,save :: mpi_comm_row, mpi_comm_col ! sub communication domains for rows and columns in 2D processor grid

!!!!!!!!!! FFT plans !!!!!!!!!!

integer*8,save :: my_plan_x_f, my_plan_y_f, my_plan_z_f
integer*8,save :: my_plan_x_b, my_plan_y_b, my_plan_z_b
integer*8,save :: my_plan_y_fr

!!!!!!!!!! input parameters !!!!!!!!!!

real,save :: re, pr(ns) ! Reynolds number, Prandtl number
real,save :: arf, bat, dt ! streamwise minimum wave number, spanwise minimum wave number, time step length
real,save :: fr0, dp0 ! prescribed flow rate, precribed pressure gradient
integer,save :: nt, nsp, nsm, nsutp, nsprb ! number of total time steps to be computed, number of steps interval to write: instantaneous uvwp field, mid files for continuing computing, history file of XUTP and XFR, probe data
integer,save :: ifc, iff ! if continue computation, if fix flow rate or fix pressure gradient
integer,save :: ifscl, ifsgs, ifmfu, ifprb ! number of scalars to be computed, if use LES (not developed yet), if compute MFU (remove spanwise uniform eddies),  number of probes (slices to be high-time-resolved recorded)
integer,save :: ifscalarbound(ns), ifscalarsgs(ns)
integer,save :: j_probes(ndy+1) ! y positions of probes

!!!!!!!!!! mathematical parameters !!!!!!!!!!

! mathematical constants
real,parameter :: pai = 3.14159265359

! finite difference coefficients for 1st order derivative
real,save :: &
    firdera(1:ndy), &   ! lhs operator for x = ( d/dy ) f
    firderb(0:ndy), &   ! lhs operator for x = ( d/dy ) f
    firderc(0:ndy-1), & ! lhs operator for x = ( d/dy ) f
    firrhsc(0:ndy, 1:5) ! rhs operator for x = ( d/dy ) f
!real,save :: &
!    flt(2:ndy+1), &
!    ut (1:ndy+1), &
!    ct (1:ndy)
      
! finite difference coefficients for 2nd order derivative
real,save :: &
    secdera(1:ndy), &   ! rhs operator for ( d^2/dy^2 ) x = f
    secderb(0:ndy), &   ! rhs operator for ( d^2/dy^2 ) x = f
    secderc(0:ndy-1), & ! rhs operator for ( d^2/dy^2 ) x = f
    cofm1(0:ndxhnpr-1, 0:ndznpc-1), & ! rhs operator on lower boundary for ( d^2/dy^2 ) x = f, applying Neumann BC
    cofm2(0:ndxhnpr-1, 0:ndznpc-1)    ! rhs operator on upper boundary for ( d^2/dy^2 ) x = f, applying Neumann BC
real,save :: &
    ztv(3:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 - \gamma_0 * Re / dt ) x = f, applying Dirichlet BC
    gmv(2:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 - \gamma_0 * Re / dt ) x = f, applying Dirichlet BC
    afv(1:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 - \gamma_0 * Re / dt ) x = f, applying Dirichlet BC
    btv(1:ndy  , 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 - \gamma_0 * Re / dt ) x = f, applying Dirichlet BC
    qtv(1:ndy-1, 0:ndxhnpr-1, 0:ndznpc-1)    ! lhs operator for ( \nabla^2 - \gamma_0 * Re / dt ) x = f, applying Dirichlet BC
real,save :: &
    ztp(3:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Neumann BC
    gmp(2:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Neumann BC
    afp(1:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Neumann BC
    btp(1:ndy  , 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Neumann BC
    qtp(1:ndy-1, 0:ndxhnpr-1, 0:ndznpc-1)    ! lhs operator for ( \nabla^2 ) x = f, applying Neumann BC
real,save :: &
    ztp1(3:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Dirichlet BC
    gmp1(2:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Dirichlet BC
    afp1(1:ndy+1, 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Dirichlet BC
    btp1(1:ndy  , 0:ndxhnpr-1, 0:ndznpc-1), & ! lhs operator for ( \nabla^2 ) x = f, applying Dirichlet BC
    qtp1(1:ndy-1, 0:ndxhnpr-1, 0:ndznpc-1)    ! lhs operator for ( \nabla^2 ) x = f, applying Dirichlet BC
    
!!!!!!!!!! fields !!!!!!!!!!

! mesh grid
real,save :: &
    x (0:ndx-1), &
    y (0:ndy), &
    z (0:ndz-1)

! bulk parameters
real,save :: dp ! mean pressure gradient. only assigned in (0,0) core
real,save :: fr ! flow rate. only assigned in (0,0) core

! fields to be computed of current time step
complex,save :: &
    u (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! velocity U in spectral space
    v (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! velocity V in spectral space
    w (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! velocity W in spectral space
    p (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! total pressure PI ( = p + 0.5(u^2+v^2+w^2) ) in spectral space
    fx(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! non-linear term Fx in spectral space
    fy(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), & ! non-linear term Fy in spectral space
    fz(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)    ! non-linear term Fz in spectral space

! fields of past time steps
complex,save :: &
    u0 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    v0 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    w0 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    p0 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    u1 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    v1 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    w1 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    p1 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    u2 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    v2 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    w2 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    p2 (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    fx1(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    fy1(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    fz1(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    fx2(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    fy2(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    fz2(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

! pressure boundary condition
complex,save :: &
    pb0(0:ndxhnpr-1, 0:ndznpc-1, 2), &
    pb1(0:ndxhnpr-1, 0:ndznpc-1, 2), &
    pb2(0:ndxhnpr-1, 0:ndznpc-1, 2)

! scalar fields
complex,save :: &
    h (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    h0(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    h1(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    h2(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns)

! nonlinear terms of scalar fields
complex,save :: &
    fh (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    fh1(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns), &
    fh2(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1, ns)

! constant fields for divergence correction
complex,save :: & ! pressure correction
    pplus (0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    pminus(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)
complex,save :: & ! velocity correction
    ump(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    vmp(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    wmp(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    umm(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    vmm(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1), &
    wmm(0:ndxhnpr-1, 0:ndy, 0:ndznpc-1)

end module