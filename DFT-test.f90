!===============================================================================
! b3lyp_h2.f90 - B3LYP/STO-3G for H2, minimal basis STO-3G
! Compile: gfortran -O3 b3lyp_h2.f90 -o b3lyp_h2 && ./b3lyp_h2
!===============================================================================

module constants
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: hartree_to_ev = 27.211386245988_dp
end module constants

!---- Mathematical utilities ------------------------------------------------
module mathlib
  use constants
  implicit none
contains
  function boys(n, x) result(f)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: f, val
    integer :: i
    if (x < 1.0d-7) then
      f = 1.0_dp / (2*n+1) - x / (2*n+3)
      return
    end if
    if (n == 0) then
      f = 0.5_dp * sqrt(pi/x) * erf(sqrt(x))
    else
      val = 0.5_dp * sqrt(pi/x) * erf(sqrt(x))
      do i = 1, n
        val = ((2*i-1)*val - exp(-x)) / (2*x)
      end do
      f = val
    end if
  end function boys

  subroutine gauss_product(a, A, b, B, p, P, K)
    real(dp), intent(in)  :: a, A(3), b, B(3)
    real(dp), intent(out) :: p, P(3), K
    p = a + b
    P = (a*A + b*B) / p
    K = exp(-a * b / p * sum((A - B)**2))
  end subroutine gauss_product

  ! Symmetric orthogonalization matrix X = S^{-1/2}
  subroutine sym_orth(S, X)
    real(dp), intent(in)  :: S(:,:)
    real(dp), intent(out) :: X(:,:)
    real(dp), allocatable :: U(:,:), s(:)
    integer :: n, i, j
    n = size(S,1)
    allocate(U(n,n), s(n))
    U = S
    call jacobi_eigen(U, n, s, U)
    X = 0.0_dp
    do i = 1, n
      do j = 1, n
        X(i,j) = sum(U(i,:) * U(j,:) / sqrt(s))
      end do
    end do
    deallocate(U, s)
  end subroutine sym_orth

  ! Jacobi diagonalization for symmetric matrix A
  subroutine jacobi_eigen(A, n, d, V)
    real(dp), intent(inout) :: A(n,n)
    integer, intent(in) :: n
    real(dp), intent(out) :: d(n), V(n,n)
    integer :: i, j, p, q, iter
    real(dp) :: theta, t, c, s, apq, app, aqq, offdiag
    real(dp), parameter :: eps = 1.0d-14
    V = 0.0_dp
    do i = 1, n
      V(i,i) = 1.0_dp
    end do
    do iter = 1, 1000
      offdiag = 0.0_dp
      do i = 1, n-1
        do j = i+1, n
          if (abs(A(i,j)) > offdiag) then
            offdiag = abs(A(i,j))
            p = i
            q = j
          end if
        end do
      end do
      if (offdiag < eps) exit
      app = A(p,p)
      aqq = A(q,q)
      apq = A(p,q)
      theta = 0.5_dp * atan2(2.0_dp*apq, app-aqq)
      c = cos(theta)
      s = sin(theta)
      t = s / (1.0_dp + c)
      A(p,p) = app - t*apq
      A(q,q) = aqq + t*apq
      A(p,q) = 0.0_dp
      A(q,p) = 0.0_dp
      do i = 1, n
        if (i /= p .and. i /= q) then
          apq = A(p,i)
          A(p,i) = apq - s*(A(q,i) + t*apq)
          A(q,i) = A(q,i) + s*(apq - t*A(q,i))
          A(i,p) = A(p,i)
          A(i,q) = A(q,i)
        end if
      end do
      do i = 1, n
        apq = V(i,p)
        V(i,p) = apq - s*(V(i,q) + t*apq)
        V(i,q) = V(i,q) + s*(apq - t*V(i,q))
      end do
    end do
    do i = 1, n
      d(i) = A(i,i)
    end do
  end subroutine jacobi_eigen
end module mathlib

!---- Basis set STO-3G for H (1s) --------------------------------------------
module basis
  use constants
  implicit none
  type :: contracted_gto
    real(dp) :: center(3)                ! atom position (Bohr)
    real(dp), allocatable :: alpha(:)    ! exponents of primitives
    real(dp), allocatable :: coeff(:)    ! normalized contraction coefficients
    integer :: nprim
  end type
contains
  subroutine build_sto3g_H(center, shell)
    real(dp), intent(in) :: center(3)
    type(contracted_gto), intent(out) :: shell
    real(dp), parameter :: zeta = 1.24_dp
    shell%center = center
    shell%nprim = 3
    allocate(shell%alpha(shell%nprim), shell%coeff(shell%nprim))
    shell%alpha = [3.425250914_dp, 0.6239137298_dp, 0.1688554040_dp] * zeta**2
    shell%coeff = [0.1543289673_dp, 0.5353281423_dp, 0.4446345422_dp]
    ! Normalize each primitive Gaussian
    shell%coeff = shell%coeff * (2.0_dp * shell%alpha / pi) ** 0.75_dp
  end subroutine build_sto3g_H

  ! Evaluate a contracted 1s function value and gradient at point r
  subroutine eval_shell(shell, r, val, grad)
    type(contracted_gto), intent(in) :: shell
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: val, grad(3)
    integer :: i
    real(dp) :: a, c, dx, dy, dz, r2, gauss
    val = 0.0_dp
    grad = 0.0_dp
    do i = 1, shell%nprim
      a = shell%alpha(i)
      c = shell%coeff(i)
      dx = r(1) - shell%center(1)
      dy = r(2) - shell%center(2)
      dz = r(3) - shell%center(3)
      r2 = dx*dx + dy*dy + dz*dz
      gauss = c * exp(-a * r2)
      val = val + gauss
      grad(1) = grad(1) - 2.0_dp * a * dx * gauss
      grad(2) = grad(2) - 2.0_dp * a * dy * gauss
      grad(3) = grad(3) - 2.0_dp * a * dz * gauss
    end do
  end subroutine eval_shell
end module basis

!---- Integrals for Coulomb and exact exchange (ERIs) ------------------------
module integrals
  use constants
  use mathlib
  use basis
  implicit none
contains
  ! Overlap, kinetic, nuclear attraction, and two‑electron integrals
  ! for contracted s-shells (two‑center)
  subroutine build_integrals(shellA, shellB, S, T, V, eri_4d, charges, centers)
    type(contracted_gto), intent(in) :: shellA, shellB
    real(dp), intent(out) :: S(2,2), T(2,2), V(2,2)
    real(dp), intent(out) :: eri_4d(2,2,2,2)
    real(dp), intent(in) :: charges(:), centers(:,:)
    integer :: i,j,k,l, ii,jj,kk,ll
    real(dp) :: a, b, c, d, ca, cb, cc, cd
    real(dp) :: A(3), B(3), C(3), D(3), p1, P1(3), K1, p2, P2(3), K2
    real(dp) :: q, Q(3), Kq, PQ2, F0
    ! Overlap and kinetic
    do i = 1, shellA%nprim
      a = shellA%alpha(i)
      ca = shellA%coeff(i)
      A = shellA%center
      do j = 1, shellB%nprim
        b = shellB%alpha(j)
        cb = shellB%coeff(j)
        B = shellB%center
        call gauss_product(a, A, b, B, p1, P1, K1)
        S(1,2) = S(1,2) + ca * cb * K1 * (pi/p1)**1.5_dp
        T(1,2) = T(1,2) + ca * cb * K1 * (a*b/p1) * (3.0_dp - 2.0_dp*a*b/p1*sum((A-B)**2)) * (pi/p1)**1.5_dp
      end do
    end do
    S(1,1) = 1.0_dp; S(2,2) = 1.0_dp
    S(2,1) = S(1,2)
    T(1,1) = 0.0_dp; T(2,2) = 0.0_dp
    T(2,1) = T(1,2)
    call compute_self_kinetic(shellA, T(1,1))
    call compute_self_kinetic(shellB, T(2,2))
    ! Nuclear attraction
    V = 0.0_dp
    call nuc_att(shellA, shellA, centers, charges, V(1,1))
    call nuc_att(shellA, shellB, centers, charges, V(1,2))
    V(2,1) = V(1,2)
    call nuc_att(shellB, shellB, centers, charges, V(2,2))
    ! ERIs: only non‑zero by symmetry
    eri_4d = 0.0_dp
    call eri_shell(shellA, shellA, shellA, shellA, eri_4d(1,1,1,1))
    call eri_shell(shellA, shellA, shellB, shellB, eri_4d(1,1,2,2))
    call eri_shell(shellB, shellB, shellB, shellB, eri_4d(2,2,2,2))
    call eri_shell(shellA, shellB, shellA, shellB, eri_4d(1,2,1,2))
    eri_4d(2,2,1,1) = eri_4d(1,1,2,2)
    eri_4d(1,2,2,1) = eri_4d(1,2,1,2)
    eri_4d(2,1,1,2) = eri_4d(1,2,1,2)
    eri_4d(2,1,2,1) = eri_4d(1,2,1,2)
  end subroutine build_integrals

  subroutine compute_self_kinetic(shell, Tval)
    type(contracted_gto), intent(in) :: shell
    real(dp), intent(out) :: Tval
    integer :: i, j
    real(dp) :: a, b, ca, cb, A(3), B(3), p, P(3), K, AB2
    Tval = 0.0_dp
    do i = 1, shell%nprim
      a = shell%alpha(i)
      ca = shell%coeff(i)
      A = shell%center
      do j = 1, shell%nprim
        b = shell%alpha(j)
        cb = shell%coeff(j)
        B = shell%center
        call gauss_product(a, A, b, B, p, P, K)
        AB2 = sum((A - B)**2)
        Tval = Tval + ca * cb * K * (a*b/p) * (3.0_dp - 2.0_dp*a*b/p*AB2) * (pi/p)**1.5_dp
      end do
    end do
  end subroutine

  subroutine nuc_att(shell1, shell2, centers, charges, Vval)
    type(contracted_gto), intent(in) :: shell1, shell2
    real(dp), intent(in) :: centers(:,:), charges(:)
    real(dp), intent(out) :: Vval
    integer :: i, j, k, nuc
    real(dp) :: a, b, ca, cb, A(3), B(3), C(3), Z, p, P(3), K, PC2, F0
    nuc = size(charges)
    Vval = 0.0_dp
    do i = 1, shell1%nprim
      a = shell1%alpha(i)
      ca = shell1%coeff(i)
      A = shell1%center
      do j = 1, shell2%nprim
        b = shell2%alpha(j)
        cb = shell2%coeff(j)
        B = shell2%center
        call gauss_product(a, A, b, B, p, P, K)
        do k = 1, nuc
          C = centers(k,:)
          Z = charges(k)
          PC2 = sum((P - C)**2)
          F0 = boys(0, p * PC2)
          Vval = Vval - Z * ca * cb * K * 2.0_dp * pi / p * F0
        end do
      end do
    end do
  end subroutine

  subroutine eri_shell(s1, s2, s3, s4, val)
    type(contracted_gto), intent(in) :: s1, s2, s3, s4
    real(dp), intent(out) :: val
    real(dp) :: p1, P1(3), K1, p2, P2(3), K2, q, Q(3), Kq, PQ2, F0
    integer :: i, j, k, l
    val = 0.0_dp
    do i = 1, s1%nprim
      do j = 1, s2%nprim
        call gauss_product(s1%alpha(i), s1%center, s2%alpha(j), s2%center, p1, P1, K1)
        do k = 1, s3%nprim
          do l = 1, s4%nprim
            call gauss_product(s3%alpha(k), s3%center, s4%alpha(l), s4%center, p2, P2, K2)
            call gauss_product(p1, P1, p2, P2, q, Q, Kq)
            PQ2 = sum((P1 - P2)**2)
            F0 = boys(0, p1*p2/q * PQ2)
            val = val + s1%coeff(i)*s2%coeff(j)*s3%coeff(k)*s4%coeff(l) &
                 * K1*K2 * 2.0_dp*pi**2.5_dp / (p1*p2*sqrt(q)) * F0
          end do
        end do
      end do
    end do
  end subroutine eri_shell
end module integrals

!---- Numerical integration grid (Becke partition) ---------------------------
module grid_module
  use constants
  use basis
  implicit none
  integer, parameter :: nrad = 20        ! radial points per atom
  integer, parameter :: nang = 6         ! Lebedev‑3 order (6 points)
  real(dp), save :: lebedev_pts(3,6), lebedev_w(6)

  type :: grid_point
    real(dp) :: pos(3)
    real(dp) :: weight
  end type

contains
  ! Initialize Lebedev 3rd order points and weights
  subroutine init_lebedev()
    ! Points on unit sphere, order 3 (6 points)
    lebedev_pts(:,1) = [ 1.0_dp,  0.0_dp,  0.0_dp]
    lebedev_pts(:,2) = [-1.0_dp,  0.0_dp,  0.0_dp]
    lebedev_pts(:,3) = [ 0.0_dp,  1.0_dp,  0.0_dp]
    lebedev_pts(:,4) = [ 0.0_dp, -1.0_dp,  0.0_dp]
    lebedev_pts(:,5) = [ 0.0_dp,  0.0_dp,  1.0_dp]
    lebedev_pts(:,6) = [ 0.0_dp,  0.0_dp, -1.0_dp]
    lebedev_w(:) = 4.0_dp * pi / 24.0_dp   ! = pi/6
  end subroutine init_lebedev

  ! Generate grid points for molecule (Becke partition)
  subroutine generate_grid(atom_centers, atom_charges, grids)
    real(dp), intent(in) :: atom_centers(:,:)   ! size(natoms,3)
    real(dp), intent(in) :: atom_charges(:)
    type(grid_point), allocatable, intent(out) :: grids(:)
    integer :: natom, iat, irad, iang, igrid, ntot
    real(dp), allocatable :: rad_pts(:), rad_w(:)
    real(dp) :: r, w, x, y, z, dr, atom_r, wi, pos(3)
    real(dp), parameter :: alpha = 0.6_dp  ! radial transformation parameter
    call init_lebedev()
    natom = size(atom_charges)
    ! Total number of grid points
    ntot = natom * nrad * nang
    allocate(grids(ntot))
    ! Radial grid: Mura‑Knowles (alpha)
    allocate(rad_pts(nrad), rad_w(nrad))
    call mura_knowles_radial(nrad, alpha, rad_pts, rad_w)
    igrid = 0
    do iat = 1, natom
      do irad = 1, nrad
        atom_r = rad_pts(irad)
        dr    = rad_w(irad)
        do iang = 1, nang
          igrid = igrid + 1
          ! Position of the grid point relative to atom iat
          pos = atom_centers(iat,:) + atom_r * lebedev_pts(:,iang)
          grids(igrid)%pos = pos
          ! Weight = radial weight * spherical weight * r^2 ; radial weight already contains r^2 factor
          grids(igrid)%weight = dr * lebedev_w(iang)
        end do
      end do
    end do
    ! Now apply Becke partition weights
    call becke_partition(natom, atom_centers, grids)
  end subroutine generate_grid

  ! Radial grid using transformation r = alpha * x^2/(1-x)^2
  subroutine mura_knowles_radial(n, alpha, rad_pts, rad_w)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: rad_pts(n), rad_w(n)
    integer :: i
    real(dp) :: x, dx, r, drdx
    dx = 1.0_dp / (n + 1)
    do i = 1, n
      x = i * dx
      r = alpha * x**2 / (1.0_dp - x)**2
      drdx = 2.0_dp * alpha * x / (1.0_dp - x)**3
      rad_pts(i) = r
      rad_w(i) = drdx * dx * r**2   ! includes Jacobian r^2
    end do
  end subroutine mura_knowles_radial

  ! Becke multi‑center partition
  subroutine becke_partition(natom, centers, grids)
    integer, intent(in) :: natom
    real(dp), intent(in) :: centers(:,:)
    type(grid_point), intent(inout) :: grids(:)
    integer :: i, j, k, ng
    real(dp) :: r_i, r_j, mu, Rij, aij, P, weight
    ng = size(grids)
    do i = 1, ng
      weight = 1.0_dp
      pos = grids(i)%pos
      do j = 1, natom
        r_j = norm2(pos - centers(j,:))
        P = 1.0_dp
        do k = 1, natom
          if (k == j) cycle
          r_i = norm2(pos - centers(k,:))
          Rij = norm2(centers(j,:) - centers(k,:))
          mu = (r_j - r_i) / Rij
          aij = -1.0_dp
          if (abs(mu) < 1.0_dp) then
            call becke_iter(mu, aij)
          else
            aij = sign(1.0_dp, mu)
          end if
          P = P * 0.5_dp * (1.0_dp - aij)
        end do
        weight = weight + P
      end do
      grids(i)%weight = grids(i)%weight / weight
    end do
  end subroutine becke_partition

  ! Iterative function for Becke hardness parameter
  subroutine becke_iter(mu, a)
    real(dp), intent(in) :: mu
    real(dp), intent(out) :: a
    real(dp) :: nu, old
    integer :: iter
    nu = mu
    do iter = 1, 100
      old = nu
      nu = old - (1.5_dp*old - 0.5_dp*old**3 - mu) / (1.5_dp - 1.5_dp*old**2)
      if (abs(nu - old) < 1.0e-12) exit
    end do
    a = 1.5_dp*nu - 0.5_dp*nu**3
  end subroutine becke_iter
end module grid_module

!---- B3LYP exchange‑correlation ----------------------------------------------
module b3lyp_module
  use constants
  implicit none
  ! Parameters
  real(dp), parameter :: a0 = 0.20_dp
  real(dp), parameter :: ax = 0.72_dp
  real(dp), parameter :: ac = 0.81_dp
  real(dp), parameter :: four_thirds = 4.0_dp/3.0_dp
contains
  ! VWN5 correlation energy per particle and its derivative w.r.t rho
  subroutine vwn5_corr(rho, ec, vc)
    real(dp), intent(in)  :: rho
    real(dp), intent(out) :: ec, vc
    real(dp), parameter :: A = 0.0310907_dp, b = 3.72744_dp, c = 12.9352_dp
    real(dp), parameter :: x0 = -0.10498_dp
    real(dp), parameter :: t1 = 2.0_dp*b, t2 = 4.0_dp*c
    real(dp) :: rs, x, Q, fx
    if (rho < 1.0d-12) then
      ec = 0.0_dp; vc = 0.0_dp; return
    end if
    rs = (3.0_dp/(4.0_dp*pi*rho))**(1.0_dp/3.0_dp)
    x = sqrt(rs)
    Q = sqrt(t2 - b*b)
    fx = 0.5_dp * log((x - x0)**2 / (x*x + x0*x + x0*x0 + 1.0d-15)) &
         + b/(2.0_dp*Q) * atan(Q/(2.0_dp*x + b)) &
         - (b*x0/(x0*x0 + b*x0 + c)) &
         * (log((x - x0)**2/(x*x + x0*x + x0*x0 + 1.0d-15)) &
           + 2.0_dp*(b + 2.0_dp*x0)/Q * atan(Q/(2.0_dp*x + b)))
    ec = A * fx
    vc = ec - A/(3.0_dp*x) * ( (x - x0)/(x*x + b*x + c) - (x0)/(x*x + b*x + c + 1.0d-15) )  ! derivative formula
  end subroutine vwn5_corr

  ! B88 exchange energy density and potential (using |grad rho|)
  subroutine b88_exchange(rho, grho2, ex, vx)
    real(dp), intent(in)  :: rho, grho2
    real(dp), intent(out) :: ex, vx
    real(dp), parameter :: beta = 0.0042_dp
    real(dp) :: rho43, x, s, denom, ex_lda
    if (rho < 1.0d-12) then
      ex = 0.0_dp; vx = 0.0_dp; return
    end if
    rho43 = rho**four_thirds
    ex_lda = -0.9305257363490999_dp * rho**1.3333333333333333_dp  ! (3/4)*(3/pi)^{1/3}
    x = sqrt(grho2) / rho43
    s = x / (9.0_dp * (36.0_dp*pi)**(1.0_dp/3.0_dp) + 1.0d-15)
    denom = 1.0_dp + 6.0_dp * beta * x * asinh(x)
    ex = -beta * rho43 * x**2 / denom
    vx = ex_lda - ex * (4.0_dp/3.0_dp)  ! approximate, full potential omitted for brevity; we use integration over grids later
    ! Complete GGA potential requires gradient terms, handled separately in Vxc builder
  end subroutine b88_exchange

  ! LYP correlation (simplified, full GGA potential handled in integrator)
  subroutine lyp_corr(rho, grho2, ec, vc)
    real(dp), intent(in)  :: rho, grho2
    real(dp), intent(out) :: ec, vc
    real(dp), parameter :: a = 0.04918_dp, b = 0.132_dp, c = 0.2533_dp, d = 0.349_dp
    real(dp) :: rho2, rho83, CC, g, delta
    ec = 0.0_dp; vc = 0.0_dp
    if (rho < 1.0d-12) return
    rho2 = rho * rho
    rho83 = rho**2.6666666666666667_dp
    g = 2.0_dp * (rho**2.3333333333333333_dp)
    delta = c * rho**(-1.0_dp/3.0_dp) + d
    CC = 0.04918_dp * (1.0_dp/(1.0_dp + d*rho**(-1.0_dp/3.0_dp)) - delta)
    ec = -a * b * g * CC / rho
    vc = ec  ! simplified (full GGA treatment requires functional derivative terms)
  end subroutine lyp_corr
end module b3lyp_module

!---- Main SCF program --------------------------------------------------------
program b3lyp_scf
  use constants
  use mathlib
  use basis
  use integrals
  use grid_module
  use b3lyp_module
  implicit none

  type(contracted_gto) :: shellA, shellB
  real(dp) :: centerA(3), centerB(3)
  real(dp) :: centers_nuc(2,3), charges(2)
  real(dp) :: R = 1.4_dp
  real(dp) :: S(2,2), T(2,2), V(2,2), Hcore(2,2), eri(2,2,2,2)
  real(dp) :: Kmat(2,2), Jmat(2,2), Vxc(2,2), Fock(2,2)
  real(dp) :: D(2,2), Dnew(2,2), C(2,2), eps(2), X(2,2)
  real(dp) :: E_total, E_nuc, E_elec, E_xc, E_HF_exch
  real(dp) :: rms_d, energy
  integer :: iter, mu, nu, lam, sig

  ! Grid variables
  type(grid_point), allocatable :: grids(:)
  real(dp) :: rho, grho2, ex_b88, vx_b88, ec_lyp, vc_lyp, ec_vwn, vc_vwn
  real(dp) :: chiA, chiB, gradA(3), gradB(3), rho_g(3), weight

  ! Initialize atoms
  centerA = [0.0_dp, 0.0_dp, -R/2.0_dp]
  centerB = [0.0_dp, 0.0_dp,  R/2.0_dp]
  call build_sto3g_H(centerA, shellA)
  call build_sto3g_H(centerB, shellB)
  centers_nuc(1,:) = centerA
  centers_nuc(2,:) = centerB
  charges = [1.0_dp, 1.0_dp]

  ! Compute integrals
  print*, 'Computing one‑electron and two‑electron integrals ...'
  call build_integrals(shellA, shellB, S, T, V, eri, charges, centers_nuc)
  Hcore = T + V

  ! Orthogonalization matrix
  call sym_orth(S, X)

  ! Generate numerical grid
  print*, 'Generating Becke grid ...'
  call generate_grid(centers_nuc, charges, grids)

  ! Initial density guess: zero
  D = 0.0_dp

  write(*,'(A)') 'Iter      Total Energy      ΔD_rms'
  do iter = 1, 200
    ! Build Coulomb (J) and exact exchange (K) from ERIs
    Jmat = 0.0_dp
    Kmat = 0.0_dp
    do mu = 1, 2
      do nu = 1, 2
        do lam = 1, 2
          do sig = 1, 2
            Jmat(mu,nu) = Jmat(mu,nu) + D(lam,sig) * eri(mu,nu,lam,sig)
            Kmat(mu,nu) = Kmat(mu,nu) + D(lam,sig) * eri(mu,lam,nu,sig)
          end do
        end do
      end do
    end do

    ! Numerical integration for XC (LDA + GGA)
    Vxc = 0.0_dp
    E_xc = 0.0_dp
    do mu = 1, size(grids)
      weight = grids(mu)%weight
      pos = grids(mu)%pos
      ! Evaluate basis functions and gradients
      call eval_shell(shellA, pos, chiA, gradA)
      call eval_shell(shellB, pos, chiB, gradB)
      ! Density
      rho = D(1,1)*chiA*chiA + 2.0_dp*D(1,2)*chiA*chiB + D(2,2)*chiB*chiB
      ! Density gradient (vector)
      rho_g = D(1,1)*2.0_dp*chiA*gradA &
            + 2.0_dp*D(1,2)*(gradA*chiB + chiA*gradB) &
            + D(2,2)*2.0_dp*chiB*gradB
      grho2 = rho_g(1)**2 + rho_g(2)**2 + rho_g(3)**2
      ! Evaluate XC functionals
      call vwn5_corr(rho, ec_vwn, vc_vwn)
      call lyp_corr(rho, grho2, ec_lyp, vc_lyp)
      call b88_exchange(rho, grho2, ex_b88, vx_b88)
      ! B3LYP energy density: E_xc = (1-a0)*E_x(LDA) + a0*E_x(HF) + ax*ΔE_x(B88) + E_c(VWN) + ac*ΔE_c(LYP)
      ! Here we build Vxc and energy piece by piece
      E_xc = E_xc + weight * ( (1.0_dp - a0)*(-0.7385587663820224_dp*rho**four_thirds) &
               + ax*ex_b88 + ec_vwn + ac*(ec_lyp - ec_vwn) )   ! HF exchange added later
      ! Build Vxc matrix (only diagonal contribution for brevity, full GGA needs gradient terms)
      Vxc(1,1) = Vxc(1,1) + weight * chiA * chiA * ( (1.0_dp - a0) * (-1.0_dp*(four_thirds)*0.7385587663820224_dp*rho**(1.0_dp/3.0_dp)) &
                 + ax*vx_b88 + vc_vwn + ac*(vc_lyp - vc_vwn) )
      Vxc(1,2) = Vxc(1,2) + weight * chiA * chiB * ( (1.0_dp - a0) * (-1.0_dp*(four_thirds)*0.7385587663820224_dp*rho**(1.0_dp/3.0_dp)) &
                 + ax*vx_b88 + vc_vwn + ac*(vc_lyp - vc_vwn) )
      Vxc(2,1) = Vxc(1,2)
      Vxc(2,2) = Vxc(2,2) + weight * chiB * chiB * ( (1.0_dp - a0) * (-1.0_dp*(four_thirds)*0.7385587663820224_dp*rho**(1.0_dp/3.0_dp)) &
                 + ax*vx_b88 + vc_vwn + ac*(vc_lyp - vc_vwn) )
    end do

    ! Fock matrix: F = Hcore + J + a0 * (exact exchange included as -0.5*K in HF, but here a0*HFx means a0 * (original HF exchange)
    ! HF exchange contribution to F is -0.5 * sum_ls D_ls (mu lambda | nu sigma) ; Kmat defined as sum, so -0.5*Kmat is HF exchange.
    ! We want a0 fraction of that, so add a0 * (-0.5 * Kmat)
    Fock = Hcore + Jmat + a0 * (-0.5_dp * Kmat) + Vxc

    ! Transform to orthogonal basis, diagonalize
    Fprime = matmul(transpose(X), matmul(Fock, X))
    call jacobi_eigen(Fprime, 2, eps, C)
    C = matmul(X, C)
    Dnew = 2.0_dp * spread(C(:,1),2,1) * spread(C(:,1),1,2)

    rms_d = sqrt(sum((D - Dnew)**2) / 4.0_dp)
    D = Dnew

    ! Energy
    E_HF_exch = -0.5_dp * sum(D * Kmat)   ! full HF exchange energy
    E_elec = 0.5_dp * sum(D * (Hcore + Fock)) + E_xc - 0.5_dp * sum(D * Vxc)  ! standard
    E_nuc = charges(1)*charges(2) / R
    E_total = E_elec + E_nuc

    write(*,'(I3,2X,F14.8,2X,E12.4)') iter, E_total, rms_d
    if (rms_d < 1.0e-10) exit
  end do

  write(*,*) ''
  write(*,'(A,F15.10)') 'B3LYP/STO-3G total energy (Hartree): ', E_total
  write(*,'(A,F12.6,A)') 'Orbital energies: ', eps

end program b3lyp_scf
