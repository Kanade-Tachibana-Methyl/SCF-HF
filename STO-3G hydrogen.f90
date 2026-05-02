!===============================================================================
! rhf_h2.f90 - Restricted Hartree-Fock for H2 molecule, STO-3G basis
! Compile: gfortran -O3 rhf_h2.f90 -o rhf_h2 && ./rhf_h2
!===============================================================================

module const
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp
end module const

!-------------------------------------------------------------------------------
module mathlib
  use const
  implicit none
contains
  !---------------------------------------------------------------------------
  function boys(n, x) result(f)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp) :: f
    real(dp) :: val
    integer :: i
    if (x < 1.0d-7) then
      f = 1.0d0 / (2*n+1) - x / (2*n+3)
      return
    end if
    if (n == 0) then
      f = 0.5d0 * sqrt(pi/x) * erf(sqrt(x))
    else
      val = 0.5d0 * sqrt(pi/x) * erf(sqrt(x))  ! F0
      do i = 1, n
         val = ((2*i-1)*val - exp(-x)) / (2*x)
      end do
      f = val
    end if
  end function boys

  !---------------------------------------------------------------------------
  ! Gauss product: g_a(A) * g_b(B) = K * g_{p}(P)
  subroutine gauss_product(a, A, b, B, p, P, K)
    real(dp), intent(in)  :: a, A(3), b, B(3)
    real(dp), intent(out) :: p, P(3), K
    real(dp) :: AB2
    p = a + b
    P = (a*A + b*B) / p
    AB2 = sum((A - B)**2)
    K = exp(-a*b/p * AB2)
  end subroutine gauss_product

  !---------------------------------------------------------------------------
  ! Symmetric orthogonalisation matrix X = S^{-1/2}
  subroutine sym_orth(S, X)
    real(dp), intent(in)  :: S(:,:)
    real(dp), intent(out) :: X(:,:)
    real(dp), allocatable :: U(:,:), s(:), work(:)
    integer :: n, info, lwork
    n = size(S,1)
    allocate(U(n,n), s(n))
    U = S
    ! LAPACK call for eigen-decomposition (use internal Jacobi instead)
    call jacobi_eigen(U, n, s, U)      ! U now eigenvectors, s eigenvalues
    ! X = U * s^{-1/2} * U^T
    X = 0.0d0
    call dsyr2k('u', 'n', 1.0d0, U, spread(1.0d0/sqrt(s),2,n), 0.0d0, X)  ! Not a standard approach, use explicit loop
    deallocate(U, s)
    ! Explicit construction
    allocate(U(n,n), s(n))
    U = S
    call jacobi_eigen(U, n, s, U)
    X = 0.0d0
    block
      integer :: i, j
      do i = 1, n
        do j = 1, n
          X(i,j) = sum(U(i,:) * U(j,:) / sqrt(s))
        end do
      end do
    end block
    deallocate(U, s)
  end subroutine sym_orth

  !---------------------------------------------------------------------------
  ! Jacobi diagonalisation for symmetric matrix A -> eigenvalues d, eigenvectors V
  subroutine jacobi_eigen(A, n, d, V)
    real(dp), intent(inout) :: A(n,n)
    integer, intent(in) :: n
    real(dp), intent(out) :: d(n), V(n,n)
    integer :: i, j, p, q, iter
    real(dp) :: theta, t, c, s, tau, apq, app, aqq, offdiag
    real(dp), parameter :: eps = 1.0d-14
    ! Initialise V to identity
    V = 0.0d0
    do i = 1, n
      V(i,i) = 1.0d0
    end do
    do iter = 1, 1000
      ! Find maximum off-diagonal element
      offdiag = 0.0d0
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
      theta = 0.5d0 * atan2(2.0d0*apq, app-aqq)
      c = cos(theta)
      s = sin(theta)
      t = s / (1.0d0 + c)

      ! Update matrix A
      A(p,p) = app - t*apq
      A(q,q) = aqq + t*apq
      A(p,q) = 0.0d0
      A(q,p) = 0.0d0
      do i = 1, n
        if (i /= p .and. i /= q) then
          apq = A(p,i)
          A(p,i) = apq - s*(A(q,i) + t*apq)
          A(q,i) = A(q,i) + s*(apq - t*A(q,i))
          A(i,p) = A(p,i)
          A(i,q) = A(q,i)
        end if
      end do
      ! Update eigenvectors
      do i = 1, n
        apq = V(i,p)
        V(i,p) = apq - s*(V(i,q) + t*apq)
        V(i,q) = V(i,q) + s*(apq - t*V(i,q))
      end do
    end do
    ! Eigenvalues
    do i = 1, n
      d(i) = A(i,i)
    end do
  end subroutine jacobi_eigen

end module mathlib

!-------------------------------------------------------------------------------
module hf_integrals
  use const
  use mathlib
  implicit none
contains
  !---------------------------------------------------------------------------
  ! H 1s STO-3G shell: exponent and contraction coefficient arrays
  subroutine sto3g_h(center, alpha, coeff, nprim)
    real(dp), intent(in)  :: center(3)
    real(dp), allocatable, intent(out) :: alpha(:), coeff(:)
    integer, intent(out) :: nprim
    real(dp), parameter :: zeta = 1.24_dp
    real(dp) :: a(3), d(3)
    a = [3.425250914_dp, 0.6239137298_dp, 0.1688554040_dp] * zeta**2
    d = [0.1543289673_dp, 0.5353281423_dp, 0.4446345422_dp]
    nprim = 3
    allocate(alpha(nprim), coeff(nprim))
    alpha = a
    ! Normalise each primitive
    coeff = d * (2.0_dp * alpha / pi) ** (0.75_dp)
  end subroutine sto3g_h

  !---------------------------------------------------------------------------
  subroutine overlap_matrix(alpha1, coeff1, center1, &
       alpha2, coeff2, center2, S)
    real(dp), intent(in)  :: alpha1(:), coeff1(:), center1(3)
    real(dp), intent(in)  :: alpha2(:), coeff2(:), center2(3)
    real(dp), intent(out) :: S
    integer :: i, j
    real(dp) :: p, P(3), K
    S = 0.0_dp
    do i = 1, size(alpha1)
      do j = 1, size(alpha2)
         call gauss_product(alpha1(i), center1, alpha2(j), center2, p, P, K)
         S = S + coeff1(i)*coeff2(j) * K * (pi / p) ** 1.5_dp
      end do
    end do
  end subroutine overlap_matrix

  subroutine kinetic_matrix(alpha1, coeff1, center1, &
       alpha2, coeff2, center2, T)
    real(dp), intent(in)  :: alpha1(:), coeff1(:), center1(3)
    real(dp), intent(in)  :: alpha2(:), coeff2(:), center2(3)
    real(dp), intent(out) :: T
    integer :: i, j
    real(dp) :: p, P(3), K, AB2, term
    T = 0.0_dp
    do i = 1, size(alpha1)
      do j = 1, size(alpha2)
         call gauss_product(alpha1(i), center1, alpha2(j), center2, p, P, K)
         AB2 = sum((center1 - center2)**2)
         term = alpha1(i)*alpha2(j)/p * (3.0_dp - 2.0_dp*alpha1(i)*alpha2(j)/p * AB2)
         T = T + coeff1(i)*coeff2(j) * K * term * (pi / p) ** 1.5_dp
      end do
    end do
  end subroutine kinetic_matrix

  subroutine nuclear_attraction(alpha1, coeff1, center1, &
       alpha2, coeff2, center2, centers_nuc, charges, V)
    real(dp), intent(in)  :: alpha1(:), coeff1(:), center1(3)
    real(dp), intent(in)  :: alpha2(:), coeff2(:), center2(3)
    real(dp), intent(in)  :: centers_nuc(:,:), charges(:)
    real(dp), intent(out) :: V
    integer :: i, j, k
    real(dp) :: p, P(3), K, PC2, F0
    V = 0.0_dp
    do i = 1, size(alpha1)
      do j = 1, size(alpha2)
         call gauss_product(alpha1(i), center1, alpha2(j), center2, p, P, K)
         do k = 1, size(charges)
            PC2 = sum((P - centers_nuc(k,:))**2)
            F0 = boys(0, p * PC2)
            V = V - charges(k) * coeff1(i)*coeff2(j) * K * 2.0_dp*pi/p * F0
         end do
      end do
    end do
  end subroutine nuclear_attraction

  subroutine eri_shells(alpha1, coeff1, c1, &
       alpha2, coeff2, c2, alpha3, coeff3, c3, &
       alpha4, coeff4, c4, eri_val)
    real(dp), intent(in) :: alpha1(:), coeff1(:), c1(3)
    real(dp), intent(in) :: alpha2(:), coeff2(:), c2(3)
    real(dp), intent(in) :: alpha3(:), coeff3(:), c3(3)
    real(dp), intent(in) :: alpha4(:), coeff4(:), c4(3)
    real(dp), intent(out) :: eri_val
    integer :: i, j, k, l
    real(dp) :: p1, P1(3), K1, p2, P2(3), K2
    real(dp) :: q, Q(3), Kq, PQ2, F0
    eri_val = 0.0_dp
    do i = 1, size(alpha1)
      do j = 1, size(alpha2)
         call gauss_product(alpha1(i), c1, alpha2(j), c2, p1, P1, K1)
         do k = 1, size(alpha3)
            do l = 1, size(alpha4)
               call gauss_product(alpha3(k), c3, alpha4(l), c4, p2, P2, K2)
               call gauss_product(p1, P1, p2, P2, q, Q, Kq)
               PQ2 = sum((P1 - P2)**2)
               F0 = boys(0, p1*p2/q * PQ2)
               eri_val = eri_val + coeff1(i)*coeff2(j)*coeff3(k)*coeff4(l) &
                    * K1*K2 * 2.0_dp*pi**2.5_dp / (p1*p2*sqrt(q)) * F0
            end do
         end do
      end do
    end do
  end subroutine eri_shells

end module hf_integrals

!-------------------------------------------------------------------------------
program rhf_h2
  use const
  use mathlib
  use hf_integrals
  implicit none

  real(dp) :: centerA(3), centerB(3)
  real(dp), allocatable :: alphaA(:), coeffA(:), alphaB(:), coeffB(:)
  integer :: nprimA, nprimB
  real(dp) :: S(2,2), T(2,2), V(2,2), Hcore(2,2)
  real(dp) :: eri(2,2,2,2)
  real(dp) :: X(2,2), D(2,2), Dnew(2,2), Fock(2,2), G(2,2)
  real(dp) :: C(2,2), eps(2), Fprime(2,2)
  real(dp) :: E_nuc, E_elec, E_total, rms_diff
  integer :: iter, mu, nu, lam, sig
  real(dp) :: centers_nuc(2,3), charges(2)
  ! H2 bond length in Bohr
  real(dp), parameter :: R = 1.4_dp

  ! ---------- Define molecule ----------
  centerA = [0.0_dp, 0.0_dp, -R/2.0_dp]
  centerB = [0.0_dp, 0.0_dp,  R/2.0_dp]
  call sto3g_h(centerA, alphaA, coeffA, nprimA)
  call sto3g_h(centerB, alphaB, coeffB, nprimB)

  centers_nuc(1,:) = centerA
  centers_nuc(2,:) = centerB
  charges = [1.0_dp, 1.0_dp]

  ! ---------- Compute integrals ----------
  write(*,*) 'Computing integrals ...'
  call overlap_matrix(alphaA, coeffA, centerA, alphaA, coeffA, centerA, S(1,1))
  call overlap_matrix(alphaA, coeffA, centerA, alphaB, coeffB, centerB, S(1,2))
  S(2,1) = S(1,2)
  call overlap_matrix(alphaB, coeffB, centerB, alphaB, coeffB, centerB, S(2,2))

  call kinetic_matrix(alphaA, coeffA, centerA, alphaA, coeffA, centerA, T(1,1))
  call kinetic_matrix(alphaA, coeffA, centerA, alphaB, coeffB, centerB, T(1,2))
  T(2,1) = T(1,2)
  call kinetic_matrix(alphaB, coeffB, centerB, alphaB, coeffB, centerB, T(2,2))

  call nuclear_attraction(alphaA, coeffA, centerA, alphaA, coeffA, centerA, centers_nuc, charges, V(1,1))
  call nuclear_attraction(alphaA, coeffA, centerA, alphaB, coeffB, centerB, centers_nuc, charges, V(1,2))
  V(2,1) = V(1,2)
  call nuclear_attraction(alphaB, coeffB, centerB, alphaB, coeffB, centerB, centers_nuc, charges, V(2,2))

  Hcore = T + V

  call eri_shells(alphaA, coeffA, centerA, alphaA, coeffA, centerA, alphaA, coeffA, centerA, alphaA, coeffA, centerA, eri(1,1,1,1))
  call eri_shells(alphaA, coeffA, centerA, alphaA, coeffA, centerA, alphaB, coeffB, centerB, alphaB, coeffB, centerB, eri(1,1,2,2))
  call eri_shells(alphaB, coeffB, centerB, alphaB, coeffB, centerB, alphaB, coeffB, centerB, alphaB, coeffB, centerB, eri(2,2,2,2))
  call eri_shells(alphaA, coeffA, centerA, alphaB, coeffB, centerB, alphaA, coeffA, centerA, alphaB, coeffB, centerB, eri(1,2,1,2))
  eri(2,2,1,1) = eri(1,1,2,2)
  eri(1,2,2,1) = eri(1,2,1,2)
  eri(2,1,1,2) = eri(1,2,1,2)
  eri(2,1,2,1) = eri(1,2,1,2)
  eri(1,1,1,2)=0.d0; eri(1,1,2,1)=0.d0; eri(1,2,1,1)=0.d0; eri(2,1,1,1)=0.d0
  eri(2,2,2,1)=0.d0; eri(2,2,1,2)=0.d0; eri(2,1,2,2)=0.d0; eri(1,2,2,2)=0.d0
  ! Other symmetries remain zero.

  ! ---------- Orthogonalizer X = S^{-1/2} ----------
  call sym_orth(S, X)

  ! ---------- Initial guess ----------
  D = 0.0_dp

  ! ---------- SCF loop ----------
  write(*,*) 'Starting SCF ...'
  do iter = 1, 200
    ! Form G matrix
    G = 0.0_dp
    do mu = 1, 2
      do nu = 1, 2
        do lam = 1, 2
          do sig = 1, 2
            G(mu,nu) = G(mu,nu) + D(lam,sig) * (2.0_dp*eri(mu,nu,lam,sig) - eri(mu,lam,nu,sig))
          end do
        end do
      end do
    end do

    Fock = Hcore + G

    ! Transform to orthogonal basis
    Fprime = matmul(transpose(X), matmul(Fock, X))

    ! Diagonalize
    call jacobi_eigen(Fprime, 2, eps, C)

    ! Back-transform eigenvectors
    C = matmul(X, C)

    ! New density (1 occupied orbital, closed-shell)
    Dnew = 2.0_dp * spread(C(:,1),2,1) * spread(C(:,1),1,2)  ! outer product

    rms_diff = sqrt(sum((D - Dnew)**2) / 4.0_dp)
    D = Dnew

    ! Energy
    E_elec = 0.5_dp * sum(D * (Hcore + Fock))
    E_nuc = charges(1)*charges(2) / R
    E_total = E_elec + E_nuc

    write(*,'(I4, 2X, F15.10, 2X, E12.4)') iter, E_total, rms_diff

    if (rms_diff < 1.0d-10) exit
  end do

  write(*,*) ''
  write(*,'(A,F15.10)') 'Final total energy (Hartree): ', E_total
  write(*,'(A,F15.10)') 'Nuclear repulsion energy:    ', E_nuc
  write(*,'(A,F15.10)') 'Electronic energy:           ', E_elec
  write(*,'(A,2F12.6)') 'Orbital energies:            ', eps
end program rhf_h2
