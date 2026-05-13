
module screening_module
  use constants
  implicit none
  private

  public :: screening_potential
  public :: ion_sphere_radius
  public :: gamma_coupling
  public :: h0_fit
  public :: u_mean_field
  public :: H_mean_field

  real(kind=dbl), parameter :: e2_mev_fm   = 1.4397_dbl
  real(kind=dbl), parameter :: kb_mev_per_k = 8.617333262145d-11
  real(kind=dbl), parameter :: cm_to_fm    = 1.0d13
  real(kind=dbl), parameter :: tiny_val    = 1.0d-30

! ============================================================
! SCREENING PROJECT MODULE
!
! This module contains the screening physics added during the
! electron / plasma screening project.
!
! Two screening levels are included here:
!
!   (1) screening_potential(...)
!       Original toy linear screening model used for early
!       prototype barrier testing and visualization.
!
!   (2) H_mean_field(...)
!       Main mean-field plasma screening potential used in the
!       final screened barrier implementation.
!
! Supporting helpers added for the mean-field model:
!   - ion_sphere_radius
!   - gamma_coupling
!   - h0_fit
!   - u_mean_field
!
! The linear screening model is a prototype tool.
! The mean-field model is the main physics-based implementation.
! ============================================================

contains

  !=========================================================
  ! Original toy linear screening model
  !=========================================================
! ------------------------------------------------------------
! TOY LINEAR SCREENING MODEL
!
! Purpose:
!   Early prototype screening model used to test whether a
!   screening correction could be inserted into the barrier.
!
! Usage:
!   Primarily used for:
!     - early barrier plots
!     - code-development checks
!     - comparison against later mean-field screening
!
! Status:
!   Prototype only; not the final physics model.
! ------------------------------------------------------------

  real(kind=dbl) function screening_potential(r, V_screen1, V_screen2)
    implicit none
    real(kind=dbl), intent(in) :: r
    real(kind=dbl), intent(in) :: V_screen1
    real(kind=dbl), intent(in) :: V_screen2
    real(kind=dbl) :: r_min, r_max, m, b

    r_min = 0.1_dbl
    r_max = 6.0_dbl

    m = (V_screen2 - V_screen1) / (r_max - r_min)
    b = V_screen1 - m * r_min

    screening_potential = m * r + b
  end function screening_potential

  !=========================================================
  ! Ion-sphere radius a in fm for an OCP of identical nuclei
  !
  ! ni = rho / (A * mu)
  ! a  = (3 / 4 pi ni)^(1/3)
  !=========================================================
! ------------------------------------------------------------
! MEAN-FIELD SUPPORT FUNCTION
! Computes the ion-sphere radius a(rho) used in the
! one-component plasma screening model.
! ------------------------------------------------------------

  real(kind=dbl) function ion_sphere_radius(rho_g_cm3, A_int)
    implicit none
    real(kind=dbl), intent(in) :: rho_g_cm3
    integer, intent(in)        :: A_int

    real(kind=dbl) :: ni_cm3
    real(kind=dbl) :: a_cm
    real(kind=dbl) :: pi_local

    pi_local = acos(-1.0_dbl)

    if (rho_g_cm3 <= 0.0_dbl .or. A_int <= 0) then
      ion_sphere_radius = 0.0_dbl
      return
    end if

    ni_cm3 = rho_g_cm3 / (real(A_int, dbl) * amu)
    a_cm   = (3.0_dbl / (4.0_dbl * pi_local * ni_cm3))**(1.0_dbl / 3.0_dbl)

    ion_sphere_radius = a_cm * cm_to_fm
  end function ion_sphere_radius

  !=========================================================
  ! Coulomb coupling parameter Gamma = Z^2 e^2 / (a k_B T)
  !
  ! with a in fm, k_B T in MeV
  !=========================================================
! ------------------------------------------------------------
! MEAN-FIELD SUPPORT FUNCTION
! Computes the Coulomb coupling parameter Gamma(rho, T),
! which controls the strength of plasma screening.
! ------------------------------------------------------------

  real(kind=dbl) function gamma_coupling(rho_g_cm3, temp_k, A_int, Z_int)
    implicit none
    real(kind=dbl), intent(in) :: rho_g_cm3
    real(kind=dbl), intent(in) :: temp_k
    integer, intent(in)        :: A_int
    integer, intent(in)        :: Z_int

    real(kind=dbl) :: a_fm
    real(kind=dbl) :: kbt_mev

    if (temp_k <= 0.0_dbl .or. A_int <= 0 .or. Z_int <= 0) then
      gamma_coupling = 0.0_dbl
      return
    end if

    a_fm    = ion_sphere_radius(rho_g_cm3, A_int)
    kbt_mev = kb_mev_per_k * temp_k

    if (a_fm <= 0.0_dbl .or. kbt_mev <= 0.0_dbl) then
      gamma_coupling = 0.0_dbl
      return
    end if

    gamma_coupling = real(Z_int * Z_int, dbl) * e2_mev_fm / (a_fm * kbt_mev)
  end function gamma_coupling

  !=========================================================
  ! Chugunov et al. fit for h0(Gamma)
  !
  ! h0_fit(Gamma) =
  !   Gamma^(3/2) [ A1/sqrt(A2 + Gamma) + A3/(1 + Gamma) ]
  ! + B1 Gamma^2/(B2 + Gamma)
  ! + B3 Gamma^2/(B4 + Gamma^2)
  !=========================================================
! ------------------------------------------------------------
! CHUGUNOV / MEAN-FIELD FIT
! Implements the analytic fit used for h0(Gamma) in the
! mean-field screening model.
! ------------------------------------------------------------

  real(kind=dbl) function h0_fit(gamma)
    implicit none
    real(kind=dbl), intent(in) :: gamma

    real(kind=dbl), parameter :: A1 = 2.7822_dbl
    real(kind=dbl), parameter :: A2 = 98.34_dbl
    real(kind=dbl), parameter :: A3 = 1.4515_dbl
    real(kind=dbl), parameter :: B1 = -1.7476_dbl
    real(kind=dbl), parameter :: B2 = 66.07_dbl
    real(kind=dbl), parameter :: B3 = 1.12_dbl
    real(kind=dbl), parameter :: B4 = 65.0_dbl

    if (gamma <= 0.0_dbl) then
      h0_fit = 0.0_dbl
      return
    end if

    h0_fit = gamma**1.5_dbl * (A1 / sqrt(A2 + gamma) + A3 / (1.0_dbl + gamma)) &
           + B1 * gamma**2 / (B2 + gamma) &
           + B3 * gamma**2 / (B4 + gamma**2)
  end function h0_fit

  !=========================================================
  ! Chugunov et al. fit for u(x), x = r / a
  !
  ! u(x) = alpha0 * sqrt( numerator / denominator )
  !
  ! alpha0 = h0_fit(gamma) / gamma
  !
  ! numerator   = 1 - C4 x - (2 C1 / alpha0) x^2 + C3 x^4 + C2 x^8
  ! denominator = 1 + C2 alpha0^2 x^10
  !
  ! This is the fit that reproduces the small-r expansion,
  ! including alpha2 ~ -1/4 at strong coupling.
  !=========================================================
! ------------------------------------------------------------
! DIMENSIONLESS MEAN-FIELD PROFILE
! Computes u(x), where x = r/a, for the screened mean-field
! plasma potential shape.
! ------------------------------------------------------------

  real(kind=dbl) function u_mean_field(x, gamma)
    implicit none
    real(kind=dbl), intent(in) :: x
    real(kind=dbl), intent(in) :: gamma

    real(kind=dbl) :: alpha0
    real(kind=dbl) :: C1, C2, C3, C4
    real(kind=dbl) :: numerator, denominator, ratio

    if (gamma <= tiny_val) then
      u_mean_field = 0.0_dbl
      return
    end if

    alpha0 = h0_fit(gamma) / gamma

    if (abs(alpha0) <= tiny_val) then
      u_mean_field = 0.0_dbl
      return
    end if

    C1 = 0.25_dbl - 0.267_dbl * gamma**(-1.44_dbl)
    C2 = 0.05_dbl
    C3 = 0.084_dbl - 0.144_dbl * gamma**(-1.7_dbl)
    C4 = 0.434_dbl * gamma**(-1.2_dbl)

    numerator   = 1.0_dbl - C4 * x - (2.0_dbl * C1 / alpha0) * x**2 &
                + C3 * x**4 + C2 * x**8

    denominator = 1.0_dbl + C2 * alpha0**2 * x**10

    ratio = numerator / denominator

    if (ratio < 0.0_dbl) ratio = 0.0_dbl

    u_mean_field = alpha0 * sqrt(ratio)
  end function u_mean_field

  !=========================================================
  ! Mean-field screening potential
  !
  ! H(r) = k_B T * Gamma * u(r/a)
  !
  ! returns H in MeV
  !=========================================================
! ------------------------------------------------------------
! FINAL MEAN-FIELD SCREENING POTENTIAL
!
! Purpose:
!   Compute H(r), the mean-field plasma screening potential
!   used in the screened barrier:
!
!       U(r) = Vcoulomb(r) - H(r)
!
!   This is the main screening function used in the final
!   screened WKB and screened rate comparisons.
! ------------------------------------------------------------

  real(kind=dbl) function H_mean_field(r_fm, rho_g_cm3, temp_k, A_int, Z_int)
    implicit none
    real(kind=dbl), intent(in) :: r_fm
    real(kind=dbl), intent(in) :: rho_g_cm3
    real(kind=dbl), intent(in) :: temp_k
    integer, intent(in)        :: A_int
    integer, intent(in)        :: Z_int

    real(kind=dbl) :: a_fm
    real(kind=dbl) :: gamma
    real(kind=dbl) :: x
    real(kind=dbl) :: u
    real(kind=dbl) :: kbt_mev

    if (temp_k <= 0.0_dbl .or. rho_g_cm3 <= 0.0_dbl) then
      H_mean_field = 0.0_dbl
      return
    end if

    a_fm = ion_sphere_radius(rho_g_cm3, A_int)
    if (a_fm <= 0.0_dbl) then
      H_mean_field = 0.0_dbl
      return
    end if

    gamma   = gamma_coupling(rho_g_cm3, temp_k, A_int, Z_int)
    kbt_mev = kb_mev_per_k * temp_k
    x       = r_fm / a_fm
    u       = u_mean_field(x, gamma)

    H_mean_field = kbt_mev * gamma * u
  end function H_mean_field

end module screening_module

