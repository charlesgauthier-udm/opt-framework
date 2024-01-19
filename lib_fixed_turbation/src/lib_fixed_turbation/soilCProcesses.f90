!> \file
!> Central module for all soil C processes involving movement of soil C up or down in soil column
module soilCProcesses

  ! J. Melton. May 30 2016

  implicit none

  ! Subroutines contained in this module:
  public  :: turbation
  public  :: tridiag

contains

  !=============================================================================================================

  !> \ingroup soilcprocesses_turbation
  !! Simulation of soil C movement due to turbation processes (presently only cryo).
  !!
  !! Modelled simply as a diffusion process.
  !!
  !> @author Joe Melton
  !! @{

  subroutine turbation (il1, il2, ilg, ignd, iccp2, iccp1, cryodiffus,biodiffus, kterm, zbotw, isand, &
                        actlyr, spinfast, ipeatland, & ! In
                        litter, soilC,litter_out,soilC_out) ! In/Out

    !use classicParams, only : icc, ilg, ignd, iccp2, iccp1, zero, tolrance, deltat, &
    !                            cryodiffus, biodiffus, kterm

    implicit none


    ! Arguments:
    integer, intent(in) :: il1       !< il1=1
    integer, intent(in) :: il2       !< il2=ilg
    integer, intent(in) :: ilg
    integer, intent(in) :: ignd
    integer, intent(in) :: iccp2
    integer, intent(in) :: iccp1
    real, parameter :: zero     = 1.0e-20                  !< Defintion of zerozbotw
    real, parameter  :: tolrance = 0.00025
    real, parameter :: deltat   = 1.0       !< CTEM's time step in days
    real, intent(in) :: cryodiffus
    real, intent(in) :: biodiffus
    real, intent(in) :: kterm
    integer, intent(in) :: spinfast  !< spinup factor for soil carbon whose default value is 1.
    !!  as this factor increases the soil c pool will come
    !! into equilibrium faster. Reasonable value for spinfast is between 5 and 10.
    integer, intent(in) ::  ipeatland(:)     !< Peatland flag. 0 non-peatlands, 1 = bog, 2 = fen.
    real, intent(in) :: zbotw(:,:)               !< Bottom of soil layers [m]
    integer, intent(in) :: isand(:,:)            !< flag for non-permeable layers
    real, intent(in) :: actlyr(:)                   !< active layer depth [m]

    real, intent(in) :: litter(ilg,iccp2,ignd)   !< litter mass for the pfts + bare [ \f$kg C/m^2\f$ ]
    real, intent(in) :: soilC(ilg,iccp2,ignd)   !< soil carbon mass for the pfts + bare [ \f$kg C/m^2\f$ ]
    real, intent(out) :: litter_out(ilg,iccp2,ignd)
    real, intent(out) :: soilC_out(ilg,iccp2,ignd)

    ! Local variables:
    real, allocatable, dimension(:) :: avect            !< vectors for tridiagonal solver, subdiagonal a
    real, allocatable, dimension(:) :: bvect            !< diagonal b
    real, allocatable, dimension(:) :: cvect            !< superdiagonal c
    real, allocatable, dimension(:) :: rvect_lt         !< Righthand side of Crank-Nicholson solver for litter
    real, allocatable, dimension(:) :: rvect_sc         !< Righthand side of Crank-Nicholson solver for soilc
    real, allocatable, dimension(:) :: soilcinter       !< soil carbon array for the tridiagonal solver[ \f$kg C/m^2\f$ ]
    real, allocatable, dimension(:) :: littinter        !< litter array for the tridiagonal solver [ \f$kg C/m^2\f$ ]
    real, allocatable, dimension(:) :: depthinter       !< soil depth at layer interfaces [m]
    real, allocatable, dimension(:) :: effectiveD       !< effective diffusion coefficient for that soil layer

    integer :: i, l, j, k                                    !< counters
    integer :: botlyr                                   !< index of deepest permeable soil layer
    integer :: turblyrbot                               !< index of deepest layer where turbation occurs
    real :: dzm                                         !< temp var
    real :: psoilc                                      !< total soil C mass before turbation [ \f$kg C/m^2\f$ ]
    real :: asoilc                                      !< total soil C mass after turbation [ \f$kg C/m^2\f$ ]
    real :: plit                                        !< total litter mass before turbation [ \f$kg C/m^2\f$ ]
    real :: alit                                        !< total litter mass after turbation [ \f$kg C/m^2\f$ ]
    real :: kactlyr                                     !< temp var
    real :: termr                                       !< temp var
    real :: amount                                      !< temp var
    real :: botthick                                    !< temp var
    real :: bturblyr                                    !< temp var
    real :: kbturblyr                                   !< temp var
    real :: soilCchange(ilg,iccp1,ignd)
    real :: litterchange(ilg,iccp1,ignd)
    real :: diffus                                      !< diffusion coefficient used (either cryodiffus or biodiffus)
    real :: initLitter(ilg,iccp2,ignd)                  !< Initial litter pool
    real :: initSoilC(ilg,iccp2,ignd)                   !< Initial soil C pool
    !----

    initLitter = litter
    initSoilC = soilC

    do i = il1,il2
      if (ipeatland(i) /= 1 .or. ipeatland(i) /= 2) then ! turbation only occurs in mineral soils (so not bogs/fens)

        !> First, find the bottom of the permeable soil column by looking at the isand flags.
        botlyr = 0
        do l = 1,ignd
          if (isand(i,l) == - 3 .or. isand(i,l) == - 4) exit
          botlyr = l
        end do

        !> Next if the grid cell has some permeable soil (i.e. is not bedrock or ice) then take what was found to
        !! be the bottom soil layer and
        if (botlyr > 0) then

          botlyr = botlyr + 2 !> increment two more for our boundary conditions.

          !> Next,allocate the vectors for the tridiagional solver.
          allocate(avect(botlyr))
          allocate(bvect(botlyr))
          allocate(cvect(botlyr))
          allocate(rvect_sc(botlyr))
          allocate(rvect_lt(botlyr))
          allocate(soilcinter(botlyr))
          allocate(littinter(botlyr))
          allocate(depthinter(botlyr))
          allocate(effectiveD(botlyr))

          do j = 1,iccp1     !> We assume that the LUC pools do not migrate from the first soil layer so
            !! we don't include them in our calculations (we don't calculate over position iccp2).

            !> At the start,store the size of the present pool for later comparison
            psoilc = sum(soilC(i,j,:))
            plit   = sum(litter(i,j,:))

            if (psoilc > zero) then !> and only do turbation for PFTs with some soil C. (assume soil C is always ~> lit).
              !> Next set up the tridiagional solver soil and litter C. In this the first interface is
              !! the soil surface (so soilC/litter is 0). The position at the bottom (botlyr) is
              !! where the soil meets the bedrock so it is also 0 there.

              ! The first interface is at the soil surface so the soilC/litter is 0.
              soilcinter(1) = 0.
              littinter(1)  = 0.
              depthinter(1) = 0.

              ! Put the soil C/litter into the soilcinter/littinter array where
              ! position 1 is the soil surface and position botlyr is the bedrock bottom.
              soilcinter(2:botlyr - 1) = soilC(i,j,1:botlyr - 2)
              littinter(2:botlyr - 1)  = litter(i,j,1:botlyr - 2)
              depthinter(2:botlyr - 1) = zbotw(i,1:botlyr - 2)

              ! The final interface is at the bottom of the soil column and the bedrock, so the soilC/litter is 0 there.
              soilcinter(botlyr) = 0.
              littinter(botlyr) = 0.

              ! Check for special case where the soil is permeable all the way to the bottom
              ! so botlyr = ignd. In that case botlyr is now greater than ignd.
              if (botlyr <= ignd) then
                botthick = zbotw(i,botlyr)
              else
                botthick = zbotw(i,ignd)
              end if

              depthinter(botlyr) = botthick

              !> Next, determine the effective diffusion coefficient for each soil layer. This follows \cite Koven2011-796
              !! It is assumed that there is a linear dependence on depth beyond the active layer
              !! with cryoturbation terminating about 3 times the active layer depth (so an active layer
              !! of 40 cm would have cryoturbation cease at 1.2 m.
              !! As bioturbation captures the influence of biota (worms, bugs, macrofauna, etc.) moving through the soil and
              !! thereby mixing it up, it only occurs at shallow depths down to around 30 cm. Thus, we reduce the diffusion
              !! coefficient below 10 cm (bturblyr) and assume no bioturbation occurrs below 30 cm (kturblyr), similar to how
              !! cryoturbation is treated.

              kactlyr = actlyr(i) * kterm
              bturblyr = 0.1
              kbturblyr = 0.3

              if (actlyr(i) <= 1.) then !> Cryoturbation is assumed to act only if the present active layer is shallower than one metre. FLAG
                ! Cryoturbation is dominant so set the diffusion coefficient to the cryoturbation values
                diffus = cryodiffus
              else  ! bioturbation is dominant
                diffus = biodiffus
              end if

              if (actlyr(i) <= 1.) then ! Cryoturbation is dominant
                do l = 1,botlyr
                  if (depthinter(l) <= actlyr(i)) then ! shallow so vigorous cryoturb.
                    effectiveD(l) = diffus * real(spinfast)
                  else if (depthinter(l) > actlyr(i) .and. depthinter(l) <= kactlyr) then ! linear reduction in diffusion coef.
                    effectiveD(l) = diffus * (1. - (depthinter(l) - actlyr(i)) &
                                    / ((kterm - 1.) * actlyr(i))) * real(spinfast)
                    turblyrbot = l  ! deepest layer in which turbation occurs
                  else ! too deep, no cryoturbation assumed
                    effectiveD(l) = 0.
                  end if
                end do  
              else
                do l = 1,botlyr
                  if (depthinter(l) <= bturblyr) then ! shallow so vigorous bioturb.
                    effectiveD(l) = diffus * real(spinfast)
                  else if (depthinter(l) > bturblyr .and. depthinter(l) <= kbturblyr) then ! linear reduction in diffusion coef.
                    effectiveD(l) = diffus * (1. - (depthinter(l) - bturblyr) &
                                    / (kbturblyr - bturblyr)) * real(spinfast)
                    turblyrbot = l  ! deepest layer in which turbation occurs
                  else ! too deep, no bioturbation assumed
                    effectiveD(l) = 0.
                  end if
                end do  
              end if

              !> Then set up the coefficients for the tridiagonal matrix solver.
              !! Derivation here: https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method#Example:_1D_diffusion

              ! upper boundary condition
              avect(1) = 0.0 ! not defined in this element
              bvect(1) = 1.0
              cvect(1) = 0.0
              rvect_sc(1) = 0.
              rvect_lt(1) = 0.

              ! soil layers
              do l = 2,botlyr - 1

                dzm = depthinter(l) - depthinter(l - 1)

                termr = effectiveD(l) * deltat / dzm ** 2
                avect(l) = - termr
                bvect(l) = 2. * (1. + termr)
                cvect(l) = - termr
                rvect_sc(l) = termr * soilcinter(l - 1) + 2. * (1. - termr) &
                              * soilcinter(l) +  termr * soilcinter(l + 1)
                rvect_lt(l) = termr * littinter(l - 1) + 2. * (1. - termr) &
                              * littinter(l) +  termr * littinter(l + 1)

              end do

              ! bottom boundary condition
              avect(botlyr) = 0.
              bvect(botlyr) = 1.
              cvect(botlyr) = 0. ! not defined in this element
              rvect_sc(botlyr) = 0
              rvect_lt(botlyr) = 0.

              !> Once set up, call tridiagonal solver for the soil C mass:
              call tridiag(avect, bvect, cvect, rvect_sc, soilcinter)

              !> Next call tridiagonal solver for the litter mass:
              call tridiag(avect, bvect, cvect, rvect_lt, littinter)
              
              !> Lastly, put the soil C/litter back into their arrays
              soilC_out(i,j,1:botlyr - 2) = soilcinter(2:botlyr - 1)
              litter_out(i,j,1:botlyr - 2) = littinter(2:botlyr - 1)

              !> We also ensure we have conservation of C by comparing the total C
              !! before and after the diffusion and distributing/removing
              !! any excess/deficit over the layers we operated on.
              asoilc = sum(soilC_out(i,j,:))
              alit   = sum(litter_out(i,j,:))

              !> (positive amount means we lost some C, negative we gained.)
              amount = psoilc - asoilc
              ! soilC(i,j,1:botlyr - 2) = soilC(i,j,1:botlyr - 2) + amount/real(botlyr - 2)
              soilC_out(i,j,1:turblyrbot) = soilC_out(i,j,1:turblyrbot) + amount/real(turblyrbot)

              amount = plit - alit
              ! litter(i,j,1:botlyr - 2) = litter(i,j,1:botlyr - 2) + amount/real(botlyr - 2)
              litter_out(i,j,1:turblyrbot) = litter_out(i,j,1:turblyrbot) + amount/real(turblyrbot)

            end if ! test if any soil c present
          end do ! icc

          !> Finish by deallocate the vectors for the tridiagional solver
          deallocate(avect)
          deallocate(bvect)
          deallocate(cvect)
          deallocate(rvect_sc)
          deallocate(rvect_lt)
          deallocate(soilcinter)
          deallocate(littinter)
          deallocate(depthinter)
          deallocate(effectiveD)

        end if ! must be some soil (botlyr > 0)
      end if ! mineral vs peat soils 
    end do ! i

  end subroutine turbation
  !! @}
  !> \ingroup soilcprocesses_tridiag
  !! Subroutine to solve triadiagonal system of equations
  !!
  !! Solves for a vector u of size N the tridiagonal linear set using given by equation 2.4.1 in
  !! Numerical recipes in Fortran 90 (\cite Press2007-bp) using a serial algorithm. Input vectors b (diagonal elements)
  !! and r (right-hand side) have size N, while a and c (off-diagonal elements) are not defined
  !! in the first and last elements, respectively.
  !!
  !! @author Joe Melton
  subroutine tridiag (a, b, c, r, u)

    ! Subroutine to solve triadiagonal system of equations
    ! Solves for a vector u of size N the tridiagonal linear set using given by equation 2.4.1 in
    ! Numerical recipes in Fortran 90 using a serial algorithm. Input vectors b (diagonal elements)
    ! and r (right-hand side) have size N, while a and c (off-diagonal elements) are not defined
    ! in the first and last elements, respectively.

    ! Based on Numerical Recipes in Fortran 90

    ! J. Melton, May 30 2016

    implicit none

    real, dimension(:), intent(in) :: a   !< Components of tridiag solver
    real, dimension(:), intent(in) :: b   !< Components of tridiag solver
    real, dimension(:), intent(in) :: c   !< Components of tridiag solver
    real, dimension(:), intent(in) :: r   !< Components of tridiag solver
    real, dimension(:), intent(out) :: u  !< Components of tridiag solver

    integer :: n
    integer :: k
    real :: bet
    real, dimension(size(b)) :: gam

    !----

    n = size(b)

    bet = b(1)

    u(1) = r(1) / bet

    ! decomposition and forward substitution

    do k = 2,n

      gam(k) = c(k - 1) / bet

      bet = b(k) - a(k) * gam(k)

      u(k) = (r(k) - a(k) * u(k - 1)) / bet

    end do

    ! backsubstitution

    do k = n - 1,1, - 1
      u(k) = u(k) - gam(k + 1) * u(k + 1)
    end do

  end subroutine tridiag
  !! @}
  !> \namespace soilcprocesses
  !! Central module for all soil C processes involving movement of soil C up or down in soil column.
  !!

end module soilCProcesses
