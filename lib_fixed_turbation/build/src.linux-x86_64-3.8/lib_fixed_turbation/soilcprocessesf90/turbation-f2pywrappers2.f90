!     -*- f90 -*-
!     This file is autogenerated with f2py (version:1.24.4)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_soilcprocesses_turbation (il1, il2, ilg, ignd,&
     & iccp2, iccp1, cryodiffus, biodiffus, kterm, zbotw, isand, actlyr,&
     & spinfast, ipeatland, litter, soilc, litter_out, soilc_out, f2py_z&
     &botw_d0, f2py_zbotw_d1, f2py_isand_d0, f2py_isand_d1, f2py_actlyr_&
     &d0, f2py_ipeatland_d0)
      use soilcprocesses, only : turbation
      integer il1
      integer il2
      integer ilg
      integer ignd
      integer iccp2
      integer iccp1
      real cryodiffus
      real biodiffus
      real kterm
      integer spinfast
      integer f2py_zbotw_d0
      integer f2py_zbotw_d1
      integer f2py_isand_d0
      integer f2py_isand_d1
      integer f2py_actlyr_d0
      integer f2py_ipeatland_d0
      real zbotw(f2py_zbotw_d0,f2py_zbotw_d1)
      integer isand(f2py_isand_d0,f2py_isand_d1)
      real actlyr(f2py_actlyr_d0)
      integer ipeatland(f2py_ipeatland_d0)
      real litter(ilg,iccp2,ignd)
      real soilc(ilg,iccp2,ignd)
      real litter_out(ilg,iccp2,ignd)
      real soilc_out(ilg,iccp2,ignd)
      call turbation(il1, il2, ilg, ignd, iccp2, iccp1, cryodiffus, biod&
     &iffus, kterm, zbotw, isand, actlyr, spinfast, ipeatland, litter, s&
     &oilc, litter_out, soilc_out)
      end subroutine f2pywrap_soilcprocesses_turbation
      subroutine f2pywrap_soilcprocesses_tridiag (a, b, c, r, u, f2py_a_&
     &d0, f2py_b_d0, f2py_c_d0, f2py_r_d0, f2py_u_d0)
      use soilcprocesses, only : tridiag
      integer f2py_a_d0
      integer f2py_b_d0
      integer f2py_c_d0
      integer f2py_r_d0
      integer f2py_u_d0
      real a(f2py_a_d0)
      real b(f2py_b_d0)
      real c(f2py_c_d0)
      real r(f2py_r_d0)
      real u(f2py_u_d0)
      call tridiag(a, b, c, r, u)
      end subroutine f2pywrap_soilcprocesses_tridiag
      
      subroutine f2pyinitsoilcprocesses(f2pysetupfunc)
      interface 
      subroutine f2pywrap_soilcprocesses_turbation (il1, il2, ilg, ignd,&
     & iccp2, iccp1, cryodiffus, biodiffus, kterm, zbotw, isand, actlyr,&
     & spinfast, ipeatland, litter, soilc, litter_out, soilc_out, f2py_z&
     &botw_d0, f2py_zbotw_d1, f2py_isand_d0, f2py_isand_d1, f2py_actlyr_&
     &d0, f2py_ipeatland_d0)
      integer il1
      integer il2
      integer ilg
      integer ignd
      integer iccp2
      integer iccp1
      real cryodiffus
      real biodiffus
      real kterm
      integer spinfast
      integer f2py_zbotw_d0
      integer f2py_zbotw_d1
      integer f2py_isand_d0
      integer f2py_isand_d1
      integer f2py_actlyr_d0
      integer f2py_ipeatland_d0
      real zbotw(f2py_zbotw_d0,f2py_zbotw_d1)
      integer isand(f2py_isand_d0,f2py_isand_d1)
      real actlyr(f2py_actlyr_d0)
      integer ipeatland(f2py_ipeatland_d0)
      real litter(ilg,iccp2,ignd)
      real soilc(ilg,iccp2,ignd)
      real litter_out(ilg,iccp2,ignd)
      real soilc_out(ilg,iccp2,ignd)
      end subroutine f2pywrap_soilcprocesses_turbation 
      subroutine f2pywrap_soilcprocesses_tridiag (a, b, c, r, u, f2py_a_&
     &d0, f2py_b_d0, f2py_c_d0, f2py_r_d0, f2py_u_d0)
      integer f2py_a_d0
      integer f2py_b_d0
      integer f2py_c_d0
      integer f2py_r_d0
      integer f2py_u_d0
      real a(f2py_a_d0)
      real b(f2py_b_d0)
      real c(f2py_c_d0)
      real r(f2py_r_d0)
      real u(f2py_u_d0)
      end subroutine f2pywrap_soilcprocesses_tridiag
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_soilcprocesses_turbation,f2pywrap_soil&
     &cprocesses_tridiag)
      end subroutine f2pyinitsoilcprocesses


