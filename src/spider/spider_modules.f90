module sp_parameters

   implicit none

   integer, parameter :: nrp=256, ntp=256, nr1p=nrp-1, &
      nt1p=ntp-1, nr2p=nrp-2, nt2p=ntp-2
   integer, parameter :: neqp=nr2p*ntp + 1, neq1p=neqp + 1, &
      lp=9*neqp+ntp, nspp=3*neqp + 7*lp + 200000
   integer, parameter :: nloopp=128, nprobp=128
   integer, parameter :: njlim=1550, npfc0=200, nplim=400, nkp=njlim
   integer, parameter :: NiLIM=NJLIM - NpLIM, NCLIM=NJLIM, &
      NNLIM=NJLIM*NJLIM, NSP=90000
   integer, parameter :: nip=3, njp=3, ni2p=nip-2, nj2p = njp-2, &
      nbndp=2*(nip+njp), nbndp2=nbndp*2, nbndp4=nbndp2+4, nbndp6=nbndp4*6
   integer, parameter :: nblmp=600, neqp_blc=ni2p*nj2p, lp_blc=10
   integer, parameter :: nstep_p=5001
   integer, parameter :: nc_p=5, nf_p=3, ncf_p=nc_p+nf_p+1, &
      nnc_p=25, mmc_p=25
   integer, parameter :: nursp=4000, n_ursp=1000
   integer, parameter :: nekp=npfc0+nplim
   real*8, parameter :: pi=3.14159265359d0, amu0=0.4d0*pi, twopi = 2.d0*pi

end module sp_parameters

!--------------------
module rectgrid
  
  implicit none

  integer :: nr, nz
  double precision, allocatable, dimension(:) :: Rrect, Zrect
  
end module rectgrid

!--------------------
module iopath

   implicit none

   character(len=120) :: path

end module iopath

!--------------------
module comtim

   use sp_parameters, only: nstep_p

   implicit none

   real*8, dimension(nstep_p) :: torcur, rm_t, zm_t, rxp_t, zxp_t, &
      volpla_t, alp_t, betpol_t, bettor_t, zli3_t, psim_t, psib_t, time_t

end module comtim

!--------------------
module comrec

   use sp_parameters, only: nip, njp, nekp

   implicit none

   integer :: ni, nj, ni1, nj1, ni2, nj2, nbnd, &
      imax, jmax, ix1, jx1, ix2, jx2, ipr(nip, njp)
   real*8 :: xmax, xmin, ymax, ymin, ux0, ux1, ux2, up, um, xm, ym, &
             xx0, yx0, xx1, yx1, xx2, yx2, xx10, yx10, xx20, yx20
   real*8, dimension(nip) :: x, dx
   real*8, dimension(njp) :: y, dy
   real*8, dimension(nip, njp) :: u, ue, un
   real*4, dimension(nip, njp, nekp) :: zaindk
end module comrec

!--------------------
module parcur

   use sp_parameters, only: nc_p, nf_p, nnc_p, mmc_p

   implicit none

   integer :: n_cc, m_cc
   real*8 :: s_wght, rx_p, zx_p, dpsx_r, dpsx_z, psip_x, psip_a, &
             r_ax, z_ax, alph, Lfit, Lpre, dpsx_rz, dpsx_zz, dpsx_rr, &
             c_wght, psip_x2, rx2_p,zx2_p, dpsx2_r, dpsx2_z
   real*8, dimension(nc_p) :: rpre, zpre, psipre, curref, Gx_r, Gx_z, &
             Ginda, Gindx, Gx_rz, Gx_zz, Gx_rr, Gx2_r, Gx2_z
   real*8, dimension(nf_p) :: rfit, zfit, psifit, w_wght, d_wght
   real*8, dimension(nc_p, nf_p) :: Gindk, Gindp
   real*8, dimension(nnc_p) :: g_wght
   real*8, dimension(nnc_p, mmc_p) :: CxC, Cpir

end module parcur

!--------------------
module comevl

   use sp_parameters, only: npfc0, njlim

   implicit none

   integer :: NEQUI, NPFC, NLOC(0: npfc0)
   integer, dimension(npfc0) :: NEPFC, NDIV
   real*8 :: RM_EF, ZM_EF, TSTEP_EF, DT_EF, TIME_EF, RXPNT_EF, ZXPNT_EF
   real*8, dimension(npfc0) :: PFVOL1, PFVOL2, DPFVDT, PFCEQW, &
      PFCUR1, PFCUR2, PFCW1, PFCW2, PFCD1, PFCD2, WEPFC, NTURN
   real*8, dimension(npfc0, npfc0) :: PFRES
   real*8, dimension(njlim, njlim) :: PPIND, DPSIDJ

end module comevl

!--------------------
module comblc

   use sp_parameters, only: nkp, nbndp, nip, njp, nblmp, lp_blc, neqp_blc

   implicit none

   integer :: nblm, nctrl, numlim, ni, nj, ni1, nj1, ni2, nj2, &
     nbnd, nkin, nkout, neq, nnz, iter, itin, nrun, icont, &
     Nitl, Nitin, imax, jmax, ix1, jx1, ix2, jx2, iterbf, nnstpp

   integer, dimension(nkp) :: itok, jtok
   integer, dimension(nblmp) :: iblm, jblm
   integer, dimension(lp_blc) :: ja
   integer, dimension(neqp_blc+1) :: ia
   integer, dimension(nip, njp) :: ipr

   real*8 :: ublmax, alp, alpnew, tok, tokn, cnor, qcen, ucen, & 
     b0ax, r0ax, rmax, rmin, zmax, zmin, erru, errx, errm, eps, & 
     ux0, ux1, ux2, up, um, rm, zm, rx0, zx0, rx1, zx1, rx2, zx2, & 
     rx10, zx10, rx20, zx20, rm0, zm0, psi_bon, rcondzzz, & 
     rl, zl, clr, clz, f_cur

   real*8, dimension(lp_blc) :: a
   real*8, dimension(nblmp) :: rblm, zblm
   real*8, dimension(neqp_blc) :: right
   real*8, dimension(nbndp) :: dgdn
   real*8, dimension(nbndp, nbndp) :: binadg
   real*8, dimension(nkp, nbndp) :: pinadg
   real*8, dimension(nip) :: r, dr, dri, r12
   real*8, dimension(njp) :: z, dz, dzj
   real*8, dimension(nip, njp) :: u, ue, un, ui, g, curf, f, q
   real*8, dimension(nip, njp, 500) :: ue_k

end module comblc

!--------------------
module compol

   use sp_parameters, only: nrp, ntp, lp, neqp, neq1p

   implicit none

   integer :: nr, nt, nr1, nt1, nr2, nt2, iplas, iplas1, neq, nnz, neqpla, &
      iter, itin, nrun, icont, Nitl, Nitin, nitdel, nitbeg, ngav, iswtch
   integer, dimension(lp) :: ja
   integer, dimension(neq1p) :: ia
! Careful, names overlapping comblc
   real*8 :: tok, tokp, cnor, qcen, b0ax, r0ax, erru, errx, errm, eps, &
      psipla, psip, psim, rm, zm, psiax, psibon, psibon0, psi_eav, &
      toksfi, fvac, flucfm, tokff, tokpp, tokww
   real*8, dimension(lp) :: a, aop0, daop, app0, dapp
   real*8, dimension(neqp) :: right
   real*8, dimension(nrp) :: psia, dpsda, f, q, dwdpsi, dpdpsi, dfdpsi, flx_fi
   real*8, dimension(nrp, ntp) :: r, z, ro, tta, ronor, dlr, dlt, &
      sr, st, vol, s, sq1, sq2, sq3, sq4, vol1, vol2, vol3, vol4, &
      cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4, psi, psin, cur

end module compol

!--------------------
module compol_add

   use sp_parameters, only: nrp, ntp, neqp, nloopp, nkp, nprobp

   implicit none

   integer :: ich, itrmax, Nitmax, ixp1, jxp1, ixp2, jxp2, nk_out, &
      nlop_out, nprob_out
   integer, dimension(nkp) :: iprcon
   integer, dimension(nloopp) :: iprlop
   integer, dimension(nprobp) :: iprprob
   real*8 :: alp, alpnew, pscen, psix0, psix1, psxi2, rx0, zx0, rx1, zx1, &
      rx2, zx2, rl, zl, clr, clz, rolim, ron_max, ron_max_g(500), fpv
   real*8, dimension(neqp) :: zpro
   real*8, dimension(ntp) :: dgdn
   real*8, dimension(nrp, ntp) :: psie, psii, g, binadg, aex
   real*8, dimension(nkp, ntp) :: pinadg
   real*8, dimension(nloopp, ntp) :: adginl
   real*8, dimension(nprobp, ntp) :: adginr, adginz

end module compol_add

!--------------------
module keys
   implicit none
   integer :: kpr, kastr, kastr2, key_0st, key_prs, key_plc, kxwx, ksnf, &
      key_fixfree, kstep, key_out, key_fixbon
end module keys

!--------------------
module e_nels
    implicit none
    real*8 :: enels
end module e_nels

!--------------------
module curpl
    implicit none
    real*8 :: cur_pl
end module curpl

!--------------------
module ndmf
    implicit none
    integer :: n_dmf
end module ndmf

!--------------------
module pres
    use sp_parameters, only: nrp
    implicit none
    real*8, dimension(nrp) :: dPdFi
end module pres

!--------------------
module sigcd
    use sp_parameters, only: nrp
    implicit none
    real*8, dimension(nrp) :: C_sig, T_el, C_bts, C_driv
end module sigcd

!--------------------
module jb
    use sp_parameters, only: nrp
    implicit none
    real*8, dimension(nrp) :: Bj_av, curfi_av
end module jb

!--------------------
module tim
    implicit none
    real*8 :: dtim, ctim
end module tim
