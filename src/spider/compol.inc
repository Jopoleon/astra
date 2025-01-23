      integer :: nr, nt, nr1, nt1, nr2, nt2, iplas, iplas1, 
     &   neq, nnz, neqpla, iter, itin, nrun, icont, Nitl, Nitin,
     &   nitdel, nitbeg, ngav, iswtch
      integer, dimension(lp) :: ja
      integer, dimension(neq1p) :: ia

      real*8 :: tok, tokp, cnor, qcen, b0ax, r0ax, 
     &   erru, errx, errm, eps, psipla, psip, psim, rm, zm, 
     &   psiax, psibon, psibon0, psi_eav, toksfi, fvac, flucfm,
     &   tokff, tokpp, tokww
      real*8, dimension(lp) :: a, aop0, daop, app0, dapp
      real*8, dimension(neqp) :: right
      real*8, dimension(nrp) :: psia, dpsda, f, q, dwdpsi, dpdpsi,
     &   dfdpsi, flx_fi
      real*8, dimension(ntp) :: teta, curs
      real*8, dimension(nrp, ntp) :: r, z, ro, ronor, 
     &   dlr, dlt, sr, st, vol, s, sq1, sq2, sq3, sq4, vol1, vol2, 
     &   vol3, vol4, cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4,
     &   psi, psin, cur

      common /com_compol/ tok, tokp, cnor, qcen, b0ax, r0ax,
     & nr, nt, nr1, nt1, nr2, nt2, iplas, iplas1,
     & r, z, ro, teta, ronor,
     & dlr, dlt, sr, st, vol, s,
     & sq1, sq2, sq3, sq4, vol1, vol2, vol3, vol4, 
     & cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4,
     & a, right, ia, ja, neq, nnz, neqpla,
     & erru, errx, errm, eps, 
     & iter, itin, nrun, icont, Nitl, Nitin,
     & psi, psin, psia, dpsda, psipla, psip, psim, 
     & rm, zm, psiax, psibon, psibon0, psi_eav,
     & cur, f, q, dwdpsi, dpdpsi, dfdpsi, curs,
     & toksfi, fvac, flucfm,
     & aop0, daop, app0, dapp, nitdel, nitbeg,
     & ngav, iswtch,
     & tokff, tokpp, tokww,
     & flx_fi
