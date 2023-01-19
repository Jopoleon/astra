      integer :: nr, nt, nr1, nt1, nr2, nt2, iplas, iplas1, 
     &   neq, nnz, neqpla, iter, itin, nrun, icont, Nitl, Nitin,
     &   nitdel, nitbeg, ngav, iswtch
      integer, dimension(lp) :: ja
      integer, dimension(neq1p) :: ia

      real*8 :: tok, tokp, cnor, qcen, b0ax, r0ax, 
     &   erru, errx, errm, eps, psipla, psip, psim, rm, zm, 
     &   psiax, psibon, psibon0, psi_eav, toksfi, fvac, flucfm,
     &   tokff, tokpp, tokww, dtim, ctim
      real*8, dimension(lp) :: a, aop0, daop, app0, dapp
      real*8, dimension(neqp) :: right
      real*8, dimension(nrp) :: psia, dpsda, f, q, dwdpsi, dpdpsi,
     &   dfdpsi, flx_fi
      real*8, dimension(ntp) :: teta, curs
      real*8, dimension(nrp, ntp) :: r, z, ro, ronor, 
     &   dlr, dlt, sr, st, vol, s, sq1, sq2, sq3, sq4, vol1, vol2, 
     &   vol3, vol4, cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4,
     &   psi, psin, cur

      common /com_cn/ tok, tokp, cnor, qcen, b0ax, r0ax
      common /com_ind/ nr, nt, nr1, nt1, nr2, nt2, iplas, iplas1
      common /com_rz/ r, z, ro, teta, ronor
      common/com_met/ dlr, dlt, sr, st, vol, s
      common/com_wr/ sq1, sq2, sq3, sq4, vol1, vol2, vol3, vol4, 
     &               cos1, cos2, cos3, cos4, sin1, sin2, sin3, sin4
      common /com_mat/ a, right, ia, ja, neq, nnz, neqpla
      common /com_itr/ erru, errx, errm, eps, 
     &                 iter, itin, nrun, icont, Nitl, Nitin
      common /com_pt/ psi, psin, psia, dpsda, psipla, psip, psim, 
     &                rm, zm, psiax, psibon, psibon0, psi_eav
      common /com_pla/ cur, f, q, dwdpsi, dpdpsi, dfdpsi, curs,
     &                 toksfi, fvac, flucfm
      common /com_ma0/ aop0, daop, app0, dapp, nitdel, nitbeg
      common /com_ctr/ ngav, iswtch
      common /com_tok/ tokff, tokpp, tokww
      common /com_flf/ flx_fi
      common /com_tim/ dtim, ctim
