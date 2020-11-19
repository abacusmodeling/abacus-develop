SUBROUTINE gradcorr( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
  IF ( .NOT. dft_is_gradient() ) RETURN

  igcc_is_lyp = (get_igcc() == 3 .or. get_igcc() == 7)
  !
  etxcgc = 0.D0
  vtxcgc = 0.D0
  !
  nspin0=nspin
  fac = 1.D0 / DBLE( nspin0 )
  !
  ALLOCATE(    h( 3, nrxx, nspin0) )
  ALLOCATE( grho( 3, nrxx, nspin0) )
  ALLOCATE( rhoout( nrxx, nspin0) )
  !
  ALLOCATE( rhogsum( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho(rho,rhoout,segni,nrxx)
     !
     ! ... bring starting rhoout to G-space
     !
     DO is = 1, nspin0
        !
        psic(:) = rhoout(:,is)
        !
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
        !
        rhogsum(:,is) = psic(nl(:))
        !
     END DO
  ELSE
     !
     rhoout(:,1:nspin0)  = rho(:,1:nspin0)
     rhogsum(:,1:nspin0) = rhog(:,1:nspin0)
     !
  ENDIF
  DO is = 1, nspin0
     !
     rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
     rhogsum(:,is) = fac * rhog_core(:) + rhogsum(:,is)
     !
     CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                   rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
     !
  END DO
  !
  DEALLOCATE( rhogsum )
  !
  IF ( nspin0 == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     !
     DO k = 1, nrxx
        !
        arho = ABS( rhoout(k,1) )
        !
        IF ( arho > epsr ) THEN
           !
           grho2(1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
           !
           IF ( grho2(1) > epsg ) THEN
              !
              segno = SIGN( 1.D0, rhoout(k,1) )
              !
              CALL gcxc( arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
              !
              ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
              !
              v(k,1) = v(k,1) + e2 * ( v1x + v1c )
              !
              ! ... h contains :
              !
              ! ...    D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
              !
              h(:,k,1) = e2 * ( v2x + v2c ) * grho(:,k,1)
              !
              vtxcgc = vtxcgc+e2*( v1x + v1c ) * ( rhoout(k,1) - rho_core(k) )
              etxcgc = etxcgc+e2*( sx + sc ) * segno
              !
           ELSE
              h(:,k,1)=0.D0
           END IF
           !
        ELSE
           !
           h(:,k,1) = 0.D0
           !
        END IF
        !
     END DO
     !
  ELSE
     !
     ! ... spin-polarised case
     !
!$omp parallel do private( rh, grho2, sx, v1xup, v1xdw, v2xup, v2xdw, rup, rdw, &
!$omp             grhoup, grhodw, grhoud, sc, v1cup, v1cdw, v2cup, v2cdw, v2cud, &
!$omp             zeta, grh2, v2c, grup, grdw  ), &
!$omp             reduction(+:etxcgc,vtxcgc)
     DO k = 1, nrxx
        !
        rh = rhoout(k,1) + rhoout(k,2)
        !
        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
        !
        CALL gcx_spin( rhoout(k,1), rhoout(k,2), grho2(1), &
                       grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
        !
        IF ( rh > epsr ) THEN
           !
           IF ( igcc_is_lyp ) THEN
              !
              rup = rhoout(k,1)
              rdw = rhoout(k,2)
              !
              grhoup = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
              grhodw = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
              !
              grhoud = grho(1,k,1) * grho(1,k,2) + &
                       grho(2,k,1) * grho(2,k,2) + &
                       grho(3,k,1) * grho(3,k,2)
              !
              CALL gcc_spin_more( rup, rdw, grhoup, grhodw, grhoud, &
                                  sc, v1cup, v1cdw, v2cup, v2cdw, v2cud )
              !
           ELSE
              !
              zeta = ( rhoout(k,1) - rhoout(k,2) ) / rh
              if (nspin.eq.4.and.domag) zeta=abs(zeta)*segni(k)
              !
              grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                     ( grho(2,k,1) + grho(2,k,2) )**2 + &
                     ( grho(3,k,1) + grho(3,k,2) )**2
              !
              CALL gcc_spin( rh, zeta, grh2, sc, v1cup, v1cdw, v2c )
              !
              v2cup = v2c
              v2cdw = v2c
              v2cud = v2c
              !
           END IF
           !
        ELSE
           !
           sc    = 0.D0
           v1cup = 0.D0
           v1cdw = 0.D0
           v2c   = 0.D0
           v2cup = 0.D0
           v2cdw = 0.D0
           v2cud = 0.D0
           !
        ENDIF
        !
        ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        v(k,1) = v(k,1) + e2 * ( v1xup + v1cup )
        v(k,2) = v(k,2) + e2 * ( v1xdw + v1cdw )
        !
        ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3
           !
           grup = grho(ipol,k,1)
           grdw = grho(ipol,k,2)
           h(ipol,k,1) = e2 * ( ( v2xup + v2cup ) * grup + v2cud * grdw )
           h(ipol,k,2) = e2 * ( ( v2xdw + v2cdw ) * grdw + v2cud * grup )
           !
        END DO
        !
        vtxcgc = vtxcgc + &
                 e2 * ( v1xup + v1cup ) * ( rhoout(k,1) - rho_core(k) * fac )
        vtxcgc = vtxcgc + &
                 e2 * ( v1xdw + v1cdw ) * ( rhoout(k,2) - rho_core(k) * fac )
        etxcgc = etxcgc + e2 * ( sx + sc )
        !
     END DO
!$omp end parallel do
     !
  END IF
  !
  DO is = 1, nspin0
     !
     rhoout(:,is) = rhoout(:,is) - fac * rho_core(:)
     !
  END DO
  !
  DEALLOCATE( grho )
  !
  ALLOCATE( dh( nrxx ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     !
     CALL grad_dot( nrx1, nrx2, nrx3, nr1, nr2, nr3, &
                    nrxx, h(1,1,is), ngm, g, nl, alat, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     vtxcgc = vtxcgc - SUM( dh(:) * rhoout(:,is) )
     !
  END DO
  !
  vtxc = vtxc + omega * vtxcgc / ( nr1 * nr2 * nr3 )
  etxc = etxc + omega * etxcgc / ( nr1 * nr2 * nr3 )

  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( rhoout )
  !
  RETURN
  !
END SUBROUTINE gradcorr
