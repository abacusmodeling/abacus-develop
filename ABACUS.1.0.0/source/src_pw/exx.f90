  subroutine vexx(lda, n, m, psi, hpsi)
  !band loop
       temppsic(:) = ( 0.D0, 0.D0 )
       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
   
        ! k point loop    
        do iqi=1,nqi
          ! G loop
             do ibnd=1,nbnd !for each band of psi
                
                !calculate rho in real space
                rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                !brings it to G-space
                CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
   
                vc(:) = ( 0.D0, 0.D0 )
                vc(nls(1:ngm)) = fac(1:ngm) * rhoc(nls(1:ngm))
                vc = vc * x_occupation(ibnd,ik) / nqs

                !brings back v in real space
                CALL cft3s( vc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 ) 

                !accumulates over bands and k points
                result(1:nrxxs)=result(1:nrxxs)+vc(1:nrxxs)*tempphic(1:nrxxs)
             end do
       end do
       !write(*,*) result(1:10)
       CALL mp_sum( result(1:nrxxs), inter_image_comm )
       !write(*,*) 'result is:  ',result(1), result(nrxxs), nrxxs, my_image_id
       !brings back result in G-space
       CALL cft3s( result, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
       !adds it to hpsi
       hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))

    deallocate (tempphic,temppsic, result, rhoc, vc, fac )
   
    call stop_clock ('vexx')

     end subroutine vexx

  function exxenergy ()
    ! This function is called to correct the deband value and have 
    ! the correct energy 
!    energy=0.d0
!    do ik=1,nks
!       current_k = ik
!       IF ( lsda ) current_spin = isk(ik)
 !      npw = ngk (ik)
 !      IF ( nks > 1 ) THEN
 !         READ( iunigk ) igk
 !         call get_buffer  (psi, nwordwfc, iunwfc, ik)
 !      ELSE
 !         psi(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
 !      END IF
 !      vxpsi(:,:) = (0.d0, 0.d0)
 !      call vexx(npwx,npw,nbnd,psi,vxpsi)
       do ibnd=1,nbnd
          energy = energy + &
                   wg(ibnd,ik) * ZDOTC(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
       end do
       if (gamma_only .and. gstart == 2) then
           do ibnd=1,nbnd
              energy = energy - &
                       0.5d0 * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd)
           end do
       end if
    end do

    if (gamma_only) energy = 2.d0 * energy


    call mp_sum( energy, intra_pool_comm )
    call mp_sum( energy, inter_pool_comm )

    exxenergy = energy

    call stop_clock ('exxenergy')
  end function exxenergy

  !-----------------------------------------------------------------------
  function exxenergy2()
  !-----------------------------------------------------------------------
    energy=0.d0

    tpiba2 = (fpi / 2.d0 / alat) **2

    allocate (tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    
    nqi=nqs/nimage

    IF ( nks > 1 ) REWIND( iunigk )
    do ikk=1,nks
       current_k = ikk
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          call get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF

       do jbnd=1, nbnd !for each band of psi (the k cycle is outside band)
          temppsic(:) = ( 0.D0, 0.D0 )
          temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)

          CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
       
          do iqi=1,nqi
             iq=iqi+nqi*my_image_id
             ikq  = index_xkq(current_k,iq)
             ik   = index_xk(ikq)
             isym = abs(index_sym(ikq))

             xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
             if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
             sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                      s(:,2,isym)*xk_cryst(2) + &
                      s(:,3,isym)*xk_cryst(3) 
             xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

             do ig=1,ngm
                q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                q(3)= xk(3,current_k) - xkq(3) + g(3,ig)
                qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )
                if (qq.gt.1.d-8) then
                   fac(ig)=e2*fpi/(tpiba2*qq + yukawa ) * grid_factor
                   if (gamma_only) fac(ig) = 2.d0 * fac(ig)
                   if (on_double_grid) fac(ig) = 0.d0
                else
                   fac(ig)= - exxdiv ! & ! or rather something else (see F.Gygi)
 !                            - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
                   if (yukawa.gt.1.d-8 .and. .not. x_gamma_extrapolation) then
                      fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
                   end if
                end if
             end do

                do ibnd=1,nbnd !for each band of psi
                   if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                   do ji=1, nrxxs
                      tempphic(ji)=exxbuff(ji,ikq,ibnd)
                   enddo
                   !calculate rho in real space
                   rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   !brings it to G-space
                   CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
                   vc = 0.D0
                   do ig=1,ngm
                      vc = vc + fac(ig) * rhoc(nls(ig)) * CONJG(rhoc(nls(ig)))
                   end do
                   vc = vc * omega * x_occupation(ibnd,ik) / nqs
                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                end do
          end do
       end do
    end do

    deallocate (tempphic, temppsic, rhoc, fac )

    call mp_sum( energy, inter_image_comm )
    call mp_sum( energy, intra_pool_comm )
    call mp_sum( energy, inter_pool_comm )

    exxenergy2 = energy

    call stop_clock ('exxen2')

  end function  exxenergy2

  function exx_divergence ()

     alpha  = 10.d0 * tpiba2 / ecutwfc

     dq1= 1.d0/DBLE(nq1)
     dq2= 1.d0/DBLE(nq2) 
     dq3= 1.d0/DBLE(nq3) 

     div = 0.d0
     do iq1=1,nq1
        do iq2=1,nq2
           do iq3=1,nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3 
              do ig=1,ngm
                 q(1)= xq(1) + g(1,ig)
                 q(2)= xq(2) + g(2,ig)
                 q(3)= xq(3) + g(3,ig)
                 qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) ) 
                 if (.not.on_double_grid) then
                    if ( qq.gt.1.d-8 .or. yukawa .gt. 1.d-8) then
                       div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                    else
                       div = div - alpha ! or maybe something else
                    end if
                 end if
              end do
           end do
        end do
     end do
     call mp_sum(  div, intra_pool_comm )

     div = div * e2 * fpi / tpiba2 / nqs
     alpha = alpha / tpiba2
     nqq = 100000
     dq = 5.0d0 / sqrt(alpha) /nqq
     aa = 0.d0
     do iq=0,  nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
        aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
     end do
     aa = aa * 8.d0 /fpi
     aa = aa + 1.d0/sqrt(alpha*0.25d0*fpi) 
     write (stdout,*) aa, 1.d0/sqrt(alpha*0.25d0*fpi)
    
     div = div - e2*omega * aa
!    div = div - e2*omega/sqrt(alpha*0.25d0*fpi)
     exx_divergence = div * nqs
     write (stdout,'(a,i4,a,3f12.4)') 'EXX divergence (',nq1,')= ', &
                                  div, alpha

     return
  end function exx_divergence 
  

end module exx
