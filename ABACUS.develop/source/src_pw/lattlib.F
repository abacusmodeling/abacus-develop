      SUBROUTINE LATTYP(A1,A2,A3,ITYP,CELDIM,PRCHAN)
!   Routine LATTYP identifies the type of bravais lattice spanned by   *
!   the primitive lattice vectors A1, A2 and A3. The lattice will be   *
!   transformed to a 'standard crystallographic setting' and the cell  *
!   dimensions will be estimated using the conventions of routine      *
!   LATGEN. The relation between 'original' and 'transformed' lattice  *
!   vectors will be given in matrix form.                              *
!                                                                      *
!      A1, A2 and A3 are the three crystallographic vectors of the     *
!        primitive unit cell (given in cartesian coordinates x,y,z)    *
!        of the real space lattice given in units of bohrs ('a').      *
!      PRCHAN defines the unit number for printing output information. *
!                                                                      *
!   Output parameters:                                                 *
!   ------------------                                                 *
!      A1, A2 and A3 contain on output the three crystallographic      *
!      vectors being transformed to some 'standard setting'.           *
!      ITYP gives the lattice type (1-14) as keyed in routine LATGEN.  *
!      CELDIM gives the cell dimensions as defined in routine LATGEN.  *
!***********************************************************************
      INTEGER PRCHAN
      DIMENSION A1(3),A2(3),A3(3),RB(3,3),XB(3,3),CELDIM(6),YB(3)
      DIMENSION CELLDM(6),SA1(3),SA2(3),SA3(3),CDMINP(6)

! First test the original vectors (what do we apparently find?):
      ISTRGE=0
! Lengths of the axes and cosines between the axes:
      ABS1=XB(1,1)*XB(1,1)+XB(2,1)*XB(2,1)+XB(3,1)*XB(3,1)
      ABS2=XB(1,2)*XB(1,2)+XB(2,2)*XB(2,2)+XB(3,2)*XB(3,2)
      ABS3=XB(1,3)*XB(1,3)+XB(2,3)*XB(2,3)+XB(3,3)*XB(3,3)
      A2XA1=XB(1,1)*XB(1,2)+XB(2,1)*XB(2,2)+XB(3,1)*XB(3,2)
      A3XA1=XB(1,1)*XB(1,3)+XB(2,1)*XB(2,3)+XB(3,1)*XB(3,3)
      A3XA2=XB(1,2)*XB(1,3)+XB(2,2)*XB(2,3)+XB(3,2)*XB(3,3)
      ABS2M1=SQRT(ABS1+ABS2-2._q*A2XA1)
      ABS3M1=SQRT(ABS1+ABS3-2._q*A3XA1)
      ABS3M2=SQRT(ABS2+ABS3-2._q*A3XA2)
      ABS2P1=SQRT(ABS1+ABS2+2._q*A2XA1)
      ABS3P1=SQRT(ABS1+ABS3+2._q*A3XA1)
      ABS3P2=SQRT(ABS2+ABS3+2._q*A3XA2)
      A1P2M3=SQRT(ABS1+ABS2+ABS3+2._q*A2XA1-2._q*A3XA1-2._q*A3XA2)
      A2P3M1=SQRT(ABS1+ABS2+ABS3+2._q*A3XA2-2._q*A2XA1-2._q*A3XA1)
      A3P1M2=SQRT(ABS1+ABS2+ABS3+2._q*A3XA1-2._q*A3XA2-2._q*A2XA1)
      ABS123=A2XA1+A3XA1+A3XA2
      ABS1=SQRT(ABS1)
      ABS2=SQRT(ABS2)
      ABS3=SQRT(ABS3)
      A2XA1=A2XA1/ABS1/ABS2
      A3XA1=A3XA1/ABS1/ABS3
      A3XA2=A3XA2/ABS2/ABS3
      A1XA2=ABS(A2XA1)
      A1XA3=ABS(A3XA1)
      A2XA3=ABS(A3XA2)
! Warning extremly strange inputs (very large/small ratios 'c/a', 'b/a'
! or cosines very close to +/-1) could sometimes lead to some problems
! due to the finite precision of the machine ...). Here we test some
! parameters and if we find 'strange parameters', we (DO NOT) stop!
      IF (A1XA2>0.995_q) ISTRGE=ISTRGE+1
      IF (A2XA3>0.995_q) ISTRGE=ISTRGE+1
      IF (A1XA3>0.995_q) ISTRGE=ISTRGE+1
      IF ((ABS3/ABS1)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS2/ABS1)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS3/ABS2)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS1/ABS3)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS1/ABS2)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS2/ABS3)<0.001_q) ISTRGE=ISTRGE+1

! Following the identification of the crystal type, transformation of
! the basis to 'crystallographic standard settings' and estimate of the
! cell dimensions (using the conventions of routine LATGEN):
      CELLDM(1)=0._q
      CELLDM(2)=0._q
      CELLDM(3)=0._q
      CELLDM(4)=0._q
      CELLDM(5)=0._q
      CELLDM(6)=0._q
      IBRAV=15
! Crystal classes with A1*A2=A1*A3=A2*A3:
      IF ((ABS(A2XA1-A3XA1)<TINY).AND. &
     &                                  (ABS(A2XA1-A3XA2)<TINY)) THEN
! Crystal classes with A1*A2=A1*A3=A2*A3 and |A1|=|A2|=|A3|:
         IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &                                    (ABS(ABS2-ABS3)<TINY)) THEN
! Orthogonal axes means simple cubic: (IBRAV=1)
            IF (A1XA2<TINY) THEN
               IBRAV=1
               CELLDM(1)=ABS1
! Cosine(alpha) equal to -1/3 means body-centered cubic: (IBRAV=2)
            ELSE IF (ABS(A2XA1+(1._q/3._q))<TINY) THEN
               IBRAV=2
               CELLDM(1)=ABS1*2._q/SQRT(3._q)
! Cosine(alpha) equal to +1/2 means face-centered cubic: (IBRAV=3)
            ELSE IF (ABS(A2XA1-0.5_q)<TINY) THEN
               IBRAV=3
               CELLDM(1)=ABS1*SQRT(2._q)
! All other angles means trigonal (rhomboedric): (IBRAV=7)
            ELSE
               IBRAV=7
               CELLDM(1)=ABS1
               CELLDM(4)=A2XA1
            END IF
! Crystal classes with inequal length of lattice vectors but also with
! A1*A2=A1*A3=A2*A3:
! Orthogonal axes:
         ELSE IF (A1XA2<TINY) THEN
! Two axes with equal lengths means simple tetragonal: (IBRAV=5)
! Adjustment: 'c-axis' shall be the special axis.
            IF (ABS(ABS1-ABS2)<TINY) THEN
               IBRAV=5
               CELLDM(1)=ABS1
               CELLDM(3)=ABS3/ABS1
! No axes with equal lengths means simple orthorhombic (IBRAV=8):
! Adjustment: Sort the axis by increasing lengths:
            ELSE IF (((ABS3-ABS2)>TINY).AND. &
     &                                   ((ABS2-ABS1)>TINY)) THEN
               IBRAV=8
               CELLDM(1)=ABS1
               CELLDM(2)=ABS2/ABS1
               CELLDM(3)=ABS3/ABS1
            END IF
         END IF
! Crystal classes with A1*A3=A2*A3=/A1*A2:
      ELSE IF (ABS(A3XA1-A3XA2)<TINY) THEN
! One axis orthogonal with respect to the other two axes:
         IF (A1XA3<TINY) THEN
! Equal length of the two nonorthogonal axes:
            IF (ABS(ABS1-ABS2)<TINY) THEN
! Cosine(alpha) equal to -1/2 means hexagonal: (IBRAV=4)
! Adjustment: 'c-axis' shall be the special axis.
               IF (ABS(A2XA1+0.5_q)<TINY) THEN
                  IBRAV=4
                  CELLDM(1)=ABS1
                  CELLDM(3)=ABS3/ABS1
! Other angles mean base-centered orthorhombic: (IBRAV=11)
! Adjustment: Cosine between A1 and A2 shall be lower than zero, the
!             'c-axis' shall be the special axis.
               ELSE IF (A2XA1<(-1._q*TINY)) THEN
                  IBRAV=11
                  CELLDM(1)=ABS2P1
                  CELLDM(2)=ABS2M1/ABS2P1
                  CELLDM(3)=ABS3/ABS2P1
               END IF
! Different length of the two axes means simple monoclinic (IBRAV=12):
! Adjustment: Cosine(gamma) should be lower than zero, special axis
!             shall be the 'b-axis'(!!) and |A1|<|A3|:
            ELSE IF ((A2XA1<(-1._q*TINY)).AND. &
     &                                       ((ABS1-ABS2)>TINY)) THEN
               IBRAV=12
               IF ((PRCHAN>=0).AND.(PRCHAN<=99)) &
     &            WRITE(PRCHAN,*) ' Warning from LATTYP: '// &
     &                'Monoclinic adjustement (A1->A3, A2->A1, A3->A2)!'
               CELLDM(1)=ABS2
               CELLDM(2)=ABS3/ABS2
               CELLDM(3)=ABS1/ABS2
               CELLDM(5)=A2XA1
               SA3(1)=XB(1,1)
               SA3(2)=XB(2,1)
               SA3(3)=XB(3,1)
               SA1(1)=XB(1,2)
               SA1(2)=XB(2,2)
               SA1(3)=XB(3,2)
               SA2(1)=XB(1,3)
               SA2(2)=XB(2,3)
               SA2(3)=XB(3,3)
            END IF
! Arbitrary angles between the axes:
         ELSE
! |A1|=|A2|=|A3| means body-centered tetragonal (IBRAV=6):
! Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
! orthogonal to one another and (adjustment!): |A1+A3|=|A2+A3|<|A1+A2|
            IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &               (ABS(ABS2-ABS3)<TINY).AND. &
     &                    (ABS(ABS3P1-ABS3P2)<TINY).AND. &
     &                        (ABS(ABS2P1-ABS3P1)>TINY).AND. &
     &                             (ABS(ABS3*ABS3+ABS123)<TINY)) THEN
               IBRAV=6
               CELLDM(1)=ABS3P1
               CELLDM(3)=ABS2P1/ABS3P1
! |A1|=|A2|=/|A3| means base-centered monoclinic (IBRAV=13):
! Adjustement: The cosine between A1 and A3 as well as the cosine
!              between A2 and A3 should be lower than zero.
            ELSE IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &                          (A3XA1<(-1._q*TINY)).AND. &
     &                                       (A3XA2<(-1._q*TINY))) THEN
               IBRAV=13
               CELLDM(1)=ABS2P1
               CELLDM(2)=ABS2M1/ABS2P1
               CELLDM(3)=ABS3/ABS2P1
! Gilles Dewijs suggested the following change to lattlib.F
! since Gilles certainly did some testing I changed the statment
! according to his suggestions
!old
!              CELLDM(5)=A3XA2
!old
               CELLDM(5)=(A3XA2*ABS2*ABS3+A3XA1*ABS1*ABS3)/ABS2P1/ABS3      

            END IF
         END IF
! Crystal classes with A1*A2=/A1*A3=/A2*A3
      ELSE
! |A1|=|A2|=|A3| means body-centered orthorhombic (IBRAV=9):
! Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
! orthogonal to one another and (adjustment!): |A1+A2|>|A1+A3|>|A2+A3|
         IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &            (ABS(ABS2-ABS3)<TINY).AND. &
     &                  ((ABS3P1-ABS3P2)>TINY).AND. &
     &                         ((ABS2P1-ABS3P1)>TINY).AND. &
     &                             (ABS(ABS3*ABS3+ABS123)<TINY)) THEN
            IBRAV=9
            CELLDM(1)=ABS3P2
            CELLDM(2)=ABS3P1/ABS3P2
            CELLDM(3)=ABS2P1/ABS3P2
! |A1|=|A2-A3| and |A2|=|A1-A3| and |A3|=|A1-A2| means face-centered
! orthorhombic (IBRAV=10):
! Adjustment: |A1+A2-A3|>|A1+A3-A2|>|A2+A3-A1|
         ELSE IF ((ABS(ABS2M1-ABS3)<TINY).AND. &
     &                 (ABS(ABS3M1-ABS2)<TINY).AND. &
     &                       (ABS(ABS3M2-ABS1)<TINY).AND. &
     &                              ((A1P2M3-A3P1M2)>TINY).AND. &
     &                                   ((A3P1M2-A2P3M1)>TINY)) THEN
            IBRAV=10
            CELLDM(1)=A2P3M1
            CELLDM(2)=A3P1M2/A2P3M1
            CELLDM(3)=A1P2M3/A2P3M1
! Now there exists only one further possibility - triclinic (IBRAV=14):
! Adjustment: All three cosines shall be greater than zero and ordered:
         ELSE IF ((A2XA1>A3XA1).AND.(A3XA1>A3XA2).AND. &
     &                                             (A3XA2>TINY)) THEN
            IBRAV=14
            CELLDM(1)=ABS1
            CELLDM(2)=ABS2/ABS1
            CELLDM(3)=ABS3/ABS1
            CELLDM(4)=A3XA2
            CELLDM(5)=A3XA1
            CELLDM(6)=A2XA1
         END IF
      END IF
! Remember the result ... :
      IBRINP=IBRAV
      CDMINP(1)=CELLDM(1)
      CDMINP(2)=CELLDM(2)
      CDMINP(3)=CELLDM(3)
      CDMINP(4)=CELLDM(4)
      CDMINP(5)=CELLDM(5)
      CDMINP(6)=CELLDM(6)

! Warning extremly strange inputs (very large/small ratios 'c/a', 'b/a'
! or cosines very close to +/-1) could sometimes lead to some problems
! due to the finite precision of the machine ...). Here we test some
! parameters and if we find 'strange parameters' we (DO NOT) stop!
      ABS1=RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1)
      ABS2=RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2)
      ABS3=RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3)
      A2XA1=RB(1,1)*RB(1,2)+RB(2,1)*RB(2,2)+RB(3,1)*RB(3,2)
      A3XA1=RB(1,1)*RB(1,3)+RB(2,1)*RB(2,3)+RB(3,1)*RB(3,3)
      A3XA2=RB(1,2)*RB(1,3)+RB(2,2)*RB(2,3)+RB(3,2)*RB(3,3)
      ABS1=SQRT(ABS1)
      ABS2=SQRT(ABS2)
      ABS3=SQRT(ABS3)
      A1XA2=ABS(A2XA1/ABS1/ABS2)
      A1XA3=ABS(A3XA1/ABS1/ABS3)
      A2XA3=ABS(A3XA2/ABS2/ABS3)
      IF (A1XA2>0.99999_q) ISTRGE=ISTRGE+1
      IF (A2XA3>0.99999_q) ISTRGE=ISTRGE+1
      IF (A1XA3>0.99999_q) ISTRGE=ISTRGE+1
      IF ((ABS3/ABS1)<0.00001_q) ISTRGE=ISTRGE+1
      IF ((ABS2/ABS1)<0.00001_q) ISTRGE=ISTRGE+1
      IF ((ABS3/ABS2)<0.00001_q) ISTRGE=ISTRGE+1
      IF ((ABS1/ABS3)<0.00001_q) ISTRGE=ISTRGE+1
      IF ((ABS1/ABS2)<0.00001_q) ISTRGE=ISTRGE+1
      IF ((ABS2/ABS3)<0.00001_q) ISTRGE=ISTRGE+1

! Search a primitive basis with shortest lattice vectors:

      ICOUNT=0
      GOTO 199

! The algorithm is the following: Take first one vector and subtract or
! add the two other vectors until you have found the shortest vector.
! Do the same with the other two vectors ... . The iterative method
! used here garantuees that the cell volume will be kept constant and
! all cell vectors occuring during the iterations span primitive cells.
! Two parameters control the flow: ILOOP checks every loop whether a
! vector could be found with smaller or with equal length compared to
! the starting vector by subtraction of an other vector. If ILOOP is not
! set you cannot lower the vector by subtraction of an other vector so
! you have to try it with addition of the vectors. Therefore every loop
! has the following structure: First try it with Ai-ICOUNT*Aj. If the
! shortest vector is found for ICOUNT=0 then try to take Ai+ICOUNT*Aj.
! The loops break if the vector length starts to increase and the vector
! Ai will be set equal to the shortest vector (Ai <-- Ai +/- ICOUNT*Aj).
! IFLAG will be set equal to one whenever a shorter vector was found
! during the total iteration step (6 loops!). If IFLAG is not set we
! found a basis with shortest lattice vectors. If not we should start
! a new iteration ... .

! First try A1-ICOUNT*A2 (loop 1a):
  100 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,1)=RB(1,1)-RB(1,2)
      RB(2,1)=RB(2,1)-RB(2,2)
      RB(3,1)=RB(3,1)-RB(3,2)
      ABSR=SQRT(RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1))
! Breakpoint of loop 1a:
      IF (ABSR>(ABS1+TINY)) GOTO 102
! Found some shorter vector ? Then set IFLAG (and store the length):
      IF (ABSR<(ABS1-TINY)) THEN
         IFLAG=1
         ABS1=ABSR
      END IF
! Vector was not longer than the starting vector. Set ILOOP ... :
      ILOOP=1
      GOTO 100
! After leaving loop 1a we have subtracted A2 one times too often ... :
  102 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,1)=RB(1,1)+RB(1,2)
      RB(2,1)=RB(2,1)+RB(2,2)
      RB(3,1)=RB(3,1)+RB(3,2)
! Nothing found ? Try A1+ICOUNT*A2 (loop 1b):
      IF (ILOOP==0) THEN
  103    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,1)=RB(1,1)+RB(1,2)
         RB(2,1)=RB(2,1)+RB(2,2)
         RB(3,1)=RB(3,1)+RB(3,2)
         ABSR=SQRT(RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1))
! Breakpoint of loop 1b:
         IF (ABSR>(ABS1+TINY)) GOTO 105
! Found some shorter vector ? Then set IFLAG (and store the length):
         IF (ABSR<(ABS1-TINY)) THEN
            IFLAG=1
            ABS1=ABSR
         END IF
         GOTO 103
! After leaving loop 1b we have added A2 one times too often ... :
  105    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,1)=RB(1,1)-RB(1,2)
         RB(2,1)=RB(2,1)-RB(2,2)
         RB(3,1)=RB(3,1)-RB(3,2)
      END IF
! Reset ILOOP before starting the next loop ... :
      ILOOP=0
! Try A1-ICOUNT*A3 (loop 2a):
  110 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,1)=RB(1,1)-RB(1,3)
      RB(2,1)=RB(2,1)-RB(2,3)
      RB(3,1)=RB(3,1)-RB(3,3)
      ABSR=SQRT(RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1))
      IF (ABSR>(ABS1+TINY)) GOTO 112
      IF (ABSR<(ABS1-TINY)) THEN
         IFLAG=1
         ABS1=ABSR
      END IF
      ILOOP=1
      GOTO 110
  112 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,1)=RB(1,1)+RB(1,3)
      RB(2,1)=RB(2,1)+RB(2,3)
      RB(3,1)=RB(3,1)+RB(3,3)
! Nothing found ? Try A1+ICOUNT*A3 (loop 2b):
      IF (ILOOP==0) THEN
  113    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,1)=RB(1,1)+RB(1,3)
         RB(2,1)=RB(2,1)+RB(2,3)
         RB(3,1)=RB(3,1)+RB(3,3)
         ABSR=SQRT(RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1))
         IF (ABSR>(ABS1+TINY)) GOTO 115
         IF (ABSR<(ABS1-TINY)) THEN
            IFLAG=1
            ABS1=ABSR
         END IF
         GOTO 113
  115    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,1)=RB(1,1)-RB(1,3)
         RB(2,1)=RB(2,1)-RB(2,3)
         RB(3,1)=RB(3,1)-RB(3,3)
      END IF
      ILOOP=0
! Try A2-ICOUNT*A1 (loop 3a):
  120 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,2)=RB(1,2)-RB(1,1)
      RB(2,2)=RB(2,2)-RB(2,1)
      RB(3,2)=RB(3,2)-RB(3,1)
      ABSR=SQRT(RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2))
      IF (ABSR>(ABS2+TINY)) GOTO 122
      IF (ABSR<(ABS2-TINY)) THEN
         IFLAG=1
         ABS2=ABSR
      END IF
      ILOOP=1
      GOTO 120
  122 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,2)=RB(1,2)+RB(1,1)
      RB(2,2)=RB(2,2)+RB(2,1)
      RB(3,2)=RB(3,2)+RB(3,1)
! Nothing found ? Try A2+ICOUNT*A1 (loop 3b):
      IF (ILOOP==0) THEN
  123    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,2)=RB(1,2)+RB(1,1)
         RB(2,2)=RB(2,2)+RB(2,1)
         RB(3,2)=RB(3,2)+RB(3,1)
         ABSR=SQRT(RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2))
         IF (ABSR>(ABS2+TINY)) GOTO 125
         IF (ABSR<(ABS2-TINY)) THEN
            IFLAG=1
            ABS2=ABSR
         END IF
         GOTO 123
  125    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,2)=RB(1,2)-RB(1,1)
         RB(2,2)=RB(2,2)-RB(2,1)
         RB(3,2)=RB(3,2)-RB(3,1)
      END IF
      ILOOP=0
! Try A2-ICOUNT*A3 (loop 4a):
  130 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,2)=RB(1,2)-RB(1,3)
      RB(2,2)=RB(2,2)-RB(2,3)
      RB(3,2)=RB(3,2)-RB(3,3)
      ABSR=SQRT(RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2))
      IF (ABSR>(ABS2+TINY)) GOTO 132
      IF (ABSR<(ABS2-TINY)) THEN
         IFLAG=1
         ABS2=ABSR
      END IF
      ILOOP=1
      GOTO 130
  132 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,2)=RB(1,2)+RB(1,3)
      RB(2,2)=RB(2,2)+RB(2,3)
      RB(3,2)=RB(3,2)+RB(3,3)
! Nothing found ? Try A2+ICOUNT*A3 (loop 4b):
      IF (ILOOP==0) THEN
  133    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,2)=RB(1,2)+RB(1,3)
         RB(2,2)=RB(2,2)+RB(2,3)
         RB(3,2)=RB(3,2)+RB(3,3)
         ABSR=SQRT(RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2))
         IF (ABSR>(ABS2+TINY)) GOTO 135
         IF (ABSR<(ABS2-TINY)) THEN
            IFLAG=1
            ABS2=ABSR
         END IF
         GOTO 133
  135    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,2)=RB(1,2)-RB(1,3)
         RB(2,2)=RB(2,2)-RB(2,3)
         RB(3,2)=RB(3,2)-RB(3,3)
      END IF
      ILOOP=0
! Try A3-ICOUNT*A1 (loop 5a):
  140 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,3)=RB(1,3)-RB(1,1)
      RB(2,3)=RB(2,3)-RB(2,1)
      RB(3,3)=RB(3,3)-RB(3,1)
      ABSR=SQRT(RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3))
      IF (ABSR>(ABS3+TINY)) GOTO 142
      IF (ABSR<(ABS3-TINY)) THEN
         IFLAG=1
         ABS3=ABSR
      END IF
      ILOOP=1
      GOTO 140
  142 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,3)=RB(1,3)+RB(1,1)
      RB(2,3)=RB(2,3)+RB(2,1)
      RB(3,3)=RB(3,3)+RB(3,1)
! Nothing found ? Try A3+ICOUNT*A1 (loop 5b):
      IF (ILOOP==0) THEN
  143    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,3)=RB(1,3)+RB(1,1)
         RB(2,3)=RB(2,3)+RB(2,1)
         RB(3,3)=RB(3,3)+RB(3,1)
         ABSR=SQRT(RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3))
         IF (ABSR>(ABS3+TINY)) GOTO 145
         IF (ABSR<(ABS3-TINY)) THEN
            IFLAG=1
            ABS3=ABSR
         END IF
         GOTO 143
  145    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,3)=RB(1,3)-RB(1,1)
         RB(2,3)=RB(2,3)-RB(2,1)
         RB(3,3)=RB(3,3)-RB(3,1)
      END IF
      ILOOP=0
! Try A3-ICOUNT*A2 (loop 6a):
  150 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,3)=RB(1,3)-RB(1,2)
      RB(2,3)=RB(2,3)-RB(2,2)
      RB(3,3)=RB(3,3)-RB(3,2)
      ABSR=SQRT(RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3))
      IF (ABSR>(ABS3+TINY)) GOTO 152
      IF (ABSR<(ABS3-TINY)) THEN
         IFLAG=1
         ABS3=ABSR
      END IF
      ILOOP=1
      GOTO 150
  152 CONTINUE
      ICOUNT=ICOUNT+1
      RB(1,3)=RB(1,3)+RB(1,2)
      RB(2,3)=RB(2,3)+RB(2,2)
      RB(3,3)=RB(3,3)+RB(3,2)
! Nothing found ? Try A3+ICOUNT*A2 (loop 6b):
      IF (ILOOP==0) THEN
  153    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,3)=RB(1,3)+RB(1,2)
         RB(2,3)=RB(2,3)+RB(2,2)
         RB(3,3)=RB(3,3)+RB(3,2)
         ABSR=SQRT(RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3))
         IF (ABSR>(ABS3+TINY)) GOTO 155
         IF (ABSR<(ABS3-TINY)) THEN
            IFLAG=1
            ABS3=ABSR
         END IF
         GOTO 153
  155    CONTINUE
         ICOUNT=ICOUNT+1
         RB(1,3)=RB(1,3)-RB(1,2)
         RB(2,3)=RB(2,3)-RB(2,2)
         RB(3,3)=RB(3,3)-RB(3,2)
      END IF
! Nothing found during loops 1a-6b ? Hurray, we have got it !!
      IF (IFLAG==0) GOTO 200
! Next iteration ... :
  199 ABS1=SQRT(RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1))
      ABS2=SQRT(RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2))
      ABS3=SQRT(RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3))
      IFLAG=0
      ILOOP=0
      GOTO 100
  200 CONTINUE

! Adjustement of the basis to right-hand-sense (by inversion of all
! three lattice vectors if necessary ...):
      CALL CELVOL(RB(1,1),RB(1,2),RB(1,3),OMEGA)
      IF (OMEGA<0) THEN
         RB(1,1)=-1._q*RB(1,1)
         RB(2,1)=-1._q*RB(2,1)
         RB(3,1)=-1._q*RB(3,1)
         RB(1,2)=-1._q*RB(1,2)
         RB(2,2)=-1._q*RB(2,2)
         RB(3,2)=-1._q*RB(3,2)
         RB(1,3)=-1._q*RB(1,3)
         RB(2,3)=-1._q*RB(2,3)
         RB(3,3)=-1._q*RB(3,3)
      END IF

! Warning extremly strange inputs (very large/small ratios 'c/a', 'b/a'
! or cosines very close to +/-1) could sometimes lead to some problems
! due to the finite precision of the machine ...). Here we test some
! parameters and if we find 'strange parameters' we (DO NOT) stop!
      ABS1=RB(1,1)*RB(1,1)+RB(2,1)*RB(2,1)+RB(3,1)*RB(3,1)
      ABS2=RB(1,2)*RB(1,2)+RB(2,2)*RB(2,2)+RB(3,2)*RB(3,2)
      ABS3=RB(1,3)*RB(1,3)+RB(2,3)*RB(2,3)+RB(3,3)*RB(3,3)
      A2XA1=RB(1,1)*RB(1,2)+RB(2,1)*RB(2,2)+RB(3,1)*RB(3,2)
      A3XA1=RB(1,1)*RB(1,3)+RB(2,1)*RB(2,3)+RB(3,1)*RB(3,3)
      A3XA2=RB(1,2)*RB(1,3)+RB(2,2)*RB(2,3)+RB(3,2)*RB(3,3)
      ABS1=SQRT(ABS1)
      ABS2=SQRT(ABS2)
      ABS3=SQRT(ABS3)
      A1XA2=ABS(A2XA1/ABS1/ABS2)
      A1XA3=ABS(A3XA1/ABS1/ABS3)
      A2XA3=ABS(A3XA2/ABS2/ABS3)
      IF (A1XA2>0.995_q) ISTRGE=ISTRGE+1
      IF (A2XA3>0.995_q) ISTRGE=ISTRGE+1
      IF (A1XA3>0.995_q) ISTRGE=ISTRGE+1
      IF ((ABS3/ABS1)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS2/ABS1)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS3/ABS2)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS1/ABS3)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS1/ABS2)<0.001_q) ISTRGE=ISTRGE+1
      IF ((ABS2/ABS3)<0.001_q) ISTRGE=ISTRGE+1

! Now construct several possible combinations of lattice vectors
! and test them all for specific criterions ... :

            !write(*,*) "RB1=",RB(1,1),RB(2,1),RB(3,1) !mohan
            !write(*,*) "RB2=",RB(1,2),RB(2,2),RB(3,2) !mohan
            !write(*,*) "RB3=",RB(1,3),RB(2,3),RB(3,3) !mohan
      ITYP=15
      COS1=1._q
      COS2=1._q
      COS3=1._q
      DO N9=-2,2
      DO N8=-2,2
      DO N7=-2,2
       DO N6=-2,2
       DO N5=-2,2
       ID1=N5*N9-N6*N8
       DO N4=-2,2
        ID2=N6*N7-N4*N9
        ID3=N4*N8-N5*N7
        DO N3=-2,2
        ID4=N3*ID3
        DO N2=-2,2
        ID5=N2*ID2+ID4
        DO N1=-2,2
         IF ((N1*ID1+ID5)==1) THEN
            XB(1,1)=N1*RB(1,1)+N2*RB(1,2)+N3*RB(1,3)
            XB(2,1)=N1*RB(2,1)+N2*RB(2,2)+N3*RB(2,3)
            XB(3,1)=N1*RB(3,1)+N2*RB(3,2)+N3*RB(3,3)
            XB(1,2)=N4*RB(1,1)+N5*RB(1,2)+N6*RB(1,3)
            XB(2,2)=N4*RB(2,1)+N5*RB(2,2)+N6*RB(2,3)
            XB(3,2)=N4*RB(3,1)+N5*RB(3,2)+N6*RB(3,3)
            XB(1,3)=N7*RB(1,1)+N8*RB(1,2)+N9*RB(1,3)
            XB(2,3)=N7*RB(2,1)+N8*RB(2,2)+N9*RB(2,3)
            XB(3,3)=N7*RB(3,1)+N8*RB(3,2)+N9*RB(3,3)

! Lengths of the axes and cosines between the axes:
            ABS1=XB(1,1)*XB(1,1)+XB(2,1)*XB(2,1)+XB(3,1)*XB(3,1)
            ABS2=XB(1,2)*XB(1,2)+XB(2,2)*XB(2,2)+XB(3,2)*XB(3,2)
            ABS3=XB(1,3)*XB(1,3)+XB(2,3)*XB(2,3)+XB(3,3)*XB(3,3)
            A2XA1=XB(1,1)*XB(1,2)+XB(2,1)*XB(2,2)+XB(3,1)*XB(3,2)
            A3XA1=XB(1,1)*XB(1,3)+XB(2,1)*XB(2,3)+XB(3,1)*XB(3,3)
            A3XA2=XB(1,2)*XB(1,3)+XB(2,2)*XB(2,3)+XB(3,2)*XB(3,3)
            ABS2M1=SQRT(ABS1+ABS2-2._q*A2XA1)
            ABS3M1=SQRT(ABS1+ABS3-2._q*A3XA1)
            ABS3M2=SQRT(ABS2+ABS3-2._q*A3XA2)
            ABS2P1=SQRT(ABS1+ABS2+2._q*A2XA1)
            ABS3P1=SQRT(ABS1+ABS3+2._q*A3XA1)
            ABS3P2=SQRT(ABS2+ABS3+2._q*A3XA2)
            A1P2M3=SQRT(ABS1+ABS2+ABS3+2._q*A2XA1-2._q*A3XA1-2._q*A3XA2)
            A2P3M1=SQRT(ABS1+ABS2+ABS3+2._q*A3XA2-2._q*A2XA1-2._q*A3XA1)
            A3P1M2=SQRT(ABS1+ABS2+ABS3+2._q*A3XA1-2._q*A3XA2-2._q*A2XA1)
            ABS123=A2XA1+A3XA1+A3XA2
            ABS1=SQRT(ABS1)
            ABS2=SQRT(ABS2)
            ABS3=SQRT(ABS3)
            A2XA1=A2XA1/ABS1/ABS2
            A3XA1=A3XA1/ABS1/ABS3
            A3XA2=A3XA2/ABS2/ABS3
            A1XA2=ABS(A2XA1)
            A1XA3=ABS(A3XA1)
            A2XA3=ABS(A3XA2)
! Warning extremly strange inputs (very large/small ratios 'c/a', 'b/a'
! or cosines very close to +/-1) could sometimes lead to some problems
! due to the finite precision of the machine ...). Here we test some
! parameters and if we find 'strange parameters', we (DO NOT) stop!

! Following the identification of the crystal type, transformation of
! the basis to 'crystallographic standard settings' and estimate of the
! cell dimensions (using the conventions of routine LATGEN):
            CELLDM(1)=0._q
            CELLDM(2)=0._q
            CELLDM(3)=0._q
            CELLDM(4)=0._q
            CELLDM(5)=0._q
            CELLDM(6)=0._q
            IBRAV=15
! Crystal classes with A1*A2=A1*A3=A2*A3:
            IF ((ABS(A2XA1-A3XA1)<TINY).AND. &
     &                                  (ABS(A2XA1-A3XA2)<TINY)) THEN
! Crystal classes with A1*A2=A1*A3=A2*A3 and |A1|=|A2|=|A3|:
               IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &                                    (ABS(ABS2-ABS3)<TINY)) THEN
! Orthogonal axes means simple cubic: (IBRAV=1)
                  IF (A1XA2<TINY) THEN
                     IBRAV=1
                     CELLDM(1)=ABS1
! Cosine(alpha) equal to -1/3 means body-centered cubic: (IBRAV=2)
                  ELSE IF (ABS(A2XA1+(1._q/3._q))<TINY) THEN
                     IBRAV=2
                     CELLDM(1)=ABS1*2._q/SQRT(3._q)
! Cosine(alpha) equal to +1/2 means face-centered cubic: (IBRAV=3)
                  ELSE IF (ABS(A2XA1-0.5_q)<TINY) THEN
                     IBRAV=3
                     CELLDM(1)=ABS1*SQRT(2._q)
! All other angles means trigonal (rhomboedric): (IBRAV=7)
                  ELSE
                     IBRAV=7
                     CELLDM(1)=ABS1
                     CELLDM(4)=A2XA1
                  END IF
! Crystal classes with inequal length of lattice vectors but also with
! A1*A2=A1*A3=A2*A3:
! Orthogonal axes:
               ELSE IF (A1XA2<TINY) THEN
! Two axes with equal lengths means simple tetragonal: (IBRAV=5)
! Adjustment: 'c-axis' shall be the special axis.
                  IF (ABS(ABS1-ABS2)<TINY) THEN
                     IBRAV=5
                     CELLDM(1)=ABS1
                     CELLDM(3)=ABS3/ABS1
! No axes with equal lengths means simple orthorhombic (IBRAV=8):
! Adjustment: Sort the axis by increasing lengths:
                  ELSE IF (((ABS3-ABS2)>TINY).AND. &
     &                                   ((ABS2-ABS1)>TINY)) THEN
                     IBRAV=8
                     CELLDM(1)=ABS1
                     CELLDM(2)=ABS2/ABS1
                     CELLDM(3)=ABS3/ABS1
                  END IF
               END IF
! Crystal classes with A1*A3=A2*A3=/A1*A2:
            ELSE IF (ABS(A3XA1-A3XA2)<TINY) THEN
! One axis orthogonal with respect to the other two axes:
               IF (A1XA3<TINY) THEN
! Equal length of the two nonorthogonal axes:
                  IF (ABS(ABS1-ABS2)<TINY) THEN
! Cosine(alpha) equal to -1/2 means hexagonal: (IBRAV=4)
! Adjustment: 'c-axis' shall be the special axis.
                     IF (ABS(A2XA1+0.5_q)<TINY) THEN
                        IBRAV=4
                        CELLDM(1)=ABS1
                        CELLDM(3)=ABS3/ABS1
! Other angles mean base-centered orthorhombic: (IBRAV=11)
! Adjustment: Cosine between A1 and A2 shall be lower than zero, the
!             'c-axis' shall be the special axis.
                     ELSE IF (A2XA1<(-1._q*TINY)) THEN
                        IBRAV=11
                        CELLDM(1)=ABS2P1
                        CELLDM(2)=ABS2M1/ABS2P1
                        CELLDM(3)=ABS3/ABS2P1
                     END IF
! Different length of the two axes means simple monoclinic (IBRAV=12):
! Adjustment: Cosine(gamma) should be lower than zero, special axis
!             shall be the 'b-axis'(!!) and |A1|<|A3|:
                  ELSE IF ((A2XA1<(-1._q*TINY)).AND. &
     &                                       ((ABS1-ABS2)>TINY)) THEN
                     IBRAV=12
                     CELLDM(1)=ABS2
                     CELLDM(2)=ABS3/ABS2
                     CELLDM(3)=ABS1/ABS2
                     CELLDM(5)=A2XA1
                     YB(1)=XB(1,3)
                     YB(2)=XB(2,3)
                     YB(3)=XB(3,3)
                     XB(1,3)=XB(1,1)
                     XB(2,3)=XB(2,1)
                     XB(3,3)=XB(3,1)
                     XB(1,1)=XB(1,2)
                     XB(2,1)=XB(2,2)
                     XB(3,1)=XB(3,2)
                     XB(1,2)=YB(1)
                     XB(2,2)=YB(2)
                     XB(3,2)=YB(3)
                  END IF
! Arbitrary angles between the axes:
               ELSE
! |A1|=|A2|=|A3| means body-centered tetragonal (IBRAV=6):
! Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
! orthogonal to one another and (adjustment!): |A1+A3|=|A2+A3|/=|A1+A2|
                  IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &                   (ABS(ABS2-ABS3)<TINY).AND. &
     &                      (ABS(ABS3P1-ABS3P2)<TINY).AND. &
     &                          (ABS(ABS2P1-ABS3P1)>TINY).AND. &
     &                             (ABS(ABS3*ABS3+ABS123)<TINY)) THEN
                     IBRAV=6
                     CELLDM(1)=ABS3P1
                     CELLDM(3)=ABS2P1/ABS3P1
! |A1|=|A2|=/|A3| means base-centered monoclinic (IBRAV=13):
! Adjustement: The cosine between A1 and A3 as well as the cosine
!              between A2 and A3 should be lower than zero.
                  ELSE IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &                                (A3XA1<(-1._q*TINY)).AND. &
     &                                       (A3XA2<(-1._q*TINY))) THEN
                     IBRAV=13
                     CELLDM(1)=ABS2P1
                     CELLDM(2)=ABS2M1/ABS2P1
                     CELLDM(3)=ABS3/ABS2P1
                     CELLDM(5)=A3XA2
                  END IF
               END IF
! Crystal classes with A1*A2=/A1*A3=/A2*A3
            ELSE
! |A1|=|A2|=|A3| means body-centered orthorhombic (IBRAV=9):
! Further additional criterions are: (A1+A2), (A1+A3) and (A2+A3) are
! orthogonal to one another and (adjustment!): |A1+A2|>|A1+A3|>|A2+A3|
               IF ((ABS(ABS1-ABS2)<TINY).AND. &
     &                (ABS(ABS2-ABS3)<TINY).AND. &
     &                    ((ABS3P1-ABS3P2)>TINY).AND. &
     &                         ((ABS2P1-ABS3P1)>TINY).AND. &
     &                             (ABS(ABS3*ABS3+ABS123)<TINY)) THEN
                  IBRAV=9
                  CELLDM(1)=ABS3P2
                  CELLDM(2)=ABS3P1/ABS3P2
                  CELLDM(3)=ABS2P1/ABS3P2
! |A1|=|A2-A3| and |A2|=|A1-A3| and |A3|=|A1-A2| means face-centered
! orthorhombic (IBRAV=10):
! Adjustment: |A1+A2-A3|>|A1+A3-A2|>|A2+A3-A1|
               ELSE IF ((ABS(ABS2M1-ABS3)<TINY).AND. &
     &                      (ABS(ABS3M1-ABS2)<TINY).AND. &
     &                          (ABS(ABS3M2-ABS1)<TINY).AND. &
     &                               ((A1P2M3-A3P1M2)>TINY).AND. &
     &                                   ((A3P1M2-A2P3M1)>TINY)) THEN
                  IBRAV=10
                  CELLDM(1)=A2P3M1
                  CELLDM(2)=A3P1M2/A2P3M1
                  CELLDM(3)=A1P2M3/A2P3M1
! Now there exists only one further possibility - triclinic (IBRAV=14):
! Adjustment: All three cosines shall be greater than zero and ordered:
               ELSE IF ((A2XA1>A3XA1).AND.(A3XA1>A3XA2).AND. &
     &                                             (A3XA2>TINY)) THEN
                  IBRAV=14
                  CELLDM(1)=ABS1
                  CELLDM(2)=ABS2/ABS1
                  CELLDM(3)=ABS3/ABS1
                  CELLDM(4)=A3XA2
                  CELLDM(5)=A3XA1
                  CELLDM(6)=A2XA1
               END IF
            END IF

! Have we found some lattice with lower index IBRAV as all previous?
            IF ((IBRAV<ITYP).OR.((IBRAV==ITYP).AND. &
     &                 ((ABS(CELLDM(4))<(COS1-1.E-9_q)).AND. &
     &                      (ABS(CELLDM(5))<(COS2-1.E-9_q)).AND. &
     &                          (ABS(CELLDM(6))<(COS3-1.E-9_q))))) THEN
            write(*,*) "test XB again" !mohan
            write(*,*) "IBRAV=",IBRAV
            write(*,*) "ITYP=",ITYP
            write(*,*) "ABS(CELLDM(4))=",ABS(CELLDM(4))
            write(*,*) "ABS(CELLDM(5))=",ABS(CELLDM(5))
            write(*,*) "ABS(CELLDM(6))=",ABS(CELLDM(6))
            write(*,*) "(COS1-1.E-9_q)=",(COS1-1.E-9_q)
            write(*,*) "(COS2-1.E-9_q)=",(COS2-1.E-9_q)
            write(*,*) "(COS3-1.E-9_q)=",(COS3-1.E-9_q)

               ITYP=IBRAV
               COS1=ABS(CELLDM(4))
               COS2=ABS(CELLDM(5))
               COS3=ABS(CELLDM(6))
               CELDIM(1)=CELLDM(1)
               CELDIM(2)=CELLDM(2)
               CELDIM(3)=CELLDM(3)
               CELDIM(4)=CELLDM(4)
               CELDIM(5)=CELLDM(5)
               CELDIM(6)=CELLDM(6)
               A1(1)=XB(1,1)
               A1(2)=XB(2,1)
               A1(3)=XB(3,1)
               A2(1)=XB(1,2)
               A2(2)=XB(2,2)
               A2(3)=XB(3,2)
               A3(1)=XB(1,3)
               A3(2)=XB(2,3)
               A3(3)=XB(3,3)

            write(*,*) N1,N2,N3
            write(*,*) N4,N5,N6
            write(*,*) N7,N8,N9 
            write(*,*) "XB1=",XB(1,1),XB(2,1),XB(3,1) !mohan
            write(*,*) "XB2=",XB(1,2),XB(2,2),XB(3,2) !mohan
            write(*,*) "XB3=",XB(1,3),XB(2,3),XB(3,3) !mohan

            END IF
         END IF
        ENDDO
        ENDDO
        ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDDO
      ENDDO
      ENDDO

      write(*,*) "LATTYP 2" !mohan
      write(*,*) "A1=",A1(1),A1(2),A1(3)
      write(*,*) "A2=",A2(1),A2(2),A2(3)
      write(*,*) "A3=",A3(1),A3(2),A3(3)

! Now we got the 'final answer' which lattice type we have. Now we have
! of course changed the original basis vector totally (interchanged or
! inverted and probably linearly combined ...?). If the original input
! vectors would have given the same answer (same IBRAV, same CELDIM) it
! would probably be nice to take them instead of a total different basis
! and so check against the result from our first trial at the beginning.
      IF (IBRINP==ITYP) THEN
         IFLAG=0
! Consistency check ... :
         IF (ABS(CELDIM(1)-CDMINP(1))>TINY) IFLAG=1
         IF (ABS(CELDIM(2)-CDMINP(2))>TINY) IFLAG=1
         IF (ABS(CELDIM(3)-CDMINP(3))>TINY) IFLAG=1
         IF (ABS(CELDIM(4)-CDMINP(4))>TINY) IFLAG=1
         IF (ABS(CELDIM(5)-CDMINP(5))>TINY) IFLAG=1
         IF (ABS(CELDIM(6)-CDMINP(6))>TINY) IFLAG=1
         IF (IFLAG==0) THEN
! Now take really the original input vectors ... :
            A1(1)=SA1(1)
            A1(2)=SA1(2)
            A1(3)=SA1(3)
            A2(1)=SA2(1)
            A2(2)=SA2(2)
            A2(3)=SA2(3)
            A3(1)=SA3(1)
            A3(2)=SA3(2)
            A3(3)=SA3(3)
         ELSE IF ((PRCHAN>=0).AND.(PRCHAN<=99)) THEN
            WRITE(PRCHAN,*) ' Warning from LATTYP: Got some problem '// &
     &                                           'with cell dimensions!'
            WRITE(PRCHAN,*) ' Tried to take original basis but could '// &
     &                                          'not verify dimensions!'
         END IF
      END IF

! Print result (or not ...)
      IF ((PRCHAN<0).OR.(PRCHAN>99)) RETURN
      END IF
      WRITE(PRCHAN,*) ' '
      WRITE(PRCHAN,*) ' Lattice vectors:'
      WRITE(PRCHAN,*) ' '
      WRITE(PRCHAN,400) ' A1 = (',A1(1),',',A1(2),',',A1(3),')'
      WRITE(PRCHAN,400) ' A2 = (',A2(1),',',A2(2),',',A2(3),')'
      WRITE(PRCHAN,400) ' A3 = (',A3(1),',',A3(2),',',A3(3),')'
  400 FORMAT(A,F15.10,A,F15.10,A,F15.10,A)

      RETURN
      END
