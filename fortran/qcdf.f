c
C     **** MODIFIED FOR ANTI-PBC  ****
C     **** ADJUSTED 7/11/87 TO USE FORTRAN 2.2.0 ****
C
C     **** 12/14/87 CHANGED QCD2A TO BREAK HAM   ****
C     **** INTO FREE AND INTERACTING PIECES TO   ****
C     **** ALLOW STORAGE AND RECALCULATION FOR   ****
C     **** DIFFERING LAMBDA; FREE PART IS IN     ****
C     **** HAM0, HNU0                            ****
C
C     **** 6/24/88 MODIFIED QCD2A2 SO THAT COLOR ****
C     **** SUMS ARE PERFORMED DIAGRAMMATICALLY   ****
C     **** RATHER THAN ITERATIVELY               ****
C
C     **** 6/25/88 MODIFIED QCDI SO THAT HAVE    ****
C     **** OPTION OF USING PAULI VILLARS TYPE    ****
C     **** CUTOFF ON STATES                      ****
c
c     **** 5/8/90 modified qcdf to include       ****
c     **** arbitrary number of quark flavors     ****
c     **** ???? denote places where modified qcdi to include
c     **** flavor, for debugging purposes; should be removed
c     **** after convinced program is working.
c
c     (C) Kent Hornbostel 1993.  All Rights Reserved.
c     The routines ESRTR8 TQR8 TRR8 are property of the authors
c     of Numerical Recipes.
C
C     MAIN PROG
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      common/masses/rmq(03), iflv(03)
      COMMON/CTOF/CUTOFF
      COMMON/NST/NUMSTA
      COMMON/HM/HAM0(03,6902,6902),HAM(6902,6902)
c     COMMON/HN/HNU0(03,6902,6902),HNU(6902,6902)
      COMMON/HN/HNU0(03,6902),HNU(6902,6902)
      COMMON/NM/HNORM(6902,6902)
      COMMON/NOBAR/IBAB
      COMMON/NOZERO/IB0,ID0
      COMMON/PARLIM/LPN
      COMMON/PRTOPS/INTPRT
      INTEGER N,NF,B,K
      integer bfr(3,2)
      real*8 rbfr(3)
      REAL*8 W(6902),Z(6902,6902)
C
c         big number:  9523 (see comments just below)
c         array size for various intermed arrays
c
c         another big number:  12552
c         array size for various intermed arrays
c         (check whether these should be tied together)
c
c         another big number:  25001
c         array size for various intermed arrays
c         (check whether these should be tied together)
c         
c
      nmstmx=6902
c     **** max number of states; fixes dimensions of ****
c     **** matrices ham, hnu, hnorm, w, z etc;       ****
c     **** to facilitate changes in this dimension,  ****
c     **** check that neither this number (nmstmx),  ****
c     **** nor the number it's changed to, appear as ****
c     **** anything other than the dimension of      ****
c     **** the appropriate matrices; ie, if changing ****
c     **** this max dim. from n02 to m02, be sure    ****
c     **** that before switching, n02 doesn't appear ****
c     **** anywhere (such as a statement number)     ****
c     **** (where 'n' and 'm' here are some integers)****
      nmxfr = 27608
c     **** four times nmstmx
      nmflmx=03
c     **** max num of flavors   ****
c
      IVCTMX = 75
c     **** max num of evects to print   ****
C
      NSONLY=1
C     **** SET NSONLY EQ 0 IF ONLY WANT NORMALIZED ****
C     **** STATES WITHOUT SOLVING HAMILTONIAN      ****
C     **** THIS WILL ALSO PRINT OUT INFO NEEDED    ****
C     **** FOR WFBIG FORTRAN                       ****
C
C
C     **** OPEN FILE IN WHICH EITHER HAM0 AND HAM   ****
C     **** ARE READ IF THE FILE EXISTS, OR STORED   ****
C     **** IF IT DOESN'T                            ****
c     OPEN(UNIT=14,ACCESS='DIRECT',FILE='qcdf.ham',RECL=70,
c    >     FORM='FORMATTED')
      jstat=0
      OPEN(unit=14,file='qcdf.ham',status='unknown')
c     OPEN(UNIT=7,FILE='/dev/console')
c     OPEN(UNIT=7,FILE='/dev/ttya')
c     **** have changed unit 7 to '*'   ****
      OPEN(UNIT=15,file='qcdf.out',status='unknown')
C
C
      WRITE(*,400)
  400 FORMAT(' DOES A HAMILTONIAN FILE EXIST? YES: 1; NO: 0')
      READ(*,*) IHFILE
C
      IF(IHFILE.NE.1) THEN
C     **** COMPUTE STATES AND HAMILTONIAN ****
C
C
C     **** INITIALIZE PARAMS                           ****
C     **** N=NUM COLORS, B=BARYON NUM, 2K=TOTAL LC MOM ****
C
C     **** INPUT N,NF,K,B,RLAMB FROM TERMINAL ****
      WRITE(*,5)
    5 FORMAT(' INPUT N:' )
      READ(*,*)N
      WRITE(*,7)
    7 FORMAT(' INPUT NF:' )
      READ(*,*)NF
      if(nf.gt.nmflmx) then
       WRITE(15,8)
    8  FORMAT(' nf exceeds nmflmx in main prog ')
       stop
      else
      endif
      WRITE(*,6)
    6 FORMAT(' INPUT B:' )
      READ(*,*)B
      WRITE(*,9)N,B
    9 FORMAT(' N,B: ',2I5)
      WRITE(*,70)
   70 FORMAT(' INPUT LAMBDA (0.0<LAMBDA<1.0):')
      READ(*,*)RLAMB
c
c  ????????????????????????????????????????????????
      rmq(1)=1.d0
      if(nf.gt.1) then
c     **** input quark masses ****
c     **** first quark mass is 1*(1-lambda**2)**.5 ; ****
c     **** the rest are rmq(i=2 to nf) * the first    ****
       do 85 nqk=2,nf
        WRITE(*,75)nqk
   75   FORMAT(' input ratio of quark mass',i3,' to quark mass  1:')
        read(*,*)rmq(nqk)
   85  continue
      else
      endif
C
      if(nf.gt.1) then
c     **** input overall flavor qm numbers  ****
       do 86 nfl=1,nf
        WRITE(*,76) nfl
   76   FORMAT(' input flavor qm number for flavor ', i3)
        read(*,*)iflv(nfl)
   86  continue
      else
      endif
C
C     WRITE(*,80)
C  80 FORMAT(' ARE THE QUARKS MASSIVE? YES: 1; NO: 0')
C     READ(*,*)MASS
      MASS=1
C
C
C     WRITE(*,60)
C  60 FORMAT(' ALLOW ZERO MOM D"S? YES: 1; NO: 0')
C     READ(*,*)ID0
      ID0=0
C
C
C     WRITE(*,62)
C  62 FORMAT(' ALLOW ZERO MOM B"S? YES: 1; NO: 0')
C     READ(*,*)IB0
      IB0=0
C
      WRITE(*,95)
   95 FORMAT(' CUTOFF/MASS SQD:  (LE.0 =  NO CUTOFF')
      READ(*,*)CUTOFF
      WRITE(*,18)CUTOFF
   18 FORMAT(' CUTOFF:',e24.16)
C
      WRITE(*,90)
   90 FORMAT(' LIMIT PARTICLE NUMBER TO:  (0 = NO LIMIT')
      READ(*,*)LPN
      WRITE(*,13)LPN
   13 FORMAT(' LPN:',I5)
C
C     **** CALC MIN K ALLOWED ****
c ?????????????????????????????????
c     IF(IB0.EQ.1) THEN
c      KMIN=N*B*(B-1)/2
c     ELSE
       IF(B.NE.0) THEN
C       KMIN=N*B*(B+1)/2
c       KMIN=N*B*B
        nqmin=n*b
        kmin=minmom(nqmin)
       ELSE
        KMIN=2
       ENDIF
c     ENDIF
C
      WRITE(*,10)KMIN
   10 FORMAT(' INPUT K  (GE ',I4,'):')
      READ(*,*)K
      WRITE(*,14)K
   14 FORMAT(' K:',I5)
C
      WRITE(15,20)N,NF,B,K,RLAMB
   20 FORMAT(' N,NF,B,K,LAMBDA:  ',4I5,e24.16)
C
C     WRITE(15,22)IB0,ID0,LPN
C  22 FORMAT(' IB0,ID0,LPN:  ',3I5)
      WRITE(15,22)LPN
   22 FORMAT(' LPN:  ',I5)
      WRITE(15,24)CUTOFF
   24 FORMAT(' CUTOFF/MASS SQD: ',e24.16)
c
      WRITE(15,1200)
 1200 FORMAT(' ratio of quark masses to first mass (rmq(nf)): ')
      WRITE(15,1210)(rmq(i),i=1,nf)
 1210 FORMAT(6e24.16)
C
      WRITE(15,1201)
 1201 FORMAT(' flavor numbers : ')
      WRITE(15,1211)(iflv(i),i=1,nf)
 1211 FORMAT(6i3)
C
      IF (KMIN.GT.K) THEN
       WRITE(15,25)
   25  FORMAT(' K IS TOO SMALL FOR N AND B INPUT ')
       STOP
      ELSE
      ENDIF
C
C     WRITE(*,30)
C  30 FORMAT(' ALLOW EXTRA BAR,ABAR PAIRS? YES: 1 ;NO: 0')
C     READ(*,*)IBAB
      IBAB=0
C
      IPTHAM=0
C     **** SET THIS .NE.0 IF WANT HAMILT. PRINTED ****
C
C     **** SET INTPRT.NE.0 IF WANT PRINTOUT OF     ****
C     **** OF INTERMEDIATE RESULTS FOR H0,H1, ETC. ****
C     WRITE(*,37)
C  37 FORMAT(' PRINT INTERMED HAM VALUES? YES: 1; NO: 0')
C     READ(*,*)INTPRT
      INTPRT=0
cf    INTPRT=1
C
C
C     ??????????????????????????????????????????????
C
      iflv(1) = N*B
      CALL QCDSTA
      WRITE(*,1000)
 1000 FORMAT('  CALLED SUB QCDSTA ')
C     **** GENERATE STATES ****
C
C     **** CHECK IF THERE ARE ANY STATES GENERATED ****
C     **** BEFORE WEEDING                          ****
      WRITE(15,200) NUMSTA
      WRITE(*,200) NUMSTA
  200 FORMAT(' NUMBER OF STATES BEFORE WEEDING: ',I6)
      IF (NUMSTA.EQ.0) THEN
       STOP
      ELSE
      ENDIF
      IF (NUMSTA.gt.nmstmx) THEN
       WRITE(15,207)
  207  FORMAT(' numsta exceeds nmstmx in sub main ')
       STOP
      ELSE
      ENDIF
C
c  ????????????????????????????????????
cf    CALL PRNTST
C     **** PRINT STATES ****
c     temp=1.d0
c     if(temp.ne.0) go to 2000
C
      CALL SLFN
      WRITE(*,1010)
 1010 FORMAT('  CALLED SUB SLFN   ')
C     **** EVAL SELFEN USED IN HAMILT. ****
C
      CALL CLRDIS(0)
      WRITE(*,1015)
 1015 FORMAT('  CALLED SUB CLRDIS(0) ')
C     **** DISTRIB COLOR AND EVAL NORM MAT ****
C
C
C$    WRITE(15,100)
C$100 FORMAT(/,' NORM MATRIX: ',/)
C$    CALL PRTMAT(NUMSTA,HNORM,0)
C     **** PRINT NORM MATRIX ****
C
C$    WRITE(15,120)
C$120 FORMAT(//,' WEEDR:')
      CALL WEEDR
      WRITE(*,1030)
 1030 FORMAT('  CALLED SUB WEEDR     ')
C     **** WEED REDUNDANT STATES ****
C
C     CALL PRNTST
C
C$    WRITE(15,101)
C$101 FORMAT(/,' NORM MATRIX: ',/)
C$    CALL PRTMAT(NUMSTA,HNORM,1)
C
C
C     **** RUN THESE AGAIN IF WANT TO ENSURE     ****
C     **** THAT WEEDR PROPERLY STORED NEW STATES ****
C
C     CALL CLRDIS(0)
C
C     WRITE(15,112)
C 112 FORMAT(/,' NORM MATRIX: ',/)
C     CALL PRTMAT(NUMSTA,HNORM,0)
C     ****
C
C
C
C     **** AFTER DIAG HNORM, WEED STATES BY ****
C     **** FINDING ZERO E-VALS              ****
C$    WRITE(15,130)
C$130 FORMAT(//,' WEEDR2:')
      WRITE(*,1040)
 1040 FORMAT('  LOOK FOR ZERO EVALS IN NORM MATRIX  ')
      IWMX=2000
      IDR=1
C
      DO 40 IW=1,IWMX
       IF(IDR.NE.0) THEN
        CALL DIAG(HNORM,W,Z,IERR)
C       **** DIAGONALIZE NORM ****
        CALL PRTEIG(W,Z,IERR,0)
C       **** PRINT EIG VALS FROM DIAG ****
        WRITE(*,1050)
 1050   FORMAT('    CALLED SUB DIAG, PRTEIG    ')
C
        CALL WEEDR2(W,Z,IDR)
        WRITE(*,1060)
 1060   FORMAT('    CALLED SUB WEEDR2          ')
       ELSE
       ENDIF
C
   40 CONTINUE
C
   45 IF(IDR.NE.0) THEN
       WRITE(15,47)IWMX
   47  FORMAT(' WEEDR2 HAS RUN ', I5,' TIMES W-OUT FINISH.')
      ELSE
      ENDIF
C
C     **** CHECK IF THERE ARE ANY STATES GENERATED ****
C     **** AFTER WEEDING                           ****
      WRITE(15,250) NUMSTA
  250 FORMAT(' NUMBER OF STATES AFTER WEEDING: ',I6)
      IF (NUMSTA.EQ.0) THEN
       STOP
      ELSE
      ENDIF
C
C
      WRITE(15,103)
  103 FORMAT(/,' NORM MATRIX: ',/)
      CALL PRTSPM(NUMSTA,HNORM,0)
C
      CALL PRNTST
      WRITE(*,1080)
 1080 FORMAT('    CALLED SUB PRNTST    ')
C
C
      CALL NUZ(W,Z)
C     **** NORMALIZE Z BY DIVID BY EIGVALS OF NORM MAT ****
      WRITE(*,1085)
 1085 FORMAT('    CALLED SUB NUZ     ')
C
      WRITE(*,1135)
 1135 FORMAT('  WRT Z FOR NEW BASIS ')
      WRITE(15,105)
  105 FORMAT(/,' MATRIX Z WHICH GIVES NEW BASIS: ')
      WRITE(15,107)
  107 FORMAT(' (ORIG. STATE NUM. , ORTHON. STATE NUM.)',/)
      CALL PRTSPM(NUMSTA,Z,1)
C     **** PRINT MATRIX WHICH NORMALIZES STATES ****
C
       CALL WFST
C      **** PRINT INFO ON STATES FOR WF FORTRAN ****
       CALL PRZ(Z)
C      **** PRINT INFO ON Z FOR WF FORTRAN ****
C
C
      IF(NSONLY.EQ.0) THEN
C      **** ONLY WANT NORMALIZED STATES ****
C      CALL WFST
C      **** PRINT INFO ON STATES FOR WF FORTRAN ****
C      CALL PRZ(Z)
C      **** PRINT INFO ON Z FOR WF FORTRAN ****
       STOP
      ELSE
      ENDIF
C
c ?????????????????????????????????????
c     temp=1.d0
c     if(temp.ne.0) go to 2000
c
c
c
C     **** NOW EVALUATE HAMILTONIAN IN ****
C     **** PROPER BASIS                ****
      WRITE(*,1070)
 1070 FORMAT('  EVALUATE HAMILTONIAN    ')
C
C
      CALL CLRDIS(1)
C     **** CALC HAM ****
      WRITE(*,1090)
 1090 FORMAT('    CALLED SUB CLRDIS(1) ')
C
C     CALL NUHAM(W,Z)
      CALL NUHAM(Z)
C     **** CALC HNU (NEW HAMILT. IN ORTHON BASIS) ****
      WRITE(*,1100)
 1100 FORMAT('    CALLED SUB NUHAM     ')
C
C
C
C     WRITE(*,1135)
C1135 FORMAT('  WRT Z FOR NEW BASIS ')
C     WRITE(15,105)
C 105 FORMAT(/,' MATRIX Z WHICH GIVES NEW BASIS: ')
C     WRITE(15,107)
C 107 FORMAT(' (ORIG. STATE NUM. , ORTHON. STATE NUM.)',/)
C     CALL PRTSPM(NUMSTA,Z,1)
C     **** PRINT MATRIX WHICH NORMALIZES STATES ****
C
C
C     **** STORE HNU0, HNU (IN QCDF OUTPUT) ****
C     WRITE(*,1110)
C1110 FORMAT('  STORE HNU0,HNU IN QCDF OUTPUT   ')
C
C     WRITE(15,675)
C 675 FORMAT('     N     B    LPN    K   NUMSTA ')
C     WRITE(15,685)N,B,LPN,K,NUMSTA
C 685 FORMAT(5I6)
C     WRITE(15,665)
C 665 FORMAT('      HNU0            HNU       ')
C
C     DO 695 IR1=1,NUMSTA
C     DO 695 IR2=1,IR1
C      WRITE(15,705) HNU0(IR1,IR2),HNU(IR1,IR2)
C 705  FORMAT(D16.10,5X,D16.10)
C 695 CONTINUE
C
C
C     **** STORE HNU0, HNU (IN QCDF HAM ) ****
      WRITE(*,1120)
 1120 FORMAT('  STORE HNU0,HNU IN QCDF HAM      ')
C
c ????????????????????????????????????
c     WRITE(14,670,REC=1)
      WRITE(14,670)
  670 FORMAT('     N     NF    B    LPN    K   NUMSTA   CUTOFF')
c     WRITE(14,680,REC=2)N,B,LPN,K,NUMSTA,CUTOFF
      WRITE(14,680)N,NF,B,LPN,K,NUMSTA,CUTOFF
  680 FORMAT(6I6,6X,e24.16)
c
c     WRITE(14,660,REC=3)
c     IREC=4
C     **** RECORD LOCATION OF FIRST ELEMENT OF HAMILT. ****
c
c     **** store HNU0; only diag elements are non-zero ****
c     **** since HNU0 gives ham0 in orthonormal basis  ****
      WRITE(14,660)
c 660 FORMAT('      HNU0            HNU       ')
c 660 FORMAT('      HNU0(nf);    HNU       ')
  660 FORMAT('      HNU0(nf):   diagonal terms for each flavor   ')
      DO 690 IR1=1,NUMSTA
c     DO 690 IR2=1,IR1
c      WRITE(14,700,REC=IREC) HNU0(IR1,IR2),HNU(IR1,IR2)
c      WRITE(14,700) (HNU0(ifl,IR1,IR2),ifl=1,nf),HNU(IR1,IR2)
       WRITE(14,700) (HNU0(ifl,IR1),ifl=1,nf)
  700  FORMAT(6(e24.16,1X))
c     IREC=IREC+1
  690 CONTINUE
c
c     **** store HNU ****
      WRITE(14,667)
  667 FORMAT('      HNU:   non-zero elements of interaction hamilt. ')
      eps=1.d-15
      nbuff=0
      DO 693 IR1=1,NUMSTA
      DO 693 IR2=1,IR1
c      **** find next three nonzero elements of hnu ;  ****
c      **** put location in buffer; write next three   ****
       if(abs(hnu(ir1,ir2)).gt.eps) then
        nbuff=nbuff+1
        bfr(nbuff,1)=ir1
        bfr(nbuff,2)=ir2
        if(nbuff.eq.3) then
         WRITE(14,704)(bfr(i,1), bfr(i,2),
     >             hnu(bfr(i,1),bfr(i,2)), i=1,3)
         nbuff=0
        else
        endif
       else
       endif
  693 CONTINUE
      if(nbuff.eq.0) then
       WRITE(14,704)-1,-1,0.d0, -1,-1,0.d0, -1,-1,0.d0
c      **** -1,-1,0.d0 indicates end of list ****
      else
       if(nbuff.eq.1) then
        WRITE(14,704)bfr(1,1),bfr(1,2),hnu(bfr(1,1),bfr(1,2)),
     >   -1,-1,0.d0, -1,-1,0.d0
       else
        if(nbuff.eq.2) then
         WRITE(14,704)bfr(1,1),bfr(1,2),hnu(bfr(1,1),bfr(1,2)),
     >   bfr(2,1),bfr(2,2),hnu(bfr(2,1),bfr(2,2)), -1,-1,0.d0
        else
        endif
       endif
      endif
c      WRITE(14,703) (HNU(IR1,IR2),ir2=1,ir1)
c 703  FORMAT(4(D16.10,1X))
      if(numsta.ge.10000) then
       WRITE(*,707)
  707  FORMAT(' allowed numsta in stmnt 704 of main prog. exceeded') 
      else
      endif
  704 FORMAT(3(i6,i6,1x,e24.16,1X))
c
      CLOSE(UNIT=14)
C
C
c   ??????????????????????????????????????????
C     **** COMBINE THE FREE AND INTERACTING HAMILTONIANS ****
C     **** FOR UNNORMALIZED STATES                       ****
      RLMSQ=RLAMB*RLAMB
c     **** scale interacting ham by coupling const ****
      DO 300 I1=1,NUMSTA
      DO 300 I2=1,I1
       HAM(I1,I2)= RLMSQ*HAM(I1,I2)
      do 305 i3=1,nf
       HAM(I1,I2)= HAM(I1,I2)+
     >        (1.0D0 - RLMSQ)*rmq(i3)*rmq(i3)*HAM0(i3,I1,I2) 
  305 continue
c      HAM(I1,I2)= (1.0D0 - RLMSQ)*HAM0(I1,I2) + RLMSQ*HAM(I1,I2)
       HAM(I2,I1)= HAM(I1,I2)
  300 CONTINUE
c
c
      IF(IPTHAM.NE.0) THEN
       WRITE(*,1130)
 1130  FORMAT('  WRT HAM IN ORIG. BASIS ')
       WRITE(15,108)
  108  FORMAT(/,' HAMILTONIAN IN ORIGINAL BASIS: ',/)
       CALL PRTMAT(NUMSTA,HAM,0)
C      **** PRINT HAM IN ORIGINAL BASIS ****
      ELSE
      ENDIF
C
C
C
      ELSE
C     **** READ IN HAMILTONIAN  ****
C
      WRITE(*,470)
  470 FORMAT(' INPUT LAMBDA (0.0<LAMBDA<1.0)')
      READ(*,*)RLAMB
C
c     READ(14,480,REC=2)N,NF,B,LPN,K,NUMSTA,CUTOFF
      READ(14,*)
      READ(14,480)N,NF,B,LPN,K,NUMSTA,CUTOFF
  480 FORMAT(6I6,6X,e24.16)
C
c
c     WRITE(15,510)
c 510 FORMAT(//'   N    NF    B   LPN   K   NUMSTA  '/)
c     WRITE(15,520)N,NF,B,LPN,K,NUMSTA
c 520 FORMAT(6I5)
c     WRITE(15,535)CUTOFF
c 535 FORMAT(/' CUTOFF:  ',D12.6//)
c     WRITE(15,530)RLAMB
c 530 FORMAT(/' LAMBDA:  ',D12.6//)
c
      WRITE(15,511)N,NF,B,K,RLAMB
  511 FORMAT(' N,NF,B,K,LAMBDA:  ',4I5,e24.16)
      WRITE(15,521)LPN
  521 FORMAT(' LPN:  ',I5)
      WRITE(15,536)CUTOFF
  536 FORMAT(' CUTOFF/MASS SQD: ',e24.16)
c
c  ????????????????????????????????????????????????
c     **** read from terminal; needed to know nf from qcdf.ham ****
      rmq(1)=1.d0
      if(nf.gt.1) then
c     **** input quark masses ****
c     **** first quark mass is 1*(1-lambda**2)**.5 ; ****
c     **** the rest are rmq(i=2 to nf) * the first    ****
       do 485 nqk=2,nf
        WRITE(*,475)nqk
  475   FORMAT(' input ratio of quark mass',i3,' to quark mass  1:')
        read(*,*)rmq(nqk)
  485  continue
      else
      endif
C
      WRITE(15,1300)
 1300 FORMAT(' ratio of quark masses to first mass (rmq(nf)): ')
      WRITE(15,1310)(rmq(i),i=1,nf)
 1310 FORMAT(6e24.16)

      WRITE(15,537)NUMSTA
  537 FORMAT(' NUMSTA:  ',I6)
C
c ??????????????????????????????????????????
c     IREC=4
C     **** IREC= RECORD LOCATION OF HAMILT ELEMENT  ****
c    
c     **** read in HNU0(nf) ****
      READ(14,*)
c     **** skip line in qcdf.ham ****
      DO 490 IR1=1,NUMSTA
c     DO 490 IR2=1,IR1
c      READ(14,500,REC=IREC) HNU0(IR1,IR2),HNU(IR1,IR2)
c      READ(14,500) (HNU0(ifl,IR1,IR2),ifl=1,nf),HNU(IR1,IR2)
       READ(14,500) (HNU0(ifl,IR1),ifl=1,nf)
c 500  FORMAT(7(D16.10,1X))       
  500  FORMAT(6(e24.16,1X))       
c     IREC=IREC+1
  490 CONTINUE
c
c     **** read in HNU ****
      READ(14,*)
      DO 493 IR1=1,NUMSTA
      DO 493 IR2=1,IR1
c      READ(14,803) (HNU(IR1,IR2),ir2=1,ir1)
       HNU(IR1,IR2)=0.d0        
c 803  FORMAT(4(D16.10,1X))
  493 CONTINUE
c     **** read in non-zero elements ****
      nzrd=0
  805 read(14,804)(bfr(i,1),bfr(i,2),rbfr(i),i=1,3)
c     WRITE(*,804)(bfr(i,1),bfr(i,2),rbfr(i),i=1,3)
      do 806 ird=1,3
       if(bfr(ird,1).ne.-1) then
        hnu(bfr(ird,1),bfr(ird,2))=rbfr(ird)
       else
        nzrd=1
c       **** end of non-zero elements ****
       endif
  806 continue
      if(nzrd.eq.0)go to 805
  804 FORMAT(3(i6,i6,1x,e24.16,1X))
c
      CLOSE(UNIT=14)
C
C
      ENDIF
C
C
      WRITE(*,1140)
 1140 FORMAT('  COMBINE FREE, INT HAMILT  ')
C     **** COMBINE THE FREE AND INTERACTING HAMILTONIANS ****
C     **** FOR ORTHONORMAL STATES                        ****
c     **** interacting hamiltonian: ****
      RLMSQ=RLAMB*RLAMB
      DO 600 I1=1,NUMSTA
      DO 600 I2=1,I1
       HNU(I1,I2)= RLMSQ*HNU(I1,I2)
c     do 605 ifl=1,nf
c      HNU(I1,I2)= (1.0D0 - RLMSQ)*rmq(ifl)*rmq(ifl)*HNU0(ifl,I1,I2) 
c    >             +HNU(I1,I2)
c 605  continue
c      HNU(I1,I2)= (1.0D0 - RLMSQ)*HNU0(I1,I2) + RLMSQ*HNU(I1,I2)
       HNU(I2,I1)= HNU(I1,I2)
  600 CONTINUE
c
c     **** diag elements: ****
      DO 605 I1=1,NUMSTA
      do 605 ifl=1,nf
       HNU(I1,I1)= (1.0D0 - RLMSQ)*rmq(ifl)*rmq(ifl)*HNU0(ifl,I1) 
     >             +HNU(I1,I1)
  605  continue
C
      IF((IHFILE.NE.1).AND.(IPTHAM.NE.0))THEN
       WRITE(15,110)
  110  FORMAT(/,' HAMILTONIAN IN ORTHONORMAL BASIS: ',/)
       CALL PRTMAT(NUMSTA,HNU,0)
C      **** PRINT HAM IN ORTHON BASIS ****
      ELSE
      ENDIF
C
C     CALL DIAG(HNU,WPR,ZPR,IERPR)
      CALL DIAG(HNU,W,Z,IERR)
C     **** DIAG. HAM MAT ****
      WRITE(*,1150)
 1150 FORMAT('  CALLED SUB DIAG  ')
C
      DO 150 IE=1,NUMSTA
C      W(IE)=K*W(IE)
       W(IE)=K*W(IE)/2.0D0
  150 CONTINUE
C     **** REPLACE HAM EVALS WITH MASS SQ EVALS ****
C
C     CALL PRTEIG(WPR,ZPR,IERPR,1)
      CALL PRTEIG(W,Z,IERR,IVCTMX)
C     **** PRINT RESULTS OF DIAG ****
      WRITE(*,1160)
 1160 FORMAT('  CALLED SUB PRTEIG  ')
      WRITE(*,1170)
 1170 FORMAT('  END OF MAIN PROG   ')
C
C
c     close(unit=7)
 2000 CLOSE(unit=15)
c
c
      STOP
      END
C
C
C
C
      SUBROUTINE QCDSTA
C     **** QCD STATE GENERATOR ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/MATS/ KPX(25001,25),KPXLOC(3,25,100,2)
      COMMON/MATSM/KPM(25001,2),KPMLOC(100,2)
      COMMON/MDAT/ NUMPRM,IX0,IXN
      COMMON/MDATM/NPRMM
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NOBAR/IBAB
      COMMON/PARLIM/LPN
      INTEGER N,NF,B,K
C
C     **** NUMSTA=STATE IDENT NUM  ****
C
      NUMSTA=0
      MXNP =25
C     **** MAX LENGTH OF STATES ****
C
C     N=2
C     B=0
C     K=04
C
C
C     **** CREATE KPX,KPXLOC,KPM,KPMLOC ****
C     KMX=08
      KMX=K
C
      CALL LPNSUB(N,NF,B,K,LNB,LND)
C     **** LNB IS THE MAX NUM OF B-ADJ'S POSSIBLE;  ****
C     **** THIS WILL ALSO BE THE MOST TERMS NEEDED  ****
C     **** IN KPX (LND ARE THE MAX NUM OF D-ADJ'S)  ****
      WRITE(*,200)
  200 FORMAT(' CALLED LPNSUB FROM QCDSTA ')
C
      IF (LNB.GT.MXNP) THEN
       WRITE(15,10)MXNP
   10  FORMAT(' MAX STATE LNGTH,',I5,',TOO SMALL;')
       WRITE(15,12)LNB
   12  FORMAT(' NEED LNGTH OF:',I5)
       STOP
      ELSE
      ENDIF
C
      IF((LPN.NE.0).AND.(LPN.LT.LNB)) THEN
       NPTMX=LPN
      ELSE
       NPTMX=LNB
      ENDIF
C     **** KMX,NPTMX ARE MAX MOM, NUM   ****
C     **** PARTICLES ALLOWED IN KPX,KPM ****
C
      IF(LPN.GT.(LNB+LND)) THEN
       LPN=LNB+LND
      ELSE
      ENDIF
C
      IX0=1
      IXN=2
C
C     CALL PRMX(KMX,NPTMX,N)
c  ???????????????????????????????????????
C     **** FOR FLAVOR, REPLACE N WITH N*NF ****
      NNF=N*NF
      CALL PRMX(KMX,NPTMX,NNF)
      WRITE(*,210)
  210 FORMAT(' CALLED PRMX FROM QCDSTA ')
C
      CALL PRMM(KMX)
      WRITE(*,220)
  220 FORMAT(' CALLED PRMM FROM QCDSTA ')
C
C     **** PRINT KPX,KPXLOC,KPM,KPMLOC ****
c     CALL PRTKPX(KMX,NPTMX,N,NF)
c     CALL PRTKPM(KMX)
C
      call ibarfl
c     **** generate baryon flavor perms ****
      WRITE(*,325)
  325 FORMAT(' CALLED ibarfl FROM QCDSTA ')
C
      call imesfl
c     **** generate meson flavor perms ****
      WRITE(*,335)
  335 FORMAT(' CALLED imesfl FROM QCDSTA ')
C
C
C     LMBF=LMBMX(N,B,K)
C     **** SWITCHED LMBMX FOR LPNSUB WHEN MODIFIED FOR FLAVOR ****
C  ????????????????????????????????????????
c     CALL LPNSUB(N,NF,B,K,LNBX,LNDX)
c  ????? have already computed this quantity (lnb) above ????
      LMBF=LNB-N*B
      DO 100 LMB=0,LMBF
C     **** LMB ('LAMBDA') IS THE NUMBER OF B-ADJ (AND ****
C     **** D-ADJ) BEYOND THE B-ADJ NEEDED TO SATISFY  ****
C     **** BARYON NUMBER B                            ****
       NXBMX=INT(LMB/N)
       IF(IBAB.EQ.0) THEN
        NXBMX=0
       ELSE
       ENDIF
C      WRITE(15,30)LMB,NXBMX
C  30  FORMAT(' LMB,NXBMX: ',3I5)
C
       DO 100 NXB=0,NXBMX
        NBAR=B+2*NXB
C       **** NBAR IS THE NUMBER OF BAR + ANTIBAR      ****
C       **** NXB IS THE NUM OF EXTRA BAR BEYOND B     ****
        NMES=LMB-N*NXB
C       WRITE(15,40)NXB,NBAR,NMES
C  40   FORMAT(' NXB,NBAR,NMES: ',3I5)
C  ????????????????????????????????????????
c       **** compute the min and max momentum carried ****
c       **** by baryons and mesons, respectively      ****
c       KBARMN=N*(B*(B-1)+2*NXB*(B-1+NXB))/2
c       KMESMN=2*INT(NMES/N)*(NMES-N*(INT(NMES/N)+1)/2)
c       **** adjusted below for flavor ****
        NBB = N*(B+NXB)
        NDB = N*NXB
        KBMN=minmom(nbb)
        KDMN=minmom(ndb)
        KMMN=minmom(nmes)
        KBARMN=KBMN+KDMN
        KMESMN=2*KMMN
c       **** double to acct for qks + antiqks in mesons ****
C
        KBARMX=K-KMESMN
        KMESMX=K-KBARMN
        IF(NBAR.EQ.0) THEN
         KMESMN=K
         KMESMX=K
         KBARMN=0
         KBARMX=0
        ELSE
        ENDIF
        IF(NMES.EQ.0) THEN
         KMESMN=0
         KMESMX=0
         KBARMN=K
         KBARMX=K
        ELSE
        ENDIF
C
C       IF((NMES.EQ.0).AND.(NBAR.EQ.0)) GO TO 100
C       IF((NMES.NE.0).OR.(NBAR.NE.0)) THEN
        IF( ((NMES.NE.0).OR.(NBAR.NE.0)) .AND.
     >      (((N*B+LMB).LE.LNB).AND.(LMB.LE.LND)) ) THEN
C
C        WRITE(15,50)KBARMN,KBARMX,KMESMN,KMESMX
C  50    FORMAT(' KBARMN,KBARMX,KMESMN,KMESMX: ',4I5)
C
         NBRB=B+NXB
         NBRD=NXB
C        **** NBRB,KBRB ARE NUM,MOM OF BAR;NBRD,KBRD OF ANTIBAR ****
C
         NOE=IODEV(N)
         NBRBOE=IODEV(NBRB)
         NBRDOE=IODEV(NBRD)
C        **** CHECK IF N, NBRB, NBRD ARE ODD OR EVEN; ****
C        **** USED TO TAKE ADVANT. OF APBC TO LIMIT   ****
C        **** POSSIBLE STATES                         ****
C
C
         DO 105 KBAR=KBARMN,KBARMX
C        **** KBAR IS THE MOM CARRIED BY BAR + ANTIBAR ****
          KMES=K-KBAR
c         KBRBMN=N*NBRB*(NBRB-1)/2
c         KBRDMN=N*NBRD*(NBRD-1)/2
          kbrbmn=minmom(n*nbrb)
          kbrdmn=minmom(n*nbrd)
c         **** modified for apbc, flavor ****
          KBRBMX=KBAR-KBRDMN
          KBRDMX=KBAR-KBRBMN
C
          IF (NBRD.EQ.0) THEN
           KBRBMN=KBAR
           KBRBMX=KBAR
           KBRDMN=0
           KBRDMX=0
          ELSE
          ENDIF
C
         DO 105 KBRB=KBRBMN,KBRBMX
          KBRD=KBAR-KBRB
C
C         WRITE(15,120) NBRB,NBRD,NMES
C         WRITE(15,140) KBRB,KBRD,KMES
C 120     FORMAT(' NBRB,NBRD,NMES: ',3I5)
C 140     FORMAT(' KBRB,KBRD,KMES: ',3I5)
C         WRITE(15,150)
C 150     FORMAT(' -------------------------------')
C
C
          KMESOE=IODEV(KMES)
          KBRBOE=IODEV(KBRB)
          KBRDOE=IODEV(KBRD)
          IF(KMESOE.EQ.1) THEN
           IF( ((NOE.EQ.1).AND.(KBRBOE.EQ.1).AND.(KBRDOE.EQ.1)).
     >         OR.((NOE.EQ.-1).AND.(KBRBOE.EQ.NBRBOE).AND.
     >            (KBRDOE.EQ.NBRDOE)) ) THEN
C         **** TAKE ADVANTAGE OF APBC TO LIMIT POSSIBLE   ****
C         **** COMBINATIONS OF MOMENTA IN BARYONS, MESONS ****
c         **** uses even*odd = even; odd*odd = odd        ****
C
            CALL BRMSGN(N,nf,NBRB,NBRD,NMES,KBRB,KBRD,KMES)
            WRITE(*,230)
  230       FORMAT(' CALLED BRMSGN FROM QCDSTA ')
C
           ELSE
           ENDIF
          ELSE
          ENDIF
C
  105    CONTINUE
C
        ELSE
        ENDIF
C
  100   CONTINUE
C
  110 RETURN
      END
C
C
C
cf    FUNCTION LMBMX0(N,B,K)
C     **** CALC MAX LAMBDA ALLOWED FOR N,B,K ****
C     **** (LAMBDA HERE REFERS TO THE NUMBER ****
C     **** OF QUARKS BEYOND THOSE NEEDED TO  ****
C     **** SATISFY BARYON NUMBER B)          ****
C      **** REPLACED WITH LPNSUB WHEN PROGRAM ****
C      **** MODIFIED TO INCLUDE FLAVOR        ****
cf    IMPLICIT REAL*8 (A-H,O-Z)
cf    INTEGER N,B,K
C     WRITE(15,20)N,B,K
C  20 FORMAT(' LMBMX;N,B,K: ',3I5)
C
cf    LMBMX=INT( (N*K-N*N*(B*B-1)/4)**(0.5) - N*(B-1)/2 )
cf 10 LMXON=INT(LMBMX/N)
cf    KPR=N*B*(B-1)/2 - N*LMXON*(LMXON+1) + LMBMX*(B+2*LMXON)
C     WRITE(15,30)LMBMX,LMXON,KPR
C  30 FORMAT(' LMBMX;LMBMX,LMXON,KPR: ',3I5)
cf    IF (KPR.GT.K) THEN
cf     LMBMX=LMBMX-1
cf     GO TO 10
cf    ELSE
cf     LBX1=LMBMX+1
cf     LMX1=INT(LBX1/N)
cf     KPR1=N*B*(B-1)/2 - N*LMX1*(LMX1+1) + LBX1*(B+2*LMX1)
cf     IF (KPR1.LE.K) THEN
C       WRITE(15,40)
C  40   FORMAT(' ERROR IN FUNCTION LMBMX ')
cf      LMBMX=LMBMX+1
cf      GO TO 10
cf     ELSE
cf     ENDIF
cf    ENDIF
C
cf    LMBMX0=LMBMX
cf    RETURN
cf    END
C
C
C
C
c   ?????????????????????????????????????
      SUBROUTINE BRMSGN(N,nf,NBRB,NBRD,NMES,KBRB,KBRD,KMES)
C     **** GIVEN NUM AND MOM OF BAR, ANTIBAR, AND MESONS, ****
C     **** COMPUTE POSSIBLE STATES                        ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      COMMON/MATS/ KPX(25001,25),KPXLOC(3,25,100,2)
      COMMON/MATSM/KPM(25001,2),KPMLOC(100,2)
      COMMON/MDAT/ NUMPRM,IX0,IXN
      COMMON/MDATM/NPRMM
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NOZERO/IB0,ID0
      COMMON/PARLIM/LPN
      COMMON/NOBAR/IBAB
      COMMON/barfl/mbarfl(1001,6),ibfdim
      COMMON/mesfl/mmesfl(101,2),imfdim
C
      INTEGER NMINB(25),NMAXB(25)
      INTEGER NMIND(25),NMAXD(25)
      INTEGER NMINM(25),NMAXM(25)
      INTEGER NMINBF(25),NMAXBF(25)
      INTEGER NMINDF(25),NMAXDF(25)
      INTEGER NMINMF(25),NMAXMF(25)
c     INTEGER ISTB(6902,25),ISTD(6902,25),ISTM(6902,25)
      INTEGER ISTB(9523,25),ISTD(9523,25),ISTM(9523,25)
      INTEGER ICHKB(25),ICHKD(25),ICHKM(25)
      INTEGER ISTATE(4,25)
C
      istmax=9523
c
c     WRITE(15,600)n,nf,nbrb,nbrd,nmes,kbrb,kbrd,kmes 
c 600 FORMAT(' n,nf,nbrb,nbrd,nmes,kbrb,kbrd,kmes: ',8i4) 
C
c     WRITE(15,730)ibfdim,imfdim
c 730 FORMAT(' ibfdim, imfdim:  ',2i4)
c
      IDIAG=0
C     **** IDIAG =1 GENERATES DIAGNOSTIC          ****
C     **** OUTPUT FOR DEBUGGING; =0 FOR NO OUTPUT ****
c
      nnf=n*nf
C
C     KMX=6
C     NPTMX=10
C
C     IX0=1
C     IXN=2
C
C     CALL PRMX(KMX,NPTMX,N)
C
C     CALL PRMM(KMX)
C
C     CALL PRTKPX(KMX,NPTMX,N,NF)
C
C     CALL PRTKPM(KMX)
C
      DO 20 IN1=1,25
       NMINB(IN1)=0
       NMAXB(IN1)=0
       NMIND(IN1)=0
       NMAXD(IN1)=0
       NMINM(IN1)=0
       NMAXM(IN1)=0
       ICHKB(IN1)=0
       ICHKD(IN1)=0
       ICHKM(IN1)=0
       nmaxbf(in1)=0
       nmaxdf(in1)=0
       nmaxmf(in1)=0
       nminbf(in1)=0
       nmindf(in1)=0
       nminmf(in1)=0
      DO 20 IN2=1,istmax
       ISTB(IN2,IN1)=0
       ISTD(IN2,IN1)=0
       ISTM(IN2,IN1)=0
   20 CONTINUE
C
C     *** PRELIM DATA; TO BE GENERATED BY CALLING ***
C     *** ROUTINE                                 ***
C     NBRB=0
C     KBRB=0
C     NBRD=0
C     KBRD=0
C     NMES=4
C     KMES=5
C     WRITE(15,25)N
C  25 FORMAT(' N: ',I4)
C     WRITE(15,30)NBRB,NBRD,NMES
C  30 FORMAT(' NUM OF BARS, ANTIBARS, MES: ',3I4)
C     WRITE(15,35)KBRB,KBRD,KMES
C  35 FORMAT(' MOM OF BARS, ANTIBARS, MES: ',3I4)
C
C     NUMSTA=0
C
      LNGSTA=N*(NBRB+NBRD)+2*NMES
C
C     **** LPN IS LIMIT IN PARTICLE NUM.;CHECK IF ****
C     **** LNGSTA EXCEEDS IT                      ****
      LPNCHK=1
      IF((LNGSTA.GT.LPN).AND.(LPN.NE.0)) THEN
       LPNCHK=0
      ELSE
      ENDIF
C
C
C     *** DETERMINE IF BARYONS ARE BOSONS OR FERMIONS;  ***
C     *** IE, IS N ODD OR EVEN. IXBAR IS THE FERM INDEX ***
      IF (((-1)**N).GT.0) THEN
c      IXBAR=N
       IXBAR=nnf
      ELSE
       IXBAR=1
      ENDIF
C
C
C
C     **** PERMUTE OVER MOM DISTR FOR EACH BAR,MES ****
      IF(NBRB.EQ.0)THEN
       IBMN=-1
       IBMX=-1
C      **** SET THESE VAL FOR MN,MX SO THAT DO 50...****
C      **** LOOP BELOW EXECUTES THIS LOOP ONCE      ****
      ELSE
       IBMN=KPXLOC(IX0,NBRB,KBRB+1,1)
       IBMX=KPXLOC(IX0,NBRB,KBRB+1,2)
c      WRITE(15,610)ibmn,ibmx
c 610  FORMAT(' ibmn,ibmx: ', 2i5)
      ENDIF
      IF(NBRD.EQ.0)THEN
       IDMN=-1
       IDMX=-1
      ELSE
       IDMN=KPXLOC(IX0,NBRD,KBRD+1,1)
       IDMX=KPXLOC(IX0,NBRD,KBRD+1,2)
c      WRITE(15,620)idmn,idmx
c 620  FORMAT(' idmn,idmx: ', 2i5)
      ENDIF
      IF(NMES.EQ.0)THEN
       IMMN=-1
       IMMX=-1
      ELSE
       IMMN=KPXLOC(IX0,NMES,KMES+1,1)
       IMMX=KPXLOC(IX0,NMES,KMES+1,2)
c      WRITE(15,630)immn,immx
c 630  FORMAT(' immn,immx: ', 2i5)
      ENDIF
C
C     IF((IBMN.EQ.0).OR.(IDMN.EQ.0).OR.(IMMN.EQ.0)) GO TO 50
      IF((IBMN.EQ.0).OR.(IDMN.EQ.0).OR.(IMMN.EQ.0)) GO TO 55
C
      DO 50 IB=IBMN,IBMX
      DO 50 ID=IDMN,IDMX
      DO 50 IM=IMMN,IMMX
C
       IF(NBRB.NE.0)THEN
        DO 60 JB=1,NBRB
         NMINB(JB)=KPXLOC(IXN,N,KPX(IB,JB)+1,1)
         NMAXB(JB)=KPXLOC(IXN,N,KPX(IB,JB)+1,2)
c        WRITE(15,650)ib,jb,ixn,n,kpx(ib,jb)
c 650    FORMAT(' ib,jb,ixn,n,kpx(ib,jb): ',5i4)
         if(nminb(jb).eq.0) go to 49
   60   CONTINUE
       ELSE
       ENDIF
C
       IF(NBRD.NE.0)THEN
        DO 70 JD=1,NBRD
         NMIND(JD)=KPXLOC(IXN,N,KPX(ID,JD)+1,1)
         NMAXD(JD)=KPXLOC(IXN,N,KPX(ID,JD)+1,2)
         if(nmind(jd).eq.0) go to 49
   70   CONTINUE
       ELSE
       ENDIF
C
       IF(NMES.NE.0)THEN
        DO 80 JM=1,NMES
         NMINM(JM)=KPMLOC(KPX(IM,JM)+1,1)
         NMAXM(JM)=KPMLOC(KPX(IM,JM)+1,2)
         if(nminm(jm).eq.0) go to 49
   80   CONTINUE
       ELSE
       ENDIF
C
C
      IF(IDIAG.NE.0) THEN
C     **** PRINT DIAGNOSTIC INFO ****
      WRITE(15,65)NBRB
   65 FORMAT(' NBRB= ',I6)
      WRITE(15,67)
   67 FORMAT(' NMINB: ')
      WRITE(15,69)(NMINB(J1),J1=1,NBRB)
      WRITE(15,68)
   68 FORMAT(' NMAXB: ')
      WRITE(15,69)(NMAXB(J2),J2=1,NBRB)
   69 FORMAT(25I6)
C
      WRITE(15,75)NBRD
   75 FORMAT(' NBRD= ',I6)
      WRITE(15,77)
   77 FORMAT(' NMIND: ')
      WRITE(15,79)(NMIND(J1),J1=1,NBRD)
      WRITE(15,78)
   78 FORMAT(' NMAXD: ')
      WRITE(15,79)(NMAXD(J2),J2=1,NBRD)
   79 FORMAT(25I6)
C
      WRITE(15,85)NMES
   85 FORMAT(' NMES= ',I6)
      WRITE(15,87)
   87 FORMAT(' NMINM: ')
      WRITE(15,89)(NMINM(J1),J1=1,NMES)
      WRITE(15,88)
   88 FORMAT(' NMAXM: ')
      WRITE(15,89)(NMAXM(J2),J2=1,NMES)
   89 FORMAT(25I6)
      ELSE
      ENDIF
C
C
C
C      *** PERMUTE OVER MOM INSIDE EACH BAR,MES ***
       NTRMSB=1
       NTRMSD=1
       NTRMSM=1
       I1=0
       I2=0
       I3=0
C
       IDBL=1
C
c
C
c  ?????????????????????????????????
c      **** adjust configs in nmax, nmin to include perms   ****
c      **** over flavor in individ. baryons and mesons; will****
c      **** iterate over a single index, then decompose     ****
c      **** it into separate indices which identify mom     ****
c      **** and flavor configs                              ****
c      **** total num of configs is (nmax-nmin+1)*ibfdim,etc****
       IF(NBRB.NE.0)THEN
        DO 265 JB=1,NBRB
         nminbf(jb)=(nminb(jb)-1)*ibfdim + 1
         nmaxbf(jb)=nmaxb(jb)*ibfdim 
  265   CONTINUE
c       WRITE(15,710)(nmaxbf(jx),jx=1,nbrb)
c 710   FORMAT(' nmaxbf:  ', 10i5)
       ELSE
       ENDIF
C
       IF(NBRD.NE.0)THEN
        DO 275 JD=1,NBRD
         nmindf(jd)=(nmind(jd)-1)*ibfdim + 1
         nmaxdf(jd)=nmaxd(jd)*ibfdim 
  275   CONTINUE
       ELSE
       ENDIF
C
       IF(NMES.NE.0)THEN
        DO 285 JM=1,NMES
         nminmf(jm)=(nminm(jm)-1)*imfdim + 1
         nmaxmf(jm)=nmaxm(jm)*imfdim 
  285   CONTINUE
       ELSE
       ENDIF
C
c
       IF(NBRB.NE.0)THEN
c       CALL GPRMX(NMINB,NMAXB,NBRB,IXBAR,NTRMSB,ISTB,IDBL)
        CALL GPRBIG(NMINBF,NMAXBF,NBRB,IXBAR,NTRMSB,ISTB,IDBL)
c       WRITE(15,780)
c 780   FORMAT(' istb: ')
c       do 770 i0=1,ntrmsb
c        WRITE(15,760)(istb(i0,i1),i1=1,nbrb)
c 760    FORMAT(10i5)
c 770   continue
       ELSE
       ENDIF
C
       IF(NBRD.NE.0)THEN
c       CALL GPRMX(NMIND,NMAXD,NBRD,IXBAR,NTRMSD,ISTD,IDBL)
        CALL GPRBIG(NMINDF,NMAXDF,NBRD,IXBAR,NTRMSD,ISTD,IDBL)
       ELSE
       ENDIF
C
       IF(NMES.NE.0)THEN
c       CALL GPRMX(NMINM,NMAXM,NMES,N,NTRMSM,ISTM,IDBL)
        CALL GPRBIG(NMINMF,NMAXMF,NMES,NNF,NTRMSM,ISTM,IDBL)
c       WRITE(15,980)
c 980   FORMAT(' istm: ')
c       do 970 i0=1,ntrmsm
c        WRITE(15,960)(istm(i0,i1),i1=1,nmes)
c 960    FORMAT(10i5)
c 970   continue
       ELSE
       ENDIF
C
       DO 100 L1=1,NTRMSB
       DO 100 L2=1,NTRMSD
       DO 100 L3=1,NTRMSM
C
       IF(NBRB.NE.0)THEN
        DO 110 J1=1,NBRB
c        **** factor out number of flavor configs: ****
         ibk=int((istb(l1,j1)-1)/ibfdim)+1
c        WRITE(15,700)j1,l1,istb(l1,j1),nminb(j1),ibk
c 700    FORMAT(' j1,l1,istb(l1,j1),nminb(j1),ibk: ',5i4)
         ibf=istb(l1,j1)-(ibk-1)*ibfdim
        DO 110 K1=1,N
         I1=N*(J1-1)+K1
         ISTATE(1,I1)=1
         ISTATE(2,I1)=0
c        ISTATE(3,I1)=KPX(ISTB(L1,J1),K1)
         ISTATE(3,I1)=KPX(ibk,K1)
c        ISTATE(4,I1)=J1
         ISTATE(4,I1)=mbarfl(ibf,k1)
c        WRITE(15,500) ibf,k1,mbarfl(ibf,k1) 
c 500    FORMAT(' ibf, k1, mbarfl(ibf,k1) :  ', 3i6)
c        **** flavor index ****
  110   CONTINUE
       ELSE
       ENDIF
C
       IF(NBRD.NE.0)THEN
        DO 120 J2=1,NBRD
c        **** factor out number of flavor configs: ****
         idk=int((istd(l2,j2)-1)/ibfdim)+1
         idf=istd(l2,j2)-(idk-1)*ibfdim
        DO 120 K2=1,N
         I2=N*NBRB+N*(J2-1)+K2
         ISTATE(1,I2)=0
         ISTATE(2,I2)=1
c        ISTATE(3,I2)=KPX(ISTD(L2,J2),K2)
         ISTATE(3,I2)=KPX(idk,K2)
c        ISTATE(4,I2)=NBRB+J2
         ISTATE(4,I2)=mbarfl(idf,k2)
c        WRITE(15,505) idf,k2,mbarfl(idf,k2) 
c 505    FORMAT(' idf, k2, mbarfl(idf,k2) :  ', 3i6)
c        **** flavor index ****
  120   CONTINUE
       ELSE
       ENDIF
C
       IF(NMES.NE.0)THEN
        DO 130 J3=1,NMES
c        **** factor out number of flavor configs: ****
         imk=int((istm(l3,j3)-1)/imfdim)+1
         imf=istm(l3,j3)-(imk-1)*imfdim
        DO 130 K3=1,2
         I3=N*(NBRB+NBRD)+2*(J3-1)+K3
         ISTATE(1,I3)=-K3+2
         ISTATE(2,I3)=K3-1
c        ISTATE(3,I3)=KPM(ISTM(L3,J3),K3)
         ISTATE(3,I3)=KPM(imk,K3)
c        ISTATE(4,I3)=NBRB+NBRD+J3
         ISTATE(4,I3)=mmesfl(imf,k3)
c        WRITE(15,510) imf,k3,mmesfl(imf,k3) 
c 510    FORMAT(' imf, k3, mmesfl(imf,k3) :  ', 3i6)
c        **** flavor index ****
  130   CONTINUE
       ELSE
       ENDIF
C
C 100  CONTINUE
C
C
c      WRITE(15,133)
c 133  FORMAT(' ISTATE: ')
c      DO 137 M1=1,4
c      WRITE(15,135)(ISTATE(M1,M2),M2=1,LNGSTA)
c 135  FORMAT(25I3)
c 137  CONTINUE
C
C
C
C      **** CHECK IF HAVE ANY EVEN MOMENTA IN STATE;  ****
C      **** IF SO, IODOK=0;DISCARD STATE              ****
       IODOK=1
       CALL ODDCHK(IODOK,LNGSTA,ISTATE)
C
C
C      **** CHECK IF B0  OR D0 ALLOWED; IF NOT, ****
C      **** SEE IF ANY OCCUR IN NEW STATE       ****
C      IZRCHK=1
C
C      IF(ID0.EQ.0) THEN
C       DO 155 I0=1,LNGSTA
C        IF((ISTATE(3,I0).EQ.0).AND.(ISTATE(1,I0).EQ.0)) THEN
C         IZRCHK=0
C        ELSE
C        ENDIF
C 155   CONTINUE
C      ELSE
C      ENDIF
C
C      IF(IB0.EQ.0) THEN
C       DO 157 J0=1,LNGSTA
C        IF((ISTATE(3,J0).EQ.0).AND.(ISTATE(2,J0).EQ.0)) THEN
C         IZRCHK=0
C        ELSE
C        ENDIF
C 157   CONTINUE
C      ELSE
C      ENDIF
C
C
C
c  ???????????????????????????????????????
C      **** CHECK IF HAVE MORE THAN N B-ADJ OR N D-ADJ            ****
C      **** WITH THE SAME MOM and flavor; IF SO, DISCARD STATE    ****
       NUMBS=N*NBRB+NMES
       NUMDS=N*NBRD+NMES
       IBCH=0
       IDCH=0
       JXOK=1
       DO 150 MCH=1,LNGSTA
        IF (ISTATE(1,MCH).EQ.1) THEN
         IBCH=IBCH+1
c        ICHKB(IBCH)=ISTATE(3,MCH)
         momb = (istate(3,mch) + 1)/2 
c         **** map odd momenta onto even integers ****
c         **** (assuming apbc's) for convenience  ****
c         **** in sub stachk                      ****
         ICHKB(IBCH)=istate(4,mch)+(momb-1)*nf
c        **** map both momentum and flavor into a single ****
c        **** number from 1 to momb*nf to use in stachk  ****
        ELSE
         IDCH=IDCH+1
c        ICHKD(IDCH)=ISTATE(3,MCH)
         momd = (istate(3,mch) + 1)/2
         ICHKD(IDCH)=istate(4,mch)+(momd-1)*nf
        ENDIF
  150  CONTINUE
       JNDX=N
       CALL STACHK(NUMBS,JNDX,ICHKB,JXOK)
       IF(JXOK.NE.0) THEN
        CALL STACHK(NUMDS,JNDX,ICHKD,JXOK)
       ELSE
       ENDIF
C
C
C      **** CHECK IF STATE OBEYS PAULI-VILLARS CUT-OFF  ****
C      **** RESTRICTION; IF NOT DISCARD STATE           ****
       CALL PAVILL(ISTATE,LNGSTA,LCTCHK)
C
c
C     **** check if state has correct flavor qm numbers ****
C     **** if so, ifchk =1; else 0                       ****
       call flvchk(ISTATE,LNGSTA,ifchk)
C
C
C      WRITE(15,153)JXOK,IZRCHK
C 153  FORMAT(' JXOK(FM STAT); IZRCHK(0 MOM STATES):',2I3)
C
C      **** STORE NEW STATE IN MSTATE; NEW ****
C      **** STATE INFO IN MSTINF           ****
C      IF((JXOK.NE.0).AND.(IZRCHK.NE.0).AND.
C      IF((JXOK.NE.0).AND.
C    >    (LPNCHK.NE.0).AND.(IODOK.NE.0)) THEN
C      IF((JXOK.NE.0).AND.(LCTCHK.NE.0).AND.
C    >    (LPNCHK.NE.0).AND.(IODOK.NE.0)) THEN
       IF((JXOK.NE.0).AND.(LCTCHK.NE.0).AND.(ifchk.NE.0).AND.
     >    (LPNCHK.NE.0).AND.(IODOK.NE.0)) THEN
        NUMSTA=NUMSTA+1
        MSTINF(NUMSTA,1)=4*(NUMSTA-1)+1
C                        * LOC IN MSTATE WHERE BEGINS *
        MSTINF(NUMSTA,2)=LNGSTA
C                        * LENGTH OF STATE *
        MSTINF(NUMSTA,3)=NMES
C                        * NUM OF MESONS *
        MSTINF(NUMSTA,4)=NBRB
C                        * NUM OF BARYONS *
        MSTINF(NUMSTA,5)=NBRD
C                        * NUM OF ANTIBARS *
        MSTINF(NUMSTA,6)=KMES
C                        * MOM OF MESONS *
        MSTINF(NUMSTA,7)=KBRB
C                        * MOM OF BARS *
        MSTINF(NUMSTA,8)=KBRD
C                        * MOM OF ANTBARS *
C
        DO 165 IW=1,4
        DO 165 IZ=1,LNGSTA
         JW=MSTINF(NUMSTA,1)-1+IW
         MSTATE(JW,IZ)=ISTATE(IW,IZ)
c       **** the fourth row is now used for the ****
c       **** flavor index                       ****
  165   CONTINUE
C
C
c       **** PRINT NEW STATE  ****
c       IF (JXOK.NE.0) THEN
c        WRITE(15,160)NUMSTA
c 160    FORMAT(' NUMSTA:  ',I4)
c        DO 170 ISP=1,4
c         WRITE(15,180)(ISTATE(ISP,JSP),JSP=1,LNGSTA)
c 180     FORMAT(25I4)
c 170    CONTINUE
c       ELSE
c       ENDIF
C
       ELSE
       ENDIF
C
  100  CONTINUE
C
   49  continue
   50 CONTINUE
C
   55 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE PRNTST
C     **** PRINT STATES STORED IN MSTATE AND ****
C     **** INFO STORED IN MSTINF             ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NST/NUMSTA
C
      WRITE(15,10)NUMSTA
   10 FORMAT('TOTAL NUM STATES:  ',I6)
      DO 20 IST=1,NUMSTA
       WRITE(15,30)IST
   30  FORMAT('STATE NUMBER:',I6)
       WRITE(15,40)
   40  FORMAT(' LNGTH,NMES,NBRB,NBRD,KMES,KBRB,KBRD')
       WRITE(15,50)(MSTINF(IST,JST),JST=2,8)
   50  FORMAT(7I5)
C
       WRITE(15,60)
   60  FORMAT(' ------- ')
       DO 70 ISU=MSTINF(IST,1),MSTINF(IST,1)+3
        WRITE(15,80)(MSTATE(ISU,JSP),JSP=1,MSTINF(IST,2))
   80   FORMAT(25I6)
   70  CONTINUE
       WRITE(15,90)
   90  FORMAT(' --------------------------------------- ')
   20 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE GPRMX(NMIN,NMAX,LNGTH,INDX,NPRMG,IST,IDBL)
C     **** GPRMX                             ****
C     **** IDENTICAL TO GPRM BUT ALLOWS FOR  ****
C     **** A FERMIONIC INDEX, "INDX"         ****
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NUTERM(25)
      INTEGER NMIN(25),NMAX(25)
      INTEGER IST(6902,25)
C
      istmax=6902
c
c     WRITE(15,5)LNGTH
c   5 FORMAT(' LNGTH= ',I4)
c     WRITE(15,7)
c   7 FORMAT(' NMIN: ')
c     WRITE(15,10)(NMIN(J1),J1=1,LNGTH)
c     WRITE(15,9)
c   9 FORMAT(' NMAX: ')
c     WRITE(15,10)(NMAX(J2),J2=1,LNGTH)
c  10 FORMAT(25I4)
C
      NPRMG=1
      IDIREC=1
      JJ=1
C
      DO 20 MM=1,25
       NUTERM(MM)=0
      DO 20 NN=1,6902
       IST(NN,MM)=0
   20 CONTINUE
C
C
c     WRITE(15,30)INDX
c  30 FORMAT(' FERMIONIC INDEX "INDX"= ',I5)
C
      DO 50 LL=1,LNGTH
       NUTERM(LL)=NMIN(LL)
   50 CONTINUE
C
  200 IF((IDIREC.EQ.1).OR.(NUTERM(JJ).LT.NMAX(JJ))) THEN
         IF(IDIREC.NE.1) NUTERM(JJ)=NUTERM(JJ)+1
C        **** PREVENT DOUBLE COUNTING IF IDBL=1 ********
         IF((IDBL.EQ.1).AND.(JJ.NE.1).AND.(NMIN(JJ).EQ.NMIN(JJ-1))
     >      .AND.((NUTERM(JJ)-1).LT.NUTERM(JJ-1))) THEN
          NUTERM(JJ)=NUTERM(JJ-1)
         ELSE
         ENDIF
C        ********
         JJ=JJ+1
         IDIREC=1
         IF(JJ.GE.LNGTH) THEN
C        **** PREVENT DOUBLE COUNTING IF IDBL=1 ********
          IF((IDBL.EQ.1).AND.(LNGTH.NE.1).AND.(NMIN(LNGTH).EQ.
     >       NMIN(LNGTH-1)).AND.((NUTERM(LNGTH)-1).LT.
     >       NUTERM(LNGTH-1))) THEN
           NUTERM(LNGTH)=NUTERM(LNGTH-1)
          ELSE
          ENDIF
C        ********
          DO 290 II=NUTERM(LNGTH),NMAX(LNGTH)
           NUTERM(LNGTH)=II
           IF(INDX.NE.0) THEN
            CALL STACHK(LNGTH,INDX,NUTERM,IXOK)
           ELSE
            IXOK=1
           ENDIF
           IF(IXOK.NE.0) THEN
c           WRITE(15,270)NPRMG
c 270       FORMAT(' NPRMG= ',I4)
c           WRITE(15,285)(NUTERM(IM),IM=1,LNGTH)
c 285       FORMAT(10I4)
            DO 300 J=1,LNGTH
             IST(NPRMG,J)=NUTERM(J)
  300       CONTINUE
            NPRMG=NPRMG+1
            if(nprmg.eq.istmax) then
             WRITE(15,500)
  500        FORMAT(' nprmg exceeds istmax in sub gprmx ')
             stop
            else
            endif
           ELSE
           ENDIF
  290     CONTINUE
            NUTERM(LNGTH)=NMIN(LNGTH)
            JJ=LNGTH-1
            IF(JJ.EQ.0) GO TO 100
            IDIREC=-1
         ELSE
C           NUTERM(JJ)=NUTERM(JJ-1)
         ENDIF
      ELSE
         NUTERM(JJ)=NMIN(JJ)
         JJ=JJ-1
         IF(JJ.EQ.0) GO TO 100
         IDIREC=-1
      ENDIF
      GO TO 200
  100 CONTINUE
      NPRMG=NPRMG-1
C
      RETURN
      END
C
C
C
C
C
C
      SUBROUTINE GPRBIG(NMIN,NMAX,LNGTH,INDX,NPRMG,IST,IDBL)
C     **** SAME SUBROUTINE AS GPRMX BUT      ****
C     **** IST DIMENSION IS LARGER; SET TO   ****
C     **** MATCH IEPSPRM IN SUB EPS          ****
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NUTERM(25)
      INTEGER NMIN(25),NMAX(25)
      INTEGER IST(9523,25)
C
      istmax=9523
c
C     WRITE(15,5)LNGTH
C   5 FORMAT(' LNGTH= ',I4)
C     WRITE(15,7)
C   7 FORMAT(' NMIN: ')
C     WRITE(15,10)(NMIN(J1),J1=1,LNGTH)
C     WRITE(15,9)
C   9 FORMAT(' NMAX: ')
C     WRITE(15,10)(NMAX(J2),J2=1,LNGTH)
C  10 FORMAT(25I4)
C
      NPRMG=1
      IDIREC=1
      JJ=1
C
      DO 20 MM=1,25
       NUTERM(MM)=0
      DO 20 NN=1,6902
       IST(NN,MM)=0
   20 CONTINUE
C
C
C     WRITE(15,30)INDX
C  30 FORMAT(' FERMIONIC INDEX "INDX"= ',I5)
C
      DO 50 LL=1,LNGTH
       NUTERM(LL)=NMIN(LL)
   50 CONTINUE
C
  200 IF((IDIREC.EQ.1).OR.(NUTERM(JJ).LT.NMAX(JJ))) THEN
         IF(IDIREC.NE.1) NUTERM(JJ)=NUTERM(JJ)+1
C        **** PREVENT DOUBLE COUNTING IF IDBL=1 ********
         IF((IDBL.EQ.1).AND.(JJ.NE.1).AND.(NMIN(JJ).EQ.NMIN(JJ-1))
     >      .AND.((NUTERM(JJ)-1).LT.NUTERM(JJ-1))) THEN
          NUTERM(JJ)=NUTERM(JJ-1)
         ELSE
         ENDIF
C        ********
         JJ=JJ+1
         IDIREC=1
         IF(JJ.GE.LNGTH) THEN
C        **** PREVENT DOUBLE COUNTING IF IDBL=1 ********
          IF((IDBL.EQ.1).AND.(LNGTH.NE.1).AND.(NMIN(LNGTH).EQ.
     >       NMIN(LNGTH-1)).AND.((NUTERM(LNGTH)-1).LT.
     >       NUTERM(LNGTH-1))) THEN
           NUTERM(LNGTH)=NUTERM(LNGTH-1)
          ELSE
          ENDIF
C        ********
          DO 290 II=NUTERM(LNGTH),NMAX(LNGTH)
           NUTERM(LNGTH)=II
           IF(INDX.NE.0) THEN
            CALL STACHK(LNGTH,INDX,NUTERM,IXOK)
           ELSE
            IXOK=1
           ENDIF
           IF(IXOK.NE.0) THEN
C           WRITE(15,270)NPRMG
C 270       FORMAT(' NPRMG= ',I4)
C           WRITE(15,285)(NUTERM(IM),IM=1,LNGTH)
C 285       FORMAT(10I4)
            DO 300 J=1,LNGTH
             IST(NPRMG,J)=NUTERM(J)
  300       CONTINUE
            NPRMG=NPRMG+1
            if(nprmg.eq.istmax) then
             WRITE(15,500)
  500        FORMAT(' nprmg exceeds istmax in sub gprbig ')
             stop
            else
            endif
           ELSE
           ENDIF
  290     CONTINUE
            NUTERM(LNGTH)=NMIN(LNGTH)
            JJ=LNGTH-1
            IF(JJ.EQ.0) GO TO 100
            IDIREC=-1
         ELSE
C           NUTERM(JJ)=NUTERM(JJ-1)
         ENDIF
      ELSE
         NUTERM(JJ)=NMIN(JJ)
         JJ=JJ-1
         IF(JJ.EQ.0) GO TO 100
         IDIREC=-1
      ENDIF
      GO TO 200
  100 CONTINUE
      NPRMG=NPRMG-1
C
      RETURN
      END
C
C
C
C
C
C
C
      SUBROUTINE PRMX(KMX,NPTMX,N)
C
C     ******* GENERATE MOM DISTR. WITH FERMIONIC STATISTICS ******
C     ******* OF INDEX = 'INDX' , WHERE INDEX GIVES THE NUM ******
C     ******* OF DEG OF FREEDOM APART FROM MOMENTUM.        ******
C     ******* MAX KMX=50, MAX NPTMX=25                      ******
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/A/KPTCL,NPTCL
      COMMON/B/NUTERM(25)
      COMMON/MATS/ KPX(25001,25),KPXLOC(3,25,100,2)
      COMMON/MDAT/ NUMPRM,IX0,IXN
      COMMON/NOZERO/IB0,ID0
C
C
      NUMPRM=1
      INMPRM=1
      LKPXMX=25001
C     **** MAX NUM OF PERMS IN KPX ****
C
      MXKMX=100
C     **** MAX TOTAL MOM OF PERM IN KPX ****
      IF (KMX.GE.MXKMX) THEN
       WRITE(15,25)
   25  FORMAT(' KMX EXCEEDS MXKMX IN SUB PRMX ')
      ELSE
      ENDIF
C
C     IX0=1
C     IXN=2
C
      DO 40 L1=1,LKPXMX
      DO 40 L2=1,25
       KPX(L1,L2)=0
   40 CONTINUE
      DO 50 L0=IX0,IXN
      DO 50 L3=1,25
      DO 50 L4=1,MXKMX
      DO 50 L5=1,2
       KPXLOC(L0,L3,L4,L5)=0
   50 CONTINUE
C
C
      DO 510 IX=IX0,IXN
       INDX=IX-1
       IF(IX.EQ.IXN) INDX=N
      DO 510 NPTCL=1,NPTMX
C
      IF( (IB0.EQ.0).AND.(ID0.EQ.0) ) THEN
C     **** NO ZERO MOMENTA ALLOWED ****
       KPTCLI=NPTCL
       KPTCLF=KMX
      ELSE
       KPTCLI=0
       KPTCLF=KMX
      ENDIF
      DO 510 KPTCL=KPTCLI,KPTCLF
C
      DO 45 L1=1,25
       IF( (IB0.EQ.0).AND.(ID0.EQ.0) ) THEN
        NUTERM(L1)=1
       ELSE
        NUTERM(L1)=0
       ENDIF
   45 CONTINUE
      KPTPR=KPTCL+1
C
      INMPRM=NUMPRM
C
      IDIREC=1
      JJ=1
  200 IF((IDIREC.EQ.1).OR.(NUTERM(JJ).LT.KMAX(JJ))) THEN
         IF(IDIREC.NE.1) NUTERM(JJ)=NUTERM(JJ)+1
          JJ=JJ+1
          IDIREC=1
          IF(JJ.GE.NPTCL) THEN
           NUTERM(NPTCL)=KPTCL-KLEFT(NPTCL)
           IF(INDX.NE.0) THEN
            CALL STACHK(NPTCL,INDX,NUTERM,IXOK)
           ELSE
            IXOK=1
           ENDIF
           IF (IXOK.NE.0) THEN
C           WRITE(15,285)(NUTERM(IM),IM=1,NPTCL)
C 285       FORMAT(10I4)
            DO 300 J=1,NPTCL
              KPX(NUMPRM,J)=NUTERM(J)
  300       CONTINUE
            KPXLOC(IX,NPTCL,KPTPR,1)=INMPRM
C            *PLACE ABOVE STMT HERE IN CASE NO NEW*
C            *PERMS GENERATED*
            KPXLOC(IX,NPTCL,KPTPR,2)=NUMPRM
            NUMPRM=NUMPRM+1
           ELSE
           ENDIF
            JJ=JJ-1
            IDIREC=-1
          ELSE
           NUTERM(JJ)=NUTERM(JJ-1)
          ENDIF
      ELSE
         JJ=JJ-1
         IF(JJ.EQ.0) GO TO 100
         IDIREC=-1
      ENDIF
      GO TO 200
  100 CONTINUE
  510 CONTINUE
      NUMPRM=NUMPRM-1
C
      IF (NUMPRM.GE.LKPXMX) THEN
       WRITE(15,520)
  520  FORMAT(' $$$$ OVERFLOW IN KPX $$$$')
       STOP
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE PRMM(KMX)
C     **** GENERATE THE MOMENTA DISTRIBUTIONS FOR ****
C     **** INDIVIDUAL MESONS                      ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MATSM/KPM(25001,2),KPMLOC(100,2)
      COMMON/MDATM/NPRMM
      COMMON/NOZERO/IB0,ID0
C
      NPRMM=0
      LKPMMX=25001
C
      IF( (IB0.EQ.0).AND.(ID0.EQ.0) ) THEN
C     **** NO ZERO MOMENTA ALLOWED ****
       K1IN=2
       K1FI=KMX
      ELSE
       K1IN=0
       K1FI=KMX
      ENDIF
      DO 50 K1=K1IN,K1FI
       KPMLOC(K1+1,1)=NPRMM+1
       IF( (IB0.EQ.0).AND.(ID0.EQ.0) ) THEN
        K2IN=1
        K2FI=K1-1
       ELSE
        K2IN=0
        K2FI=K1
       ENDIF
       DO 100 K2=1,K1-1
        NPRMM=NPRMM+1
        KPM(NPRMM,1)=K2
        KPM(NPRMM,2)=K1-K2
  100  CONTINUE
       KPMLOC(K1+1,2)=NPRMM
   50 CONTINUE
C
C     WRITE(15,150)NPRMM
C 150 FORMAT(' NPRMM= ',I5)
C
      IF (NPRMM.GE.LKPMMX) THEN
       WRITE(15,200)
  200  FORMAT(' $$$$ OVERFLOW IN KPM $$$$')
       STOP
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE PRTKPX(KMX,NPTMX,N,NF)
C     ****  PRINT CONTENTS OF KPX, KPXLOC     ****
C     ****                                    ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MATS/ KPX(25001,25),KPXLOC(3,25,100,2)
      COMMON/MDAT/ NUMPRM,IX0,IXN
C
      WRITE(15,310)N,NF,KMX,NPTMX
  310 FORMAT(' N= ',I3,' NF= ',I3,'  KMX= ',I5,'  NPTMX= ',I5)
      WRITE(15,340)
  340 FORMAT(' KPX: ')
      DO 425 KK=1,NUMPRM
      WRITE(15,400)KK,(KPX(KK,J4),J4=1,NPTMX)
  400 FORMAT(I5,' : ',25I6)
  425 CONTINUE
      WRITE(15,430)
  430 FORMAT(' KPXLOC: ')
      WRITE(15,435)
  435 FORMAT(' K:')
      WRITE(15,440)0,1,2,3,4,5,6,7,8
  440 FORMAT(25I6)
      WRITE(15,445)
  445 FORMAT('   -   -   -   -   -   -   -   -   -')
      DO 480 J0=IX0,IXN
       WRITE(15,450)J0
  450  FORMAT(' IX= ',I6)
      DO 480 J1=1,NPTMX
       WRITE(15,460)J1
  460  FORMAT(' NPTCL= ',I6)
      DO 480 J3=1,2
      WRITE(15,475)(KPXLOC(J0,J1,J2,J3),J2=1,KMX+1)
  475 FORMAT(25I6)
  480 CONTINUE
c
      RETURN
      END
C
C
C
      SUBROUTINE PRTKPM(KMX)
C     ****  PRINT CONTENTS OF KPM, KPMLOC     ****
C     ****                                    ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MATSM/ KPM(25001,2),KPMLOC(100,2)
      COMMON/MDATM/ NPRMM
C
      WRITE(15,310)KMX
  310 FORMAT('  KMX= ',I6)
      WRITE(15,340)
  340 FORMAT(' KPM: ')
      DO 425 KK=1,NPRMM
      WRITE(15,400)KK,(KPM(KK,J4),J4=1,2)
  400 FORMAT(I5,' : ',25I6)
  425 CONTINUE
      WRITE(15,430)
  430 FORMAT(' KPMLOC: ')
      WRITE(15,435)
  435 FORMAT(' K:')
      WRITE(15,440)0,1,2,3,4,5,6,7,8
  440 FORMAT(25I6)
      WRITE(15,445)
  445 FORMAT('   -   -   -   -   -   -   -   -   -')
      DO 480 J3=1,2
      WRITE(15,475)(KPMLOC(J2,J3),J2=1,KMX+1)
  475 FORMAT(25I6)
  480 CONTINUE
C
      RETURN
      END
C
C
C
      FUNCTION KMAX(MM)
C     * CALC THE MAX MOM ELEMENT MM OF NUTERM CAN CARRY *
C     * GIVEN NO. OF TERMS TO RT.,TOTAL MOM. TO LT .    *
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/A/KPTCL,NPTCL
      NRT=NPTCL-MM
      KMAX=(KPTCL-KLEFT(MM))/(NRT+1)
      RETURN
      END
C
C
C
      FUNCTION KLEFT(MM)
C     * COMPUTE MOM CARRIED BY ELEMS TO LT OF MM IN NUTERM *
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/B/NUTERM(25)
      KLEFT=0
      IF (MM.EQ.1) RETURN
      DO 60 I=1,MM-1
      KLEFT=KLEFT+NUTERM(I)
   60 CONTINUE
      RETURN
      END
C
C
C
C
c  ???????????????????????????????????????
      SUBROUTINE IBARFL
C     **** generate flavor perms for baryons ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/barfl/mbarfl(1001,6),ibfdim
      INTEGER N,NF,B,K 
      INTEGER NUTERM(6)
      INTEGER NMIN(6),NMAX(6)
c     **** max num colors = 6 **** 
c     **** max num perms = 1001 **** 
C
c     **** essentially a copy of sub gprmx: ****
C
      ncolmx=6
      nprmmx=1001
c     idbl=0
c     indx=0
C
      lngth=n
      NPRMG=1
      IDIREC=1
      JJ=1
C
c     DO 20 MM=1,ncolmx
      DO 20 MM=1,n
c      NUTERM(MM)=0
      DO 20 NN=1,nprmmx
       mbarfl(NN,MM)=0
   20 CONTINUE
C
      do 30 infl=1,n
       nmin(infl)=1
       nmax(infl)=nf
       nuterm(infl)=1
   30 continue
C
C
c     WRITE(15,500)
c 500 FORMAT(' mbarfl: ')
  200 IF((IDIREC.EQ.1).OR.(NUTERM(JJ).LT.NMAX(JJ))) THEN
         IF(IDIREC.NE.1) NUTERM(JJ)=NUTERM(JJ)+1
         JJ=JJ+1
         IDIREC=1
         IF(JJ.GE.LNGTH) THEN
          DO 290 II=NUTERM(LNGTH),NMAX(LNGTH)
           NUTERM(LNGTH)=II
            DO 300 J=1,LNGTH
             mbarfl(NPRMG,J)=NUTERM(J)
  300       CONTINUE
c           WRITE(15,400)(mbarfl(nprmg,jpr),jpr=1,lngth)
c 400       FORMAT(6i4)
            if(nprmg.eq.nprmmx) then
             WRITE(15,303)
  303        FORMAT(' nprmg exceeds nprmmx in sub ibarfl ')
             stop
            else
            endif
            NPRMG=NPRMG+1
  290     CONTINUE
            NUTERM(LNGTH)=NMIN(LNGTH)
            JJ=LNGTH-1
            IF(JJ.EQ.0) GO TO 100
            IDIREC=-1
         ELSE
         ENDIF
      ELSE
         NUTERM(JJ)=NMIN(JJ)
         JJ=JJ-1
         IF(JJ.EQ.0) GO TO 100
         IDIREC=-1
      ENDIF
      GO TO 200
  100 CONTINUE
      NPRMG=NPRMG-1
C
      ibfdim=nprmg
c     WRITE(15,410)ibfdim
c 410 FORMAT(' ibfdim: ',i6)
c
      RETURN
      END
C
C
C
c  ???????????????????????????????????????
      SUBROUTINE IMESFL
C     **** generate flavor perms for mesons ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/mesfl/mmesfl(101,2),imfdim
      INTEGER N,NF,B,K 
c
c     **** max num perms = 101 **** 
      nprmmx=101
C
      nprm=1
c
      do 10 i=1,nf
      do 10 j=1,nf
       mmesfl(nprm,1)=i
       mmesfl(nprm,2)=j
       if(nprm.eq.nprmmx) then
        WRITE(15,20)
   20   FORMAT(' nprm exceeds nprmmx in sub imesfl ')
        stop
       else
       endif
       nprm=nprm+1
   10 continue
C
      imfdim=nprm-1
C
c     WRITE(15,50)imfdim
c  50 FORMAT(' imfdim: ',i5)
c
c     do 30 id=1,imfdim
c      WRITE(15,40)(mmesfl(id,kd),kd=1,2)
c  40  FORMAT(2i4)
c  30 continue
c
      return
      end
C
C
c
c
c
C
      SUBROUTINE STACHK(NPTCL,INDX,NUTERM,IXOK)
C     * CHECK IF STATE SATISFIES FERMI STATS W/ INDEX='INDX'*
C     * IF NOT SATISFIED, IXOK=0; ELSE 1 *
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NUTERM(25)
      INTEGER MCOUNT(25001),NEGMCT(25001)
C
      DO 25 MM=1,25001
       MCOUNT(MM)=0
       NEGMCT(MM)=0
   25 CONTINUE
C
      IXOK=1
      DO 35 NN=1,NPTCL
       IF (NUTERM(NN).GE.0) THEN
        KP1=NUTERM(NN)+1
        MCOUNT(KP1)=MCOUNT(KP1)+1
        IF (MCOUNT(KP1).GT.INDX) IXOK=0
       ELSE
        NEGKP1=-NUTERM(NN) +1
        NEGMCT(NEGKP1)=NEGMCT(NEGKP1)+1
        IF (NEGMCT(NEGKP1).GT.INDX) IXOK=0
       ENDIF
   35 CONTINUE
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE ODDCHK(IODOK,LNGSTA,ISTATE)
C     * CHECK IF STATE CONTAINS ONLY ODD MOMENTA   *
C     * IF NOT SATISFIED, IODOK=0; ELSE 1 *
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISTATE(4,25)
C
      DO 10 ND=1,LNGSTA
       X=DFLOAT( ISTATE(3,ND) )
       Y=X-2.D0*DFLOAT(INT(X/2.D0))
       IF ( Y.LT.(.5D0)) THEN
        IODOK=0
        RETURN
       ELSE
       ENDIF
   10 CONTINUE
C
      RETURN
      END
C
C
C
C
      FUNCTION IODEV(I)
C     **** IF I EVEN, IODEV=1; IF ODD, -1  ****
C     **** (THIS IS USED WITH APBC)        ****
      IMPLICIT REAL*8 (A-H,O-Z)
C
      X=DFLOAT(I)
      Y=X-2.D0*DFLOAT(INT(X/2.D0))
      IF ( Y.LT.(.5D0)) THEN
       IODEV = 1
      ELSE
       IODEV =-1
      ENDIF
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE CLRDIS(NORH)
c ?????????????????????????????????????
C     **** NORH: 0 IF COMPUTING NORM MAT; 1 IF COMP HAM ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NST/NUMSTA
      COMMON/LRST/ISTRT(4,25),ISTLT(4,25)
      COMMON/STDATA/ISNPR,ISNPL,LSTRT,LSTLT
      COMMON/HM/HAM0(03,6902,6902),HAM(6902,6902)
      COMMON/NM/HNORM(6902,6902)
      COMMON/XBD/NBAJX,NBX,NDAJX,NDX
      COMMON/NBD/NBR,NDR,NBL,NDL
      COMMON/PRTOPS/INTPRT
C     **** ISNPR,L=SIGN OF STATES DUE ****
C     **** TO EPS COLOR PERMS         ****
C     INTEGER MX(4,25)
      INTEGER N,NF,B,K
c     INTEGER ITEMPL(3,25)
      INTEGER ITEMPL(4,25)
c ??????????????????????????
      REAL*8 ELEM0(03)
C
c     INTEGER mstmx = 4000
c     INTEGER matlmx = 25
c     INTEGER hmx = 4
C
C     **** RUN THROUGH MSTATE FOR LEFT AND RIGHT ****
C     **** STATES TO EVALUATE HAM. BETWEEN;      ****
C     **** WILL SET UP STATES AS ACTING TO RIGHT ****
C     **** AND CONJUGATE THE LEFT STATE BEFORE   ****
C     **** INSERTING HAM OPS                     ****
C
C     **** MAX NUM OF B,B-ADJ,D,D-ADJ, AND ****
C     **** TOTAL OPS IN HAM                ****
      NHBMX=3
      NHBJMX=3
      NHDMX=3
      NHDJMX=3
      NOPMX=4
C
      DO 8 IH1=1,numsta
      DO 8 IH2=1,numsta
      do 9 ih3=1,nf
       HAM0(ih3,IH1,IH2)=0.0D0
    9 continue
       HAM(IH1,IH2)=0.0D0
       HNORM(IH1,IH2)=0.0D0
    8 CONTINUE
c
      DO 10 I1=1,4
      DO 10 I2=1,25
       ISTRT(I1,I2)=00
       ISTLT(I1,I2)=00
   10 CONTINUE
C
      DO 1000 NSTRT=1,NUMSTA
C
C      **** LOCATION OF RT STATE IN MSTATE ****
       LOCSTR=MSTINF(NSTRT,1)
C      **** LENGTH OF RIGHT STATE ****
       LSTRT=MSTINF(NSTRT,2)
C      **** NUM OF MES, BAR, ANTIBAR, AND ****
C      **** BAR+ANTIBAR FOR RIGHT STATE   ****
       NMESR=MSTINF(NSTRT,3)
       NBRBR=MSTINF(NSTRT,4)
       NBRDR=MSTINF(NSTRT,5)
       NBARR=NBRBR+NBRDR
C      **** NUM OF B'S,D'S IN STATE ****
       NBR=N*NBRBR+NMESR
       NDR=N*NBRDR+NMESR
C
c     WRITE (*,98)mstmx, matlmx
c  98 FORMAT(' mstmx, matlmx: ', 2I6) 
c
c     WRITE (*,97)LOCSTR, LSTRT
c  97 FORMAT(' LOCSTR, LSTRT: ', 2I6) 
c
c      if ((LOCSTR+3).GE.mstmx) then
c       WRITE(15,88)
c  88   FORMAT('LOCSTR exceeds mstmx in CLRDIS')
c       stop
c      else
c      endif
c
c      if ((LSTRT*2+hmx).GE.matlmx) then
c       WRITE(15,89)
c  89   FORMAT('state length exceeds max in CLRDIS')
c       stop
c      else
c      endif
C
C      NTERMS=N*(NBARR+NBARL)+2*(NMESR+NMESL)+NOPS
C
C
C      **** ASSIGN VALUES TO ISTRT FROM NSTRT'TH  ****
C      **** TERMS IN MSTATE (EXCEPT FOR COL.)     ****
c      DO 14 J1=1,3
       DO 14 J1=1,4
       DO 14 J2=1,LSTRT
        ISTRT(J1,J2)=MSTATE(LOCSTR+J1-1,J2)
   14  CONTINUE
C
C
      DO 1000 NSTLT=NSTRT,NUMSTA
C
       IF(INTPRT.NE.0) THEN
        WRITE(15,4)NSTLT,NSTRT
    4   FORMAT(' NSTLT,NSTRT: ',2I6)
       ELSE
       ENDIF
C
C
C      **** LOCATION OF LT STATE IN MSTATE ****
       LOCSTL=MSTINF(NSTLT,1)
C      **** LENGTH OF LEFT STATE ****
       LSTLT=MSTINF(NSTLT,2)
C      **** NUM OF MES, BAR, ANTIBAR, AND ****
C      **** BAR+ANTIBAR FOR LEFT STATES   ****
       NMESL=MSTINF(NSTLT,3)
       NBRBL=MSTINF(NSTLT,4)
       NBRDL=MSTINF(NSTLT,5)
       NBARL=NBRBL+NBRDL
C      **** NUM OF B'S,D'S IN STATE ****
       NBL=N*NBRBL+NMESL
       NDL=N*NBRDL+NMESL
C
C      NTERMS=N*(NBARR+NBARL)+2*(NMESR+NMESL)+NOPS
C
C      **** ASSIGN VALUES TO ISTLT FROM NSTLT'TH  ****
C      **** TERMS IN MSTATE (EXCEPT FOR COL.)     ****
C
c      DO 18 K1=1,3
       DO 18 K1=1,4
       DO 18 K2=1,LSTLT
        ISTLT(K1,K2)=MSTATE(LOCSTL+K1-1,K2)
C       **** CREATE TEMP COPY OF ISTLT FOR INITIAL ****
C       **** CHECK ON MAT EL BETWEEN LT,RT STATES  ****
        ITEMPL(K1,K2)=ISTLT(K1,K2)
   18  CONTINUE
c
C
C
C     **** PRINT ISTRT,ISTLT ****
C     WRITE(15,360)
C 360 FORMAT(' ISTRT:')
C     DO 365 L1=1,4
C      WRITE(15,370)(ISTRT(L1,L2),L2=1,LSTRT)
C 370  FORMAT(25I4)
C 365 CONTINUE
C     WRITE(15,362)
C 362 FORMAT(' ISTLT:')
C     DO 367 L2=1,4
C      WRITE(15,372)(ISTLT(L2,L3),L3=1,LSTLT)
C 372  FORMAT(25I4)
C 367 CONTINUE
C     WRITE(15,378)
C 378 FORMAT(' ------------------------- ')
C     *************************************
C
C
C     **** SEE HOW MANY B'S,D'S DIFFER IN RT,LT STATES ****
      NBSAME=0
      NDSAME=0
      DO 400 IR=1,LSTRT
      DO 410 IL=1,LSTLT
c ????????????????????????????????????
c      IF((ISTRT(1,IR).EQ.ITEMPL(1,IL)).AND.
c    >    (ISTRT(3,IR).EQ.ITEMPL(3,IL))) THEN
       IF((ISTRT(1,IR).EQ.ITEMPL(1,IL)).AND.
     >    (ISTRT(3,IR).EQ.ITEMPL(3,IL)).and.
     >    (ISTRT(4,IR).EQ.ITEMPL(4,IL))) THEN
        IF(ITEMPL(1,IL).EQ.1) THEN
         ITEMPL(1,IL)=2
         NBSAME=NBSAME+1
         GO TO 400
        ELSE
         ITEMPL(2,IL)=2
         NDSAME=NDSAME+1
         GO TO 400
        ENDIF
       ELSE
       ENDIF
  410 CONTINUE
  400 CONTINUE
C
C     **** NUM OF EXCESS B-ADJ,B,D-ADJ,D ****
C     **** WHEN COMPARE LT, RT STATES    ****
      NBAJX=ABS(NBR-NBSAME)
      NBX=ABS(NBL-NBSAME)
      NDAJX=ABS(NDR-NDSAME)
      NDX=ABS(NDL-NDSAME)
C     **** TOTAL EXCESS ****
      NBDX=NBAJX+NBX+NDAJX+NDX
C
C     WRITE(15,415)NBSAME,NDSAME,NBAJX,NBX,NDAJX,NDX
C 415 FORMAT(' NBSAME,NDSAME,NBAJX,NBX,NDAJX,NDX:',6I3)
C
C     **** CHECK THIS MAT EL ASSUMING THE ****
C     **** MAX NUM OF OPS IN HAM          ****
      IF((NHBMX.LT.NBAJX).OR.(NHBJMX.LT.NBX).OR.
     >   (NHDMX.LT.NDAJX).OR.(NHDJMX.LT.NDX).OR.
     >   (NOPMX.LT.NBDX)) THEN
C      **** CAN'T MATCH ANY STATES WITH HAM ****
       KILL=0
      ELSE
       KILL=1
      ENDIF
C
C     **** DIAGNOSTICS ****
c     WRITE(15,425)NBSAME,NDSAME
c 425 FORMAT(' NBSAME,NDSAME:',2I3)
c     WRITE(15,435)NBAJX,NBX,NDAJX,NDX,KILL
c 435 FORMAT(' NBAJX,NBX,NDAJX,NDX,KILL:',5I3)
C
C
      IF(KILL.EQ.0) GO TO 1000
C
C
C     **** PRINT MX ****
C     WRITE(15,160)
C 160 FORMAT(' MX:')
C     DO 165 L1=1,4
C      WRITE(15,170)(MX(L1,L2),L2=1,NTERMS)
C 170  FORMAT(25I4)
C 165 CONTINUE
C     WRITE(15,175)ISNPM
C 175 FORMAT(' ISNPM= ',I4)
C     WRITE(15,178)
C 178 FORMAT(' ------------------------- ')
C     *************************************
C
C
      CALL HAMQCD(LSTRT,LSTLT,ISTRT,ISTLT,ELEM0,ELEM,NORH)
C
cf    WRITE(15,163)(ELEM0(ifl),ifl=1,nf),ELEM
cf163 FORMAT(' ELEM0,ELEM AFT. HAMQCD',7e24.16)
C
C
      IF(NORH.EQ.0)THEN
       HNORM(NSTLT,NSTRT)=HNORM(NSTLT,NSTRT)+ELEM
      ELSE
       do 1100 ifl=1,nf
        HAM0(ifl,NSTLT,NSTRT)=HAM0(ifl,NSTLT,NSTRT)+ELEM0(ifl)
 1100  continue
      HAM(NSTLT,NSTRT)=HAM(NSTLT,NSTRT)+ELEM
      ENDIF
C
C
C
C     WRITE(15,850)
C 850 FORMAT(' ========================= ')
C     IF(NORH.EQ.1) THEN
C      WRITE(15,900)NSTLT,NSTRT,HAM(NSTLT,NSTRT)
C 900  FORMAT(' HAM(',I3,',',I3,') = ',e24.16)
C     ELSE
C      WRITE(15,905)NSTLT,NSTRT,HNORM(NSTLT,NSTRT)
C 905  FORMAT(' HNORM(',I3,',',I3,') = ',e24.16)
C     ENDIF
C     WRITE(15,850)
C
c ???????????????????????????????????
      IF(NORH.EQ.1) THEN
       HAM(NSTRT,NSTLT)=HAM(NSTLT,NSTRT)
c      HAM0(NSTRT,NSTLT)=HAM0(NSTLT,NSTRT)
       do 1200 ifl=1,nf
        HAM0(ifl,NSTRT,NSTLT)=HAM0(ifl,NSTLT,NSTRT)
 1200  continue
      ELSE
       HNORM(NSTRT,NSTLT)=HNORM(NSTLT,NSTRT)
      ENDIF
C
 1000 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE OPMTCH(NHB,NHBAJ,NHD,NHDAJ,NIX)
C     **** GIVEN THE NUM OF B,B-ADJ,ETC OPS IN  ****
C     **** HAM AND NUM OF DIFF OPS WHEN COMPARE ****
C     **** LT, RT STATES, CHECK IF MAT EL ZERO  ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/XBD/NBAJX,NBX,NDAJX,NDX
      COMMON/NBD/NBR,NDR,NBL,NDL
C
      NIX=1
      NIXPR=0
C
C     **** CHECK IF STATE HAS ENOUGH CREATION ****
C     **** OPS TO ABSORB HAM OPS              ****
      IF((NHB.GT.NBR).OR.(NHBAJ.GT.NBL).OR.
     >   (NHD.GT.NDR).OR.(NHDAJ.GT.NDL)) THEN
       NIX=0
       NIXPR=1
C      **** NIXPR IS FOR DIAGNOSTICS ****
C
      ELSE
C      **** CHECK THAT TOTAL NUM OF OPS MATCH ****
       IF(((NHB+NBL).NE.(NHBAJ+NBR)).OR.
     >    ((NHD+NDL).NE.(NHDAJ+NDR))) THEN
        NIX=0
        NIXPR=2
C
       ELSE
C       **** CHECK THAT HAM HAS ENOUGH OPS TO ****
C       **** MATCH THOSE UNMATCHED WHEN LEFT  ****
C       **** AND RIGHT STATES ARE COMPARED    ****
        IF((NHB.LT.NBAJX).OR.(NHBAJ.LT.NBX).OR.
     >     (NHD.LT.NDAJX).OR.(NHDAJ.LT.NDX)) THEN
         NIX=0
         NIXPR=3
C
        ELSE
        ENDIF
       ENDIF
      ENDIF
C
C     WRITE(15,10)NHB,NHBAJ,NHD,NHDAJ
C  10 FORMAT(' OPMTCH: NHB,NHBAJ,NHD,NHDAJ:',4I3)
C
C     WRITE(15,20)NBX,NBAJX,NDX,NDAJX,NIX,NIXPR
C  20 FORMAT(' OPMTCH: NBX,NBAJX,NDX,NDAJX,NIX,NIXPR:',6I2)
C
      RETURN
      END
C
C
C
      FUNCTION NFACT(M)
C     **** COMPUTE M FACTORIAL ****
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NFACT
      NFACT=1
      IF (M.NE.0) THEN
       DO 10 J1=1,M
        NFACT=NFACT*J1
   10  CONTINUE
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE PSIGN(LNGTH,IPERM,ISGN)
C     **** CALC SIGN OF ANTISYMM PERM ****
C     **** USED IN SUB EPS            ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IPERM(25)
C
      IF(LNGTH.LE.1)THEN
C      WRITE(15,10)LNGTH
C  10  FORMAT(' ERROR IN SUBR. PSIGN; LNGTH=',I3)
       ISGN=1
       RETURN
      ELSE
      ENDIF
C
      ISGN=1
C
C     **** INIT. SO THAT SGN(1,2,3,...,N)=+1 ****
      I0=0
      DO 30 L1=2,LNGTH
       DO 30 L2=1,L1-1
        I0=I0+1
   30 CONTINUE
      IF(((-1)**I0).LT.0) THEN
       IEXP=1
      ELSE
       IEXP=0
      ENDIF
C
      DO 20 I1=2,LNGTH
       DO 20 I2=1,I1-1
        IF(IPERM(I1).GT.IPERM(I2)) THEN
         IEXP=IEXP+1
        ELSE
        ENDIF
   20 CONTINUE
      ISGN=(-1)**(IEXP)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE HAMQCD(LRT0,LLT0,MATRT,MATLT,ELEM0,ELEM,NORH)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
c     common/masses/rmq(03)
      COMMON/NST/NUMSTA
      COMMON/SE/SELFEN(100)
      COMMON/PRTOPS/INTPRT
      COMMON/NOZERO/IB0,ID0
      COMMON/ME/MX,LLT,LRT,NOPS
      INTEGER N,NF,B,K
      INTEGER MATRT(4,25),MATLT(4,25),MATLTC(4,25)
      INTEGER MX(4,25)
      INTEGER MQNPOS(4,25)
      INTEGER NBDRL(4)
      INTEGER M(6)
c ?????????????????????
      REAL*8 ELEM0(03)
      REAL*8 H80(03),H90(03)
C
      LLT=LLT0
      LRT=LRT0
C
C     **** CONJUGATE LEFT STATE ****
      CALL CONJ(LLT,MATLT,MATLTC)
C
C     WRITE(15,10)
C  10 FORMAT(' MATLTC AFTER CONJ ')
C     DO 20 IJ=1,4
C      WRITE(15,25)(MATLTC(IJ,IK),IK=1,LLT)
C  25  FORMAT(25I3)
C  20 CONTINUE
C
C
C
      do 95 ifl0=1,nf
       ELEM0(ifl0)=0.0D0
   95 continue
      ELEM=0.0D0
C
      RLMBSQ=RLAMB*RLAMB
C
      DO 100 I6=1,LRT
      DO 100 J6=1,4
       MX(J6,I6)=MATRT(J6,I6)
  100 CONTINUE
C
C
C
C     **** SET UP LIST OF POSSIBLE QUANTUM NUMBERS TO     ****
C     **** ASSIGN TO HAM. OPERATORS.                      ****
C     **** MQNPOS: 1ST INDX IDS B-ADJ,D-ADJ,B,D;          ****
C     **** 2ND IDS WHICH B-ADJ (IE 1ST IN RT ST,2ND,ETC); ****
C     **** 3RD GIVES THE MOM (INDX=1) AND COLOR (INDX=2)  ****
C     **** OF THAT B-ADJ, D-ADJ ETC.                      ****
C
      DO 200 K6=1,4
       NBDRL(K6)=0
  200 CONTINUE
C
      DO 300 I9=1,LRT
      IF (MATRT(1,I9).EQ.1) THEN
        IF(NBDRL(1).NE.0) THEN
C       **** PREVENT DOUBLE COUNTING OF MOM ****
         IDC0=0
         DO 310 IDC=1,NBDRL(1)
          IF (MQNPOS(1,IDC).EQ.MATRT(3,I9)) THEN
           IDC0=1
          ELSE
          ENDIF
  310    CONTINUE
         IF(IDC0.NE.1) THEN
          NBDRL(1)=NBDRL(1)+1
          MQNPOS(1,NBDRL(1))=MATRT(3,I9)
         ELSE
         ENDIF
        ELSE
         NBDRL(1)=NBDRL(1)+1
         MQNPOS(1,NBDRL(1))=MATRT(3,I9)
        ENDIF
      ELSE
        IF(NBDRL(2).NE.0) THEN
C       **** PREVENT DOUBLE COUNTING OF MOM ****
         IDC0=0
         DO 312 IDC=1,NBDRL(2)
          IF (MQNPOS(2,IDC).EQ.MATRT(3,I9)) THEN
           IDC0=1
          ELSE
          ENDIF
  312    CONTINUE
         IF(IDC0.NE.1) THEN
          NBDRL(2)=NBDRL(2)+1
          MQNPOS(2,NBDRL(2))=MATRT(3,I9)
         ELSE
         ENDIF
        ELSE
         NBDRL(2)=NBDRL(2)+1
         MQNPOS(2,NBDRL(2))=MATRT(3,I9)
        ENDIF
      ENDIF
  300 CONTINUE
C
      DO 400 J9=1,LLT
      IF (MATLTC(1,J9).EQ.-1) THEN
        IF(NBDRL(3).NE.0) THEN
C       **** PREVENT DOUBLE COUNTING OF MOM ****
         IDC0=0
         DO 314 IDC=1,NBDRL(3)
          IF (MQNPOS(3,IDC).EQ.MATLTC(3,J9)) THEN
           IDC0=1
          ELSE
          ENDIF
  314    CONTINUE
         IF(IDC0.NE.1) THEN
          NBDRL(3)=NBDRL(3)+1
          MQNPOS(3,NBDRL(3))=MATLTC(3,J9)
         ELSE
         ENDIF
        ELSE
         NBDRL(3)=NBDRL(3)+1
         MQNPOS(3,NBDRL(3))=MATLTC(3,J9)
        ENDIF
      ELSE
        IF(NBDRL(4).NE.0) THEN
C       **** PREVENT DOUBLE COUNTING OF MOM ****
         IDC0=0
         DO 316 IDC=1,NBDRL(4)
          IF (MQNPOS(4,IDC).EQ.MATLTC(3,J9)) THEN
           IDC0=1
          ELSE
          ENDIF
  316    CONTINUE
         IF(IDC0.NE.1) THEN
          NBDRL(4)=NBDRL(4)+1
          MQNPOS(4,NBDRL(4))=MATLTC(3,J9)
         ELSE
         ENDIF
        ELSE
         NBDRL(4)=NBDRL(4)+1
         MQNPOS(4,NBDRL(4))=MATLTC(3,J9)
        ENDIF
      ENDIF
  400 CONTINUE
C
C     WRITE(15,470)
C 470 FORMAT(' NBDRL:')
C     WRITE(15,480)(NBDRL(KK),KK=1,4)
C 480 FORMAT(4I5)
C     WRITE(15,490)
C 490 FORMAT(' MQNPOS:')
C     DO 450 II=1,4
C     WRITE(15,460)(MQNPOS(II,JJ),JJ=1,NBDRL(II))
C 460 FORMAT(8I4)
C 450 CONTINUE
C
C
C
C
      IF(NORH.EQ.0) THEN
C     **** COMPUTE NORM MAT. FOR HILB. SP. ****
C     **** SET UP H0 ****
C
C     ISKP=0
C
      NOPS0=0
      NOPS = NOPS0
      NB0=0
      ND0=0
      NBAJ0=0
      NDAJ0=0
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX0=1
      CALL OPMTCH(NB0,NBAJ0,ND0,NDAJ0,NX0)
C     **** DIAGNOSTICS ****
C     WRITE(15,48)NX0
C  48 FORMAT(' NX0:',I3)
C
      IF(NX0.NE.0)THEN
C
      LNGTH=NOPS0+LRT+LLT
      DO 50 I0=1,LLT
      DO 50 J0=1,4
       MX(J0,LRT+NOPS0+I0)=MATLTC(J0,I0)
   50 CONTINUE
C
C     **** PRINT MX ****
c     WRITE(15,53)
c  53 FORMAT(' MX (IN SUB HAMQCD):')
c     DO 58 K0=1,4
c      WRITE(15,55)(MX(K0,M0),M0=1,LNGTH)
c  55  FORMAT(25I4)
c  58 CONTINUE
C
C
C
      H0=CLFACT(0,0,0,0)
C
      ELEM=ELEM+H0
C     WRITE(15,60)ELEM
C  60 FORMAT(' ELEM IN HAMQCD',e24.16)
C
C
      ELSE
      ENDIF
C
C
C
C
      ELSE
C     **** COMPUTE HAM MAT. RATHER THAN NORM MAT ****
C
C
C
C     **** FOUR POINT VERTICES ****
C     *****************************
C
C     **** SET UP H1 ****
C
      NOPS1=4
      NOPS = NOPS1
      NB1=2
      ND1=0
      NBAJ1=2
      NDAJ1=0
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX1=1
      CALL OPMTCH(NB1,NBAJ1,ND1,NDAJ1,NX1)
C     **** DIAGNOSTICS ****
C     WRITE(15,148)NX1
C 148 FORMAT(' NX1:',I3)
C
      IF(NX1.NE.0)THEN
C
      LNGTH=NOPS1+LRT+LLT
      DO 500 I1=1,LLT
      DO 500 J1=1,4
       MX(J1,LRT+NOPS1+I1)=MATLTC(J1,I1)
  500 CONTINUE
C
      MX(1,LRT+1)=-1
      MX(2,LRT+1)=00
      MX(1,LRT+2)=-1
      MX(2,LRT+2)=00
      MX(1,LRT+3)=+1
      MX(2,LRT+3)=00
      MX(1,LRT+4)=+1
      MX(2,LRT+4)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS1.NE.0) THEN
C
       DO 525 IM=1,NOPS1
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
  525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H1=0.D0
      DO 600 L1=1,NBDRL(M(1))
      DO 600 L2=1,NBDRL(M(2))
      DO 600 L3=1,NBDRL(M(3))
      DO 600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C
       DO 550 LL=1,NOPS1
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
  550  CONTINUE
C
C      WRITE(15,560)INTH
C 560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K4-K2
         KB=K3-K1
C
c ?????????????????????????????????????
c        **** run through flavor indices of operators   ****
c        **** four pt interactions will have two pairs  ****
c        **** of quarks with the same flavor (flavor is ****
c        **** conserved along fermion lines)            ****
         do 590 if1=1,nf
         do 590 jf1=1,nf
          mx(4,lrt+1)=if1
          mx(4,lrt+2)=jf1
          mx(4,lrt+3)=if1
          mx(4,lrt+4)=jf1
c
C        **** COL IN VERTEX ****
C
         H1=H1+(1.0D0/2.0D0)*
     >         (-CLFACT(1,4,2,3)*BRACK(KA,KB))
C
  590    continue
C
       ELSE
       ENDIF
C
  600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,700)H1
  700  FORMAT(' H1= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H1
C
      ELSE
      ENDIF
C
C
C
C
C     **** SET UP H2 ****
C
      NOPS2=4
      NOPS = NOPS2
      NB2=0
      ND2=2
      NBAJ2=0
      NDAJ2=2
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX2=1
      CALL OPMTCH(NB2,NBAJ2,ND2,NDAJ2,NX2)
C     **** DIAGNOSTICS ****
C     WRITE(15,1148)NX2
C1148 FORMAT(' NX2:',I3)
C
      IF(NX2.NE.0)THEN
C
      LNGTH=NOPS2+LRT+LLT
      DO 1500 I2=1,LLT
      DO 1500 J2=1,4
       MX(J2,LRT+NOPS2+I2)=MATLTC(J2,I2)
 1500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=-1
      MX(1,LRT+2)=00
      MX(2,LRT+2)=-1
      MX(1,LRT+3)=00
      MX(2,LRT+3)=+1
      MX(1,LRT+4)=00
      MX(2,LRT+4)=+1
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS2.NE.0) THEN
C
       DO 1525 IM=1,NOPS2
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 1525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H2=0.D0
      DO 1600 L1=1,NBDRL(M(1))
      DO 1600 L2=1,NBDRL(M(2))
      DO 1600 L3=1,NBDRL(M(3))
      DO 1600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C
       DO 1550 LL=1,NOPS2
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 1550  CONTINUE
C
C      WRITE(15,1560)INTH
C1560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K2-K4
         KB=K1-K3
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
         do 1590 if2=1,nf
         do 1590 jf2=1,nf
          mx(4,lrt+1)=if2
          mx(4,lrt+2)=jf2
          mx(4,lrt+3)=if2
          mx(4,lrt+4)=jf2
c
C        **** COL IN VERTEX ****
C
         H2=H2+(1.0D0/2.0D0)*
     >         (-CLFACT(3,2,4,1)*BRACK(KA,KB))
C
 1590    continue
C
       ELSE
       ENDIF
C
 1600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,1700)H2
 1700  FORMAT(' H2= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H2
C
      ELSE
      ENDIF
C
C
C
C
C     **** SET UP H3 ****
C
      NOPS3=4
      NOPS = NOPS3
      NB3=1
      ND3=0
      NBAJ3=2
      NDAJ3=1
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX3=1
      CALL OPMTCH(NB3,NBAJ3,ND3,NDAJ3,NX3)
C     **** DIAGNOSTICS ****
C     WRITE(15,2148)NX3
C2148 FORMAT(' NX3:',I3)
C
      IF(NX3.NE.0)THEN
C
      LNGTH=NOPS3+LRT+LLT
      DO 2500 I3=1,LLT
      DO 2500 J3=1,4
       MX(J3,LRT+NOPS3+I3)=MATLTC(J3,I3)
 2500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=+1
      MX(1,LRT+2)=-1
      MX(2,LRT+2)=00
      MX(1,LRT+3)=+1
      MX(2,LRT+3)=00
      MX(1,LRT+4)=+1
      MX(2,LRT+4)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS3.NE.0) THEN
C
       DO 2525 IM=1,NOPS3
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 2525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H3=0.D0
      DO 2600 L1=1,NBDRL(M(1))
      DO 2600 L2=1,NBDRL(M(2))
      DO 2600 L3=1,NBDRL(M(3))
      DO 2600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C      **** COLOR IN 1ST LOCATION: ****
C
       DO 2550 LL=1,NOPS3
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 2550  CONTINUE
C
C      WRITE(15,2560)INTH
C2560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K4-K2
         KB=K3+K1
         KC=K4+K1
         KD=K3-K2
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
         do 2590 if3=1,nf
         do 2590 jf3=1,nf
          mx(4,lrt+1)=if3
          mx(4,lrt+2)=jf3
          mx(4,lrt+3)=jf3
          mx(4,lrt+4)=if3
c
C        **** COL IN VERTEX ****
C
C        H3=H3+(1.0D0/2.0D0)*
C    >         (-CLFACT(1,4,2,3)*BRACK(KA,KB)+
C    >           CLFACT(2,4,1,3)*BRACK(KC,KD))
C
         H3=H3+( CLFACT(2,4,1,3)*BRACK(KC,KD))
C
 2590    continue
C
C
       ELSE
       ENDIF
C
 2600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,2700)H3
 2700  FORMAT(' H3= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H3
C
      ELSE
      ENDIF
C
C
C
C     **** SET UP H4 ****
C
      NOPS4=4
      NOPS = NOPS4
      NB4=0
      ND4=1
      NBAJ4=1
      NDAJ4=2
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX4=1
      CALL OPMTCH(NB4,NBAJ4,ND4,NDAJ4,NX4)
C     **** DIAGNOSTICS ****
C     WRITE(15,3148)NX4
C3148 FORMAT(' NX4:',I3)
C
      IF(NX4.NE.0)THEN
C
      LNGTH=NOPS4+LRT+LLT
      DO 3500 I4=1,LLT
      DO 3500 J4=1,4
       MX(J4,LRT+NOPS4+I4)=MATLTC(J4,I4)
 3500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=-1
      MX(1,LRT+2)=00
      MX(2,LRT+2)=+1
      MX(1,LRT+3)=00
      MX(2,LRT+3)=+1
      MX(1,LRT+4)=+1
      MX(2,LRT+4)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS4.NE.0) THEN
C
       DO 3525 IM=1,NOPS4
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 3525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H4=0.D0
      DO 3600 L1=1,NBDRL(M(1))
      DO 3600 L2=1,NBDRL(M(2))
      DO 3600 L3=1,NBDRL(M(3))
      DO 3600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C
       DO 3550 LL=1,NOPS4
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 3550  CONTINUE
C
C      WRITE(15,3560)INTH
C3560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K1-K3
         KB=-K4-K2
         KC=K4+K3
         KD=K2-K1
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
         do 3590 if4=1,nf
         do 3590 jf4=1,nf
          mx(4,lrt+1)=if4
          mx(4,lrt+2)=jf4
          mx(4,lrt+3)=if4
          mx(4,lrt+4)=jf4
c
C        **** COL IN VERTEX ****
C
C        H4=H4+(1.0D0/2.0D0)*
C    >         ( CLFACT(2,1,3,4)*BRACK(KA,KB)-
C    >           CLFACT(2,4,3,1)*BRACK(KC,KD))
C
         H4=H4+( CLFACT(2,1,3,4)*BRACK(KA,KB))
C
 3590    continue
C
C
       ELSE
       ENDIF
C
 3600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,3700)H4
 3700  FORMAT(' H4= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H4
C
      ELSE
      ENDIF
C
C
C
C     **** SET UP H5 ****
C
      NOPS5=4
      NOPS = NOPS5
      NB5=2
      ND5=1
      NBAJ5=1
      NDAJ5=0
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX5=1
      CALL OPMTCH(NB5,NBAJ5,ND5,NDAJ5,NX5)
C     **** DIAGNOSTICS ****
C     WRITE(15,4148)NX5
C4148 FORMAT(' NX5:',I3)
C
      IF(NX5.NE.0)THEN
C
      LNGTH=NOPS5+LRT+LLT
      DO 4500 I5=1,LLT
      DO 4500 J5=1,4
       MX(J5,LRT+NOPS5+I5)=MATLTC(J5,I5)
 4500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=-1
      MX(1,LRT+2)=-1
      MX(2,LRT+2)=00
      MX(1,LRT+3)=-1
      MX(2,LRT+3)=00
      MX(1,LRT+4)=+1
      MX(2,LRT+4)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS5.NE.0) THEN
C
       DO 4525 IM=1,NOPS5
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 4525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H5=0.D0
      DO 4600 L1=1,NBDRL(M(1))
      DO 4600 L2=1,NBDRL(M(2))
      DO 4600 L3=1,NBDRL(M(3))
      DO 4600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C
       DO 4550 LL=1,NOPS5
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 4550  CONTINUE
C
C      WRITE(15,4560)INTH
C4560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K4-K3
         KB=-K1-K2
         KC=K1+K3
         KD=K2-K4
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
         do 4590 if5=1,nf
         do 4590 jf5=1,nf
          mx(4,lrt+1)=if5
          mx(4,lrt+2)=jf5
          mx(4,lrt+3)=if5
          mx(4,lrt+4)=jf5
c
C        **** COL IN VERTEX ****
C
C        H5=H5+(1.0D0/2.0D0)*
C    >         (-CLFACT(2,4,3,1)*BRACK(KA,KB)+
C    >           CLFACT(2,1,3,4)*BRACK(KC,KD))
C
         H5=H5+( CLFACT(2,1,3,4)*BRACK(KC,KD))
C
 4590    continue
C
C
       ELSE
       ENDIF
C
 4600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,4700)H5
 4700  FORMAT(' H5= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H5
C
      ELSE
      ENDIF
C
C
C
C     **** SET UP H6 ****
C
      NOPS6=4
      NOPS = NOPS6
      NB6=1
      ND6=2
      NBAJ6=0
      NDAJ6=1
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX6=1
      CALL OPMTCH(NB6,NBAJ6,ND6,NDAJ6,NX6)
C     **** DIAGNOSTICS ****
C     WRITE(15,5148)NX6
C5148 FORMAT(' NX6:',I3)
C
      IF(NX6.NE.0)THEN
C
      LNGTH=NOPS6+LRT+LLT
      DO 5500 I6=1,LLT
      DO 5500 J6=1,4
       MX(J6,LRT+NOPS6+I6)=MATLTC(J6,I6)
 5500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=-1
      MX(1,LRT+2)=00
      MX(2,LRT+2)=-1
      MX(1,LRT+3)=00
      MX(2,LRT+3)=+1
      MX(1,LRT+4)=-1
      MX(2,LRT+4)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS6.NE.0) THEN
C
       DO 5525 IM=1,NOPS6
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 5525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H6=0.D0
      DO 5600 L1=1,NBDRL(M(1))
      DO 5600 L2=1,NBDRL(M(2))
      DO 5600 L3=1,NBDRL(M(3))
      DO 5600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C
       DO 5550 LL=1,NOPS6
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 5550  CONTINUE
C
C      WRITE(15,5560)INTH
C5560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K2-K3
         KB=K1+K4
         KC=K2+K4
         KD=K1-K3
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
         do 5590 if6=1,nf
         do 5590 jf6=1,nf
          mx(4,lrt+1)=if6
          mx(4,lrt+2)=jf6
          mx(4,lrt+3)=jf6
          mx(4,lrt+4)=if6
c
C        **** COL IN VERTEX ****
C
C        H6=H6+(1.0D0/2.0D0)*
C    >         ( CLFACT(4,2,3,1)*BRACK(KA,KB)-
C    >           CLFACT(3,2,4,1)*BRACK(KC,KD))
C
         H6=H6+( CLFACT(4,2,3,1)*BRACK(KA,KB))
C
 5590    continue
C
C
       ELSE
       ENDIF
C
 5600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,5700)H6
 5700  FORMAT(' H6= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H6
C
      ELSE
      ENDIF
C
C
C
C     **** SET UP H7 ****
C
      NOPS7=4
      NOPS = NOPS7
      NB7=1
      ND7=1
      NBAJ7=1
      NDAJ7=1
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX7=1
      CALL OPMTCH(NB7,NBAJ7,ND7,NDAJ7,NX7)
C     **** DIAGNOSTICS ****
C     WRITE(15,6148)NX7
C6148 FORMAT(' NX7:',I3)
C
      IF(NX7.NE.0)THEN
C
      LNGTH=NOPS7+LRT+LLT
      DO 6500 I7=1,LLT
      DO 6500 J7=1,4
       MX(J7,LRT+NOPS7+I7)=MATLTC(J7,I7)
 6500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=-1
      MX(1,LRT+2)=00
      MX(2,LRT+2)=+1
      MX(1,LRT+3)=-1
      MX(2,LRT+3)=00
      MX(1,LRT+4)=+1
      MX(2,LRT+4)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS7.NE.0) THEN
C
       DO 6525 IM=1,NOPS7
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 6525  CONTINUE
C
      ELSE
      ENDIF
C
C
      H7=0.D0
      DO 6600 L1=1,NBDRL(M(1))
      DO 6600 L2=1,NBDRL(M(2))
      DO 6600 L3=1,NBDRL(M(3))
      DO 6600 L4=1,NBDRL(M(4))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
       MX(3,LRT+3)=MQNPOS(M(3),L3)
       MX(3,LRT+4)=MQNPOS(M(4),L4)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C      **** TO SATISFY DELTA'S IN HAM;     ****
C      **** ASSIGN B A LOWER COL INDEX,AND ****
C      **** POS. NUM FOR UPPER, NEG FOR    ****
C      **** LOWER INDEX WHEN COMPARE BELOW ****
C
       INTH=0
C
       DO 6550 LL=1,NOPS7
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 6550  CONTINUE
C
C      WRITE(15,6560)INTH
C6560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
         K3=MX(3,LRT+3)
         K4=MX(3,LRT+4)
         KA=K4-K3
         KB=K2-K1
         KC=K4+K2
         KD=-K1-K3
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
c        **** need to split these two vertices since they ****
c        **** differ in flavor flow                       ****
         do 6590 if7a=1,nf
         do 6590 jf7a=1,nf
          mx(4,lrt+1)=if7a
          mx(4,lrt+2)=if7a
          mx(4,lrt+3)=jf7a
          mx(4,lrt+4)=jf7a
c
C        **** COL IN VERTEX ****
C
c        H7=H7+(1.0D0/2.0D0)*
c    >     2.0D0*(-CLFACT(2,4,3,1)*BRACK(KA,KB)+
c    >           CLFACT(3,4,2,1)*BRACK(KC,KD))
c
         H7=H7+(1.0D0/2.0D0)*
     >     2.0D0*(-CLFACT(2,4,3,1)*BRACK(KA,KB) )
C
 6590    continue
C
c        **** run through flavor indices ****
         do 6595 if7b=1,nf
         do 6595 jf7b=1,nf
          mx(4,lrt+1)=if7b
          mx(4,lrt+2)=jf7b
          mx(4,lrt+3)=if7b
          mx(4,lrt+4)=jf7b
c
C        **** COL IN VERTEX ****
C
         H7=H7+(1.0D0/2.0D0)*
     >     2.0d0*(CLFACT(3,4,2,1)*BRACK(KC,KD))
c
 6595    continue
C
C
       ELSE
       ENDIF
C
 6600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,6700)H7
 6700  FORMAT(' H7= ',e24.16)
      ELSE
      ENDIF
C
      ELEM=ELEM+H7
C
      ELSE
      ENDIF
C
C
C
C     **** DIAGONAL VERTICES ****
C     ***************************
C
      cbreak=1.d-8
c     **** charge conjugation breaking parameter;  ****
c     **** used to split states slightly to make   ****
c     **** degenerate wavefunctions better behaved ****
c     **** will give particles a slightly higher   ****
c     **** mass than antiparticles                 ****
c
C     **** SET UP H8 ****
C
      NOPS8=2
      NOPS = NOPS8
      NB8=1
      ND8=0
      NBAJ8=1
      NDAJ8=0
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX8=1
      CALL OPMTCH(NB8,NBAJ8,ND8,NDAJ8,NX8)
C     **** DIAGNOSTICS ****
C     WRITE(15,7148)NX8
C7148 FORMAT(' NX8:',I3)
C
      IF(NX8.NE.0)THEN
C
      LNGTH=NOPS8+LRT+LLT
      DO 7500 I8=1,LLT
      DO 7500 J8=1,4
       MX(J8,LRT+NOPS8+I8)=MATLTC(J8,I8)
 7500 CONTINUE
C
      MX(1,LRT+1)=-1
      MX(2,LRT+1)=00
      MX(1,LRT+2)=+1
      MX(2,LRT+2)=00
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS8.NE.0) THEN
C
       DO 7525 IM=1,NOPS8
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 7525  CONTINUE
C
      ELSE
      ENDIF
C
C
      do 7530 iflh8=1,nf
       H80(iflh8)=0.D0
 7530 continue
      H8=0.D0
      DO 7600 L1=1,NBDRL(M(1))
      DO 7600 L2=1,NBDRL(M(2))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C
       INTH=0
C
       DO 7550 LL=1,NOPS8
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 7550  CONTINUE
C
C      WRITE(15,7560)INTH
C7560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
C
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
c        **** diagonal vertices have one flavor index ****
         do 7590 ifl8=1,nf
          mx(4,lrt+1)=ifl8
          mx(4,lrt+2)=ifl8
c
C
         IF(K1.NE.0) THEN
          IF(MASS.EQ.0) THEN
           H8=H8+CLFACT(1,2,0,0)*DFLOAT(JDELTA(K1,K2))*
     >               ( (1.0D0/2.0D0)*SELFEN(K1) )
          ELSE
           H80(ifl8)=H80(ifl8)+CLFACT(1,2,0,0)*DFLOAT(JDELTA(K1,K2))*
     >              ( (2.0D0 + cbreak)/K1 )
           H8=H8+CLFACT(1,2,0,0)*DFLOAT(JDELTA(K1,K2))*
     >              ( (1.0D0/2.0D0)*SELFEN(K1) )
          ENDIF
         ELSE
C        ****  NO CONTRIB. IF K1=0  ****
         ENDIF
C
 7590    continue
C
C
       ELSE
       ENDIF
C
 7600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,7700)(H80(ifl),ifl=1,nf),H8
 7700  FORMAT(' H80(nf),H8= ',7e24.16)
      ELSE
      ENDIF
C
      do 7750 ifle8=1,nf
       ELEM0(ifle8)=ELEM0(ifle8)+H80(ifle8)
 7750 continue
      ELEM=ELEM+H8
C
      ELSE
      ENDIF
C
C
C
C     **** SET UP H9 ****
C
      NOPS9=2
      NOPS = NOPS9
      NB9=0
      ND9=1
      NBAJ9=0
      NDAJ9=1
C
C     **** SEE IF CAN MATCH OPS IN HAM TO ****
C     **** THOSE LEFT OVER IN STATES      ****
      NX9=1
      CALL OPMTCH(NB9,NBAJ9,ND9,NDAJ9,NX9)
C     **** DIAGNOSTICS ****
C     WRITE(15,8148)NX9
C8148 FORMAT(' NX9:',I3)
C
      IF(NX9.NE.0)THEN
C
      LNGTH=NOPS9+LRT+LLT
      DO 8500 I9=1,LLT
      DO 8500 J9=1,4
       MX(J9,LRT+NOPS9+I9)=MATLTC(J9,I9)
 8500 CONTINUE
C
      MX(1,LRT+1)=00
      MX(2,LRT+1)=-1
      MX(1,LRT+2)=00
      MX(2,LRT+2)=+1
C
C
C     **** M(I) CORRESPONDS TO OP NUM I IN HAM;          ****
C     **** IF THE I'TH OP IS A B-ADJ,D-ADJ,B,OR D, THEN  ****
C     **** ASSIGN M(I) A VAL OF 1,2,3,OR 4 TO INDICATE   ****
C     **** WHERE IN NBDRL, MQNPOS TO LOOK FOR APPROPR    ****
C     **** QM NUMS FOR THAT OP                           ****
C
      IF(NOPS9.NE.0) THEN
C
       DO 8525 IM=1,NOPS9
        MU=MX(1,LRT+IM)
        MD=MX(2,LRT+IM)
        M(IM)=2*MU*MU+3*MD*MD+MU+MD
 8525  CONTINUE
C
      ELSE
      ENDIF
C
C
      do 8530 iflh9=1,nf
       H90(iflh9)=0.D0
 8530 continue
      H9=0.D0
      DO 8600 L1=1,NBDRL(M(1))
      DO 8600 L2=1,NBDRL(M(2))
       MX(3,LRT+1)=MQNPOS(M(1),L1)
       MX(3,LRT+2)=MQNPOS(M(2),L2)
C
C
C
C      **** CHECK MOMENTUM CONSERVATION    ****
C
       INTH=0
C
       DO 8550 LL=1,NOPS9
C       **** SUM MOM. INTO VERTEX ****
        INTH=INTH+MX(3,LRT+LL)*(MX(1,LRT+LL)+MX(2,LRT+LL))
 8550  CONTINUE
C
C      WRITE(15,8560)INTH
C8560  FORMAT(' MOM MATCH; INTH:',I3)
C
       IF (INTH.EQ.0) THEN
C
C
C
C        **** MOM IN VERTEX ****
         K1=MX(3,LRT+1)
         K2=MX(3,LRT+2)
C
C
c ?????????????????????????????????????
c        **** run through flavor indices ****
         do 8590 ifl9=1,nf
          mx(4,lrt+1)=ifl9
          mx(4,lrt+2)=ifl9
c
C
         IF(K1.NE.0) THEN
          IF(MASS.EQ.0) THEN
           H9=H9+CLFACT(2,1,0,0)*DFLOAT(JDELTA(K1,K2))*
     >               ( (1.0D0/2.0D0)*SELFEN(K1) )
          ELSE
           H90(ifl9)=H90(ifl9)+CLFACT(2,1,0,0)*DFLOAT(JDELTA(K1,K2))*
     >              ( (2.0D0 - cbreak)/K1 )
           H9=H9+CLFACT(2,1,0,0)*DFLOAT(JDELTA(K1,K2))*
     >              ( (1.0D0/2.0D0)*SELFEN(K1) )
          ENDIF
         ELSE
C        ****  NO CONTRIB. IF K1=0  ****
         ENDIF
C
 8590    continue
C
C
       ELSE
       ENDIF
C
 8600 CONTINUE
C
      IF(INTPRT.NE.0) THEN
       WRITE(15,8700)(H90(ifl),ifl=1,nf),H9
 8700  FORMAT(' H90(nf),H9= ',7e24.16)
      ELSE
      ENDIF
C
      do 8750 ifle9=1,nf
       ELEM0(ifle9)=ELEM0(ifle9)+H90(ifle9)
 8750 continue
      ELEM=ELEM+H9
C
      ELSE
      ENDIF
C
c
C
      ENDIF
c      **** computing ham rather than hnorm ****
C
C
      RETURN
      END
C
C
C
C
      FUNCTION JDELTA(I,J)
C     **** KRON. DELTA ****
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(I.EQ.J) THEN
       JDELTA=1
      ELSE
       JDELTA=0
      ENDIF
      RETURN
      END
C
C
C
C
C
C
C
      SUBROUTINE CONJ(LST,IST,ISTC)
C     **** CONJUGATE STATE IN IST OF LENGTH LST ****
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IST(4,25),ISTC(4,25)
C
C     WRITE(15,50)
C  50 FORMAT(' IST IN SUB CONJ BEFORE ')
C     DO 60 IJ=1,4
C      WRITE(15,65)(IST(IJ,IK),IK=1,LST)
C  65  FORMAT(25I3)
C  60 CONTINUE
C
      DO 10 I1=1,LST
      DO 10 I2=1,4
       J1=LST+1-I1
       IF(I2.LE.2) THEN
        ISTC(I2,J1)=-IST(I2,I1)
       ELSE
        ISTC(I2,J1)=IST(I2,I1)
       ENDIF
   10 CONTINUE
C
C
C     WRITE(15,80)
C  80 FORMAT(' ISTC IN SUB CONJ AFTER ')
C     DO 90 JJ=1,4
C      WRITE(15,95)(ISTC(JJ,JK),JK=1,LST)
C  95  FORMAT(25I3)
C  90 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE PRTMAT(IDIM,MAT,LLHLF)
C     **** PRINT REAL MAT (IF LLHLF=0 ****
C     **** PRINT ONLY LOWER LEFT HALF)****
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MAT(6902,6902)
C
C     **** PRINT NORM MATRIX ****
      IF(IDIM.NE.0) THEN
C      WRITE(15,750)
C 750  FORMAT(' MATRIX: ')
       DO 775 IM=1,IDIM,6
       IF(LLHLF.EQ.0) THEN
        IR1=IM
       ELSE
        IR1=1
       ENDIF
       DO 700 IR=IR1,IDIM
        LL=IM+5
        IF(LLHLF.EQ.0) THEN
          IF(LL.GT.IR) LL=IR
        ELSE
          IF(LL.GT.IDIM) LL=IDIM
        ENDIF
        WRITE(15,800)(MAT(IR,IC),IC=IM,LL)
  800   FORMAT(6e24.16)
  700  CONTINUE
       WRITE(15,760)
  760  FORMAT(' -------')
  775  CONTINUE
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE PRTSPM(IDIM,MAT,LLHLF)
C     **** PRINT REAL MAT (IF LLHLF=0 ****
C     **** PRINT ONLY LOWER LEFT HALF)****
C     **** THIS IS FOR SPARSE MATS;   ****
C     **** ONLY PRINT NON-ZERO ELEMS  ****
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MAT(6902,6902), MATNZ(10)
      INTEGER IJ(20)
C
C
      NROW=5
C     **** NUMBER OF TERMS IN EACH ROW OF OUTPUT ****
C
      IF(IDIM.NE.0) THEN
       NCOUNT=0
       EPS=1.0D-10
C      **** CUTOFF TO DETERMINE A ZERO ELEMENT ****
       DO 10 IR=1,IDIM
        IF(LLHLF.EQ.0) THEN
         IL1=IR
        ELSE
         IL1=1
        ENDIF
        DO 20 IL=IL1,IDIM
         IF(ABS(MAT(IL,IR)).GT.EPS) THEN
          NCOUNT=NCOUNT+1
          MATNZ(NCOUNT)=MAT(IL,IR)
          IJ(2*NCOUNT-1)=IL
          IJ(2*NCOUNT)=IR
          IF(NCOUNT.EQ.NROW) THEN
           WRITE(15,30)(IJ(L),L=1,2*NROW)
   30      FORMAT(4X,5(' (',I6,',',I6,')',3X))
C          **** NEED " ABOVE TO EQUAL NROW  ****
           WRITE(15,40)(MATNZ(L),L=1,NROW)
   40      FORMAT(10(e24.16))
           NCOUNT=0
           DO 50 IZ=1,NROW
            MATNZ(IZ)=0.0D0
            IJ(2*IZ-1)=0
            IJ(2*IZ)=0
   50      CONTINUE
          ELSE
          ENDIF
         ELSE
         ENDIF
   20   CONTINUE
   10  CONTINUE
C
       IF(NCOUNT.NE.0) THEN
C      **** PRINT LAST ROW ****
        WRITE(15,60)(IJ(L),L=1,2*NCOUNT)
   60   FORMAT(4X,5(' (',I6,',',I6,')',3X))
        WRITE(15,70)(MATNZ(L),L=1,NCOUNT)
   70   FORMAT(10(e24.16))
       ELSE
       ENDIF
C
      ELSE
      ENDIF
C
      WRITE(15,80)
   80 FORMAT(//)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE WEEDR
C     **** SEARCH MAT OF INNER PROD OF STATES     ****
C     **** TO WEED OUT THOSE WHICH ARE REDUNDANT; ****
C     **** REVISE BOTH THIS MAT AND MSTINF, WHICH ****
C     **** HOLDS INFO ON STATES,INCL LOCATION     ****
C     **** IN MSTATE. NOTE THAT IT IS NOT NECESS. ****
C     **** TO DROP THESE REDUND. STATES FROM      ****
C     **** MSTATE, AS STATES ARE ALWAYS ACCESSED  ****
C     **** VIA MSTINF                             ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NM/HNORM(6902,6902)
      INTEGER MDROP(6902)
C
      DO 10 II=1,6902
       MDROP(II)=0
   10 CONTINUE
C
      EPS=1.0D-4
C     **** NEED EPS NE 0 FOR COMPARING WITH SMALL NUMS ****
      NSTPR=NUMSTA
C     **** NSTPR HOLDS CURRENT NUMBER OF STATES LEFT ****
C     LOCPR=1
C     IMDR=1
C
      LOCFIN=NSTPR-1
      DO 20 LOC=1,LOCFIN
C     WRITE(15,22)LOC
C  22 FORMAT(' LOC: ',I4)
      IF(LOC.LE.(NSTPR-1)) THEN
C
       NDUMP=0
       JFIN=NSTPR
       DO 30 J=LOC+1,JFIN
C      WRITE(15,23)J
C  23  FORMAT(' J  : ',I4)
       IF(J.LE.NSTPR) THEN
C
        IF(DABS(HNORM(J,LOC)).GT.EPS) THEN
         IDUMP=1
         R=HNORM(J,LOC)/HNORM(LOC,LOC)
C
C        WRITE(15,61)HNORM(J,LOC),HNORM(LOC,LOC)
C  61    FORMAT(' N(J,LOC),N(LOC,LOC):',2F8.3)
C
         DO 60 IR=1,NSTPR
          RH=R*HNORM(LOC,IR)
          DIFF=RH-HNORM(J,IR)
C
C         WRITE(15,62)J,IR,LOC,IR,HNORM(J,IR),HNORM(LOC,IR),R,RH,DIFF
C  62     FORMAT(' H(',I3,',',I3,'),H(',I3,',',I3,
C    >             '),R,RH,DIFF:',5F7.2)
C
          IF(DABS(DIFF).GT.EPS) THEN
           IDUMP=0
          ELSE
          ENDIF
   60    CONTINUE
C
C        WRITE(15,64)IDUMP
C  64    FORMAT(' IDUMP:',I3)
C
         IF(IDUMP.EQ.1) THEN
          NDUMP=NDUMP+1
          MDROP(NDUMP)=J
C         WRITE(15,66)NDUMP, J
C  66     FORMAT(' NDUMP,J : ',2I4)
C         **** MDROP STORES THE CURRENT NUMS ****
C         **** OF THE STATES TO BE DROPPED   ****
         ELSE
         ENDIF
C
        ELSE
        ENDIF
C
       ELSE
       ENDIF
   30  CONTINUE
C
C      WRITE(15,32)LOC,NDUMP
C  32  FORMAT(' LOC,NDUMP: ',2I3)
C      IF(NDUMP.GT.0) THEN
C       WRITE(15,34)
C  34   FORMAT(' MDROP:')
C       WRITE(15,36)(MDROP(I0),I0=1,NDUMP)
C  36   FORMAT(10I3)
C      ELSE
C      ENDIF
C
C      **** REMOVE EXCESS STATE ****
       IF(NDUMP.GT.0) THEN
        CALL DROPR(NSTPR,NDUMP,MDROP)
       ELSE
       ENDIF
C
      ELSE
      ENDIF
   20 CONTINUE
C
C
C     **** UPDATE THE NUMBER OF STATES ****
      NUMSTA=NSTPR
C
C
      RETURN
      END
C
C
C
C
      SUBROUTINE WEEDR2(W,Z,IDR)
C     **** SEARCH FOR ZERO EVALS AND DISCARD EXCESS STATES ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      REAL*8 W(6902),Z(6902,6902)
      INTEGER MDROP(6902),INDCHK(6902)
C
      EPS=1.0D-4
      IDR=0
      NZER=0
      DO 70 M1=1,NUMSTA
       MDROP(M1)=0
       INDCHK(M1)=0
   70 CONTINUE
C
C     **** COUNT ZERO EIGENVALUES ****
      DO 10 L1=1,NUMSTA
       IF(DABS(W(L1)).LT.EPS) THEN
        NZER=NZER+1
       ELSE
       ENDIF
   10 CONTINUE
C
      IF(NZER.GT.0) THEN
C
       DO 20 I1=1,NZER
        IDIS=0
        DO 30 I2=NUMSTA,1,-1
C
         IF((DABS(Z(I2,I1)).GT.EPS).AND.(INDCHK(I2).EQ.0)) THEN
C
C        **** USE Z TO FIND A STATE WHICH IS PART OF THE  ****
C        **** ZERO E-VECT TO DISCARD.  ONCE A STATE IS    ****
C        **** DISCARDED, STORE THE OTHER STATES IN THAT   ****
C        **** E-VECT IN INDCHK AND DON'T ALLOW ANY MORE   ****
C        **** TO BE THROWN OUT. THIS ENSURES THAT THE NZER****
C        **** STATES EVENTUALLY DISCARDED ARE INDEPENDENT ****
C
          IF(IDIS.EQ.0) THEN
           IDR=IDR+1
           MDROP(IDR)=I2
           IDIS=1
          ELSE
          ENDIF
          INDCHK(I2)=1
C
         ELSE
         ENDIF
C
   30   CONTINUE
   20  CONTINUE
C
      ELSE
      ENDIF
C
C
   50 IF(IDR.NE.0) THEN
       WRITE(15,60)IDR
   60  FORMAT('HAVE DROPPED ',I5,' ADD. STATES AFTER DIAG')
       WRITE(15,80)
   80  FORMAT(' MDROP:')
       WRITE(15,90)(MDROP(K1),K1=1,IDR)
   90  FORMAT(10I5)
       CALL DROPR(NUMSTA,IDR,MDROP)
       WRITE(15,94)NUMSTA
   94  FORMAT(' NUMSTA= ',I6)
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE DROPR(NSTPR,NDUMP,MDROP)
C     **** GIVEN THE NUMBER OF STATES TO BE DROPPED IN ****
C     **** NDUMP, AND THEIR STATE NUM IN MDROP,REMOVE  ****
C     **** THEM FROM MSTATE, MSTINF,AND HNORM,AND GIVE ****
C     **** THE NEW NUMBER OF STATES IN NSTPR           ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
C     COMMON/NST/NUMSTA
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NM/HNORM(6902,6902)
      INTEGER MDROP(6902)
C
C     WRITE(15,5)
C   5 FORMAT(' SUB. DROPR:')
C     WRITE(15,10)NDUMP
C  10 FORMAT(' NDUMP: ',I3)
C     WRITE(15,45)NSTPR
C  45 FORMAT(' NSTPR: ',I3)
C
      IF(NDUMP.NE.0) THEN
C
       DO 70 ID=1,NDUMP
        JDRI=MDROP(ID)
        JDRF=NSTPR
C
        IF(JDRI.EQ.NSTPR) THEN
         DO 40 J2=1,NSTPR
          HNORM(NSTPR,J2)=0.0D0
          HNORM(J2,NSTPR)=0.0D0
   40    CONTINUE
C
         DO 50 K2=1,8
          MSTINF(NSTPR,K2)=0
   50    CONTINUE
C
         NSTPR=NSTPR-1
C
C        WRITE(15,45)NSTPR
C
        ELSE
C
         DO 80 JDR=JDRI,JDRF
         DO 80 KR=1,NSTPR
          HNORM(JDR,KR)=HNORM(JDR+1,KR)
          HNORM(JDR+1,KR)=0.0D0
   80    CONTINUE
C
         DO 90 LDR=JDRI,JDRF
         DO 90 MR=1,NSTPR
          HNORM(MR,LDR)=HNORM(MR,LDR+1)
          HNORM(MR,LDR+1)=0.0D0
   90    CONTINUE
C
         DO 100 NS=JDRI,JDRF
C
          DO 110 KS=1,8
           MSTINF(NS,KS)=MSTINF(NS+1,KS)
           MSTINF(NS+1,KS)=0
  110     CONTINUE
C
C
  100    CONTINUE
C
         NSTPR=NSTPR-1
C
C        WRITE(15,45)NSTPR
C
         IF(ID.NE.NDUMP) THEN
C        **** UPDATE LOC OF STATES TO BE AXED ****
         DO 130 IE=ID+1,NDUMP
          IF(MDROP(IE).GT.JDRI) THEN
           MDROP(IE)=MDROP(IE)-1
          ELSE
          ENDIF
  130    CONTINUE
         ELSE
         ENDIF
C
        ENDIF
   70  CONTINUE
C
      ELSE
      ENDIF
C
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE DIAG(RSMAT,W,Z,IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
c     REAL*8 RSMAT(6902,6902),W(6902),Z(6902,6902),FV1(6902),FV2(6902)
      REAL*8 RSMAT(6902,6902),W(6902),Z(6902,6902),e(6902)
c     INTEGER IERR,NM,N1,MATZ
      INTEGER ierr,NM,N1
C
C     WRITE(15,10)
C  10 FORMAT(' RSMAT IN SUB DIAG:')
C     CALL PRTMAT(NUMSTA,RSMAT,1)
C
C
      N1 = NUMSTA
      NM = 6902
c     MATZ = 1
      ierr=0
c     **** ierr not used in num recipes ****
c
c     CALL RS(NM,N1,RSMAT,W,MATZ,Z,FV1,FV2,IERR)
c
      do 20 iz=1,nm
      do 25 jz=1,nm
        z(iz,jz)=rsmat(iz,jz)
   25 continue
   20 continue
c
      call trr8(z,n1,nm,w,e)
c     **** Housholder reduction to tridiagonal form ****
C
      call tqr8(w,e,n1,nm,z)
c     **** diagonalize reduced matrix by QL algorithm ****
C
      call esrtr8(w,z,n1,nm )
c     **** sort by eigenvalues ****
c
c
      RETURN
      END
C
C
C
C
C
      SUBROUTINE PRTEIG(W,Z,IERR,IVECT)
C     **** PRINT EIGVALS FROM W, EIGVECTS FROM Z ****
C     **** IF IVECT=0, ONLY PRINT EIGVALS        ****
c     ****  otherwise, ivect is max mumber of evects to print ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      REAL*8 W(6902),Z(6902,6902)
      INTEGER IERR,IVECT
C
      N1 = NUMSTA
     
      N1MX = N1
      IF((IVECT.GT.0).AND.(IVECT.LT.N1)) N1MX = IVECT
c      limit output of eigenvects
C
      WRITE(15,3600) IERR
      WRITE(15,3602)
      WRITE(15,3601) (W(I), I=1,N1)
C
c     IF(IVECT.NE.0) THEN
      IF(IVECT.GT.0) THEN
       WRITE(15,3612)
c      DO 3200 L=1,N1,4
       DO 3200 L=1,N1MX,4
       DO 3100 I=1,N1
        NN=L+3
c       IF(NN.GT.N1) NN=N1
        IF(NN.GT.N1MX) NN=N1MX
 3100   WRITE(15,3611) (Z(I,J), J=L,NN)
        WRITE(15,3620)
 3200  CONTINUE
      ELSE
      ENDIF
C
 3600 FORMAT(///,' IERR =',I6,//)
 3602 FORMAT(' COMPUTED EIGENVALUES:'/)
 3601 FORMAT(3e24.16)
 3612 FORMAT(///,' COMPUTED EIGENVECTORS:',/)
 3611 FORMAT(4e24.16)
 3620 FORMAT(' ------')
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE SLFN
C     **** COMPUTE THE SELF ENERGIES AND STORE IN SELFEN ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/SE/SELFEN(100)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      INTEGER N,NF,B,K
C
      MXSLFN = 100
C     **** MAX MOM IN SELFEN (AN EVEN NUMBER;  ****
C     ****              SEE NEXT TWO DO'S)     ****
C
C     SELFEN(1)=0.5D0
C     SELFEN(1)=0.5D0 + .5D0*BRACK(0,0)
      SELFEN(1)=0.0D0
      DO 10 I=1,MXSLFN-3,2
C      SELFEN(I+1)=SELFEN(I)+.5*(1.D0/(I*I)+1.D0/((I+1.D0)*(I+1.D0)))
       SELFEN(I+2)=SELFEN(I) + 4.D0/((I+1.D0)*(I+1.D0))
   10 CONTINUE
C
      DO 15 I2=1,MXSLFN-1,2
       SELFEN(I2)=((N*N-1.0D0)/N)*SELFEN(I2)
   15 CONTINUE
C
C     DO 69 J=1,25,2
C      RLAMB2=RLAMB*RLAMB
C      SELFEN(J)=SELFEN(J)+((1.0D0-RLAMB2)/RLAMB2 )/J
C  69 CONTINUE
C
C     WRITE(15,20)
C  20 FORMAT(' SELFEN: ')
C     WRITE(15,30)(SELFEN(J),J=1,25)
C  30 FORMAT(3D16.8)
C
      RETURN
      END
C
C
C
C
C     FUNCTION BRACK(L,M)
C     **** COMPUTE THE SQUARE BRACKET OF L AND M ****
C     **** USED IN EVALUATING THE HAMILT.        ****
C     IMPLICIT REAL*8 (A-H,O-Z)
C
C     KAPPA=100
C
C     IF ((L.EQ.0).AND.(M.EQ.0)) THEN
C      BRACK=-KAPPA*KAPPA
C     ELSE
C      IF((L+M).NE.0) THEN
C      IF (((L.EQ.0).OR.(M.EQ.0)).OR.((L+M).NE.0)) THEN
C       BRACK=0.0D0
C       RETURN
C      ELSE
C       BRACK=1.D0/(L*L)
C      ENDIF
C
C     ENDIF
C
C     RETURN
C     END
C
C
C
      FUNCTION BRACK(L,M)
C     **** COMPUTE THE SQUARE BRACKET OF L AND M ****
C     **** USED IN EVALUATING THE HAMILT.        ****
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (((L.EQ.0).OR.(M.EQ.0)).OR.((L+M).NE.0)) THEN
       BRACK=0.0D0
       RETURN
      ELSE
       RL=DFLOAT(L)
C      BRACK=1.0D0/(RL*RL)
       BRACK=4.0D0/(RL*RL)
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE NUZ(W,Z)
C     **** NORMALIZE Z BY DIVIDING BY ****
C     **** SQRTS OF EIGVALUES         ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      REAL*8 W(6902),Z(6902,6902)
C
      DO 10 ML=1,NUMSTA
      DO 10 MR=1,NUMSTA
C      **** DIVIDE HNORM EVECTS BY SQRT OF ****
C      **** EVALS TO NORMALIZE STATES      ****
       Z(ML,MR)=Z(ML,MR)/(W(MR)**(.5d0))
   10 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE NUHAM(Z)
C     **** COMPUTE NEW HAM BASED ON   ****
C     **** ORTHONORMAL STATES         ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/NST/NUMSTA
      COMMON/HM/HAM0(03,6902,6902),HAM(6902,6902)
c     COMMON/HN/HNU0(03,6902,6902),HNU(6902,6902)
      COMMON/HN/HNU0(03,6902),HNU(6902,6902)
      INTEGER N,NF,B,K
C     REAL*8 W(6902),Z(6902,6902)
      REAL*8 Z(6902,6902)
C
      EPS=1.D-10
C
C     DO 50 IN1=1,NUMSTA
C     DO 50 IN2=1,NUMSTA
C      HNU(IN1,IN2)=0.0D0
c     do 50 in3=1,nf
C      HNU0(in3,IN1,IN2)=0.0D0
C  50 CONTINUE
C
C     WRITE(15,5)NUMSTA
C   5 FORMAT(' NUMSTA (IN NUHAM):',I3)
C
C     DO 10 ML=1,NUMSTA
C     DO 10 MR=1,NUMSTA
C      **** DIVIDE HNORM EVECTS BY SQRT OF ****
C      **** EVALS TO NORMALIZE STATES      ****
C      Z(ML,MR)=Z(ML,MR)/(W(MR)**(.5))
C  10 CONTINUE
C
C
c ??????????????????????????????????????
C     **** COMPUTE NEW HAM BASED ON ORTHON. STATES ****
c     **** hnu0(nf): ****
      do 20 ifl=1,nf
      DO 20 NL=1,NUMSTA
c     DO 20 NR=1,NL
       NR=NL
C
c      HLR=0.D0
       HLR0=0.D0
C
       DO 27 KL=1,NUMSTA
        ZL=Z(KL,NL)
        IF(DABS(ZL).GT.EPS) THEN
         DO 25 KR=1,NUMSTA
          ZR=Z(KR,NR)
           IF(DABS(ZR).GT.EPS) THEN
c           HLR =HLR + HAM(KL,KR)*ZL*ZR
c           HLR0=HLR0+ HAM0(KL,KR)*ZL*ZR
            HLR0=HLR0+ HAM0(ifl,KL,KR)*ZL*ZR
           ELSE
           ENDIF
   25    CONTINUE
        ELSE
        ENDIF
   27  CONTINUE
C
c      HNU(NL,NR)=HLR
c      HNU0(ifl,NL,NR)=HLR0
       HNU0(ifl,NL)=HLR0
c      HNU(NR,NL)=HLR
c      HNU0(ifl,NR,NL)=HLR0
C
c      WRITE(15,32)ifl,NL,NR,HNU0(ifl,NL,NR)
c  32  FORMAT(' HNU0(',i3,',',I3,',',I3,')= ',F8.3)
cf     WRITE(15,32)ifl,NL,HNU0(ifl,NL)
cf 32  FORMAT(' HNU0(',i3,',',I3,')= ',F8.3)
C
   20 CONTINUE
c
c     **** hnu: ****
      DO 40 NL=1,NUMSTA
      DO 40 NR=1,NL
C
       HLR=0.D0
c      HLR0=0.D0
C
       DO 47 KL=1,NUMSTA
        ZL=Z(KL,NL)
        IF(DABS(ZL).GT.EPS) THEN
         DO 45 KR=1,NUMSTA
          ZR=Z(KR,NR)
           IF(DABS(ZR).GT.EPS) THEN
            HLR =HLR + HAM(KL,KR)*ZL*ZR
c           HLR0=HLR0+ HAM0(KL,KR)*ZL*ZR
c           HLR0=HLR0+ HAM0(ifl,KL,KR)*ZL*ZR
           ELSE
           ENDIF
   45    CONTINUE
        ELSE
        ENDIF
   47  CONTINUE
C
       HNU(NL,NR)=HLR
c      HNU0(ifl,NL,NR)=HLR0
       HNU(NR,NL)=HLR
c      HNU0(ifl,NR,NL)=HLR0
C
cf     WRITE(15,50)NL,NR,HNU(NL,NR)
cf 50  FORMAT(' HNU(',I3,',',I3,')= ',F8.3)
C
   40 CONTINUE
C
      
      RETURN
      END
C
C
C
cf    SUBROUTINE LPNSUB0(N,B,K,LNB,LND)
C
C     **** GIVEN N,B,K COMPUTE THE MAX NUMBER OF ****
C     **** B-ADJ (LNB) AND D-ADJ (LND) THAT      ****
C     **** CAN APPEAR IN A STATE                 ****
C     **** THIS ROUTINE IS SPECIFICALLY FOR      ****
C     **** ANTI-PBC'S; THE QUANTA CARRY MOMENTA  ****
C     **** WHICH ARE ODD INTEGERS                ****
C      **** THIS IS MODIFIED BELOW TO INCLUDE FLAVOR ****
C
cf    IMPLICIT REAL*8 (A-H,O-Z)
cf    REAL*8 NSUBB,NSUBD,MSUBB,MSUBD
cf    INTEGER N,B,K
C
cf    RN=DFLOAT(N)
cf    RB=DFLOAT(B)
cf    RK=DFLOAT(K)
C
cf    X1=2.0D0*RN*RK - RN*RN*RB*RB
C
cf    Y= .5D0*(RN*RB + X1**.5)
C
C
cf    NSUBB=INT(Y)
cf    NSUBD=NSUBB-RB*RN
cf 70 MSUBB=INT(NSUBB/RN)
cf    MSUBD=INT(NSUBD/RN)
cf    BKMIN=RN*MSUBB*MSUBB + (NSUBB-N*MSUBB)*(2.0D0*MSUBB + 1.0D0)
cf    DKMIN=RN*MSUBD*MSUBD + (NSUBD-N*MSUBD)*(2.0D0*MSUBD + 1.0D0)
C     WRITE(15,50)BKMIN,DKMIN
C  50 FORMAT(' BKMIN= ',G12.5,'  DKMIN= ',G12.5)
cf    IF ((BKMIN+DKMIN).GT.K) THEN
cf     NSUBB=NSUBB-1.0D0
cf     NSUBD=NSUBD-1.0D0
cf     IF((NSUBB.LT.(B*N)).OR.(NSUBD.LT.0.0D0)) THEN
cf      WRITE(15,60)
cf 60   FORMAT(' INSUFFICIENT K ')
cf      STOP
cf     ELSE
cf     ENDIF
cf     GO TO 70
cf    ELSE
cf    ENDIF
C     WRITE(*,80)NSUBB,NSUBD
C     WRITE(15,80)NSUBB,NSUBD
C  80 FORMAT(' NSUBB= ',G12.5,' NSUBD= ',G12.5 )
C
cf    LNB=INT(NSUBB)
cf    LND=INT(NSUBD)
C
cf    RETURN
cf    END
C
C
C
C
      SUBROUTINE LPNSUB(N,NF,B,K,LNB,LND)
C
C     **** GIVEN N,NF,B,K COMPUTE THE MAX NUMBER ****
C     **** OF B-ADJ (LNB) AND D-ADJ (LND) THAT   ****
C     **** CAN APPEAR IN A STATE                 ****
C     **** THIS ROUTINE IS SPECIFICALLY FOR      ****
C     **** ANTI-PBC'S; THE QUANTA CARRY MOMENTA  ****
C     **** WHICH ARE ODD INTEGERS                ****
C      **** MODIFIED FOR FLAVOR ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
c     INTEGER N,NF,B,K,BPR,DPR
      INTEGER N,NF,B,K
C
C     **** ITERATE THROUGH NUMBER OF QUARKS ALLOWED ****
C     **** (NSUBB) UNTIL EXCEED MIN MOMENTUM NECESS.****
C     **** FOR THAT NUMBER OF QUARKS                ****
C
      NBMIN=N*B
      NBMAX=25
C
      DO 10 NSUBB = NBMIN, NBMAX+1
       NSUBD=NSUBB-N*B
C
       IF(NSUBB.EQ.NBMAX+1) THEN
        WRITE(*,20)
   20   FORMAT(' NSUBB EXCEEDS NBMAX IN SUB LPNSUB')
        STOP
       ELSE
       ENDIF
C
       kbmin=minmom(nsubb)
C      **** MIN MOMENTUM CARRIED BY NSUBB QUARKS ****
       kdmin=minmom(nsubd)
C      **** MIN MOMENTUM CARRIED BY NSUBD ANTIQUARKS ****
       KMIN=KDMIN + KBMIN
C      **** MIN MOMENTUM CARRIED BY BOTH ****
C
       IF (KMIN.GT.K) THEN
        LNB = NSUBB -1
        LND = LNB - N*B
        IF (LNB.LT.N*B) THEN
         WRITE(*,30)
   30    FORMAT(' K INSUFFICIENT IN SUB LPNSUB')
         STOP
        ELSE
        ENDIF
        GO TO 5
       ELSE
       ENDIF
C
   10 CONTINUE
C
    5 CONTINUE
C
      RETURN
      END
C
C
C
C
      function minmom(nquant)
c   ???????????????????????????????
c     **** given nquant particles (or antipart), ****
c     **** compute the min momentum they carry.  ****
c     **** assumes antiper. bc's and accounts    ****
c     **** for flavor                            ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      INTEGER N,NF,B,K
c   
      nnf=n*nf
      nefbr=int(nquant/nnf)
      nqrem=nquant-nnf*nefbr      
c     minmom=n*nefbr*nefbr + nqrem*(2*nefbr + 1)
      minmom=nnf*nefbr*nefbr + nqrem*(2*nefbr + 1)
C
      return
      end
C
C
c ???????????????????????????????????????????
      FUNCTION CLFACT(IC1,IC2,IC3,IC4) 
C     **** CALLED BY HAMQCD TO EVALUATE              ****
C     **** COLOR INTERACTION:                        ****
C     ****  (FOR NOPS=4)                             ****
C     **** (1/2)( DELTA(IC1,IC2)*DELTA(IC3,IC4)      ****
C     ****    -(1/N)*DELTA(IC3,IC2)*DELTA(IC1,IC4) ) ****
C     ****                                           ****
C     **** IC1,IC3 MUST BE UPPER INDICES;            ****
C     **** IC2,IC4 LOWER.                            ****
C     ****                                           ****
C     **** IF NOPS IS 2, INTER. IS DELT(IC1,IC2),W/  ****
C     **** IC1 THE UPPER INDEX.                      ****
C     ****                                           ****
C     **** THIS FUNCTION WORKS OUT THE POSSIBLE      ****
C     **** KRON. DELTAS (IN COLOR) FROM CONTRACTING  ****
C     **** OPERATORS FROM STATES AND INTERACTIONS    ****
C     **** AND THEN DIAGRAMMATICALLY CONTRACTS THEM  ****
C     **** WITH THE DELTAS AND EPSILONS FROM THE     ****
C     **** MESONIC AND BARYONIC SINGLET STATES       ****
C     ****                                           ****
C
C     NEED TO GET MX,K,LLT,LRT,NOPS,N,B
C     FROM CALLING ROUTINE
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/GL1/LPO
      COMMON/GLO/RESL0,IDEL0
      COMMON/CLR/RESLT,IDELT
      COMMON/LN/LNG,NTRMS
C     COMMON/BAR/LNKB(25)
      COMMON/BAR/LNKB(25),ICYCL(25)
      COMMON/ME/MX,LLT,LRT,NOPS
      COMMON/CLFDIAG/ICLDI
      INTEGER IDEL0(12552,25),IDELT(12552,25)
      REAL*8 RESL0(12552),RESLT(12552)
      INTEGER MX(4,25)
      INTEGER MY(4,25)
      INTEGER NCR(2,13)
      INTEGER NPERM(8,25)
c     INTEGER MOCHK(2,2,50)
      INTEGER MOCHK(2,2,50,03)
      INTEGER N,NF,B,K
C     **** THIRD INDEX OF MOCHK IS MAX MOM; ****
c     **** fourth index is flavor ****
c
      ICLDI=0
C     **** SET ICLDI NE 0 TO GET DIAGNOSTIC OUTPUT ****
C     **** FROM CLFACT AND RELATED SUBROUTINES     ****
C
c
      maxfl=03
      if(nf.gt.maxfl)then
       WRITE(15,2000)
 2000  FORMAT(' nf exceeds maxfl in sub clfact')
       stop
      else
      endif
      MXLNG=25
c     maxk=50
c      ******* 99 ********
      maxk=130
      if(k.gt.maxk) then
       WRITE(15,2010)
 2010  FORMAT(' k exceeds maxk in sub clfact')
       stop
      else
      endif
      MOMAX=K
      MXTRMS=12552
      NMSR=(LRT-N*B)/2
      NMSL=(LLT-N*B)/2
      LNG = NOPS+LLT+LRT
      IF(LNG.GT.MXLNG) THEN
       WRITE(15,7)
    7  FORMAT(' LNG EXCEEDED MXLNG IN FN CLFACT ')
       STOP
      ELSE
      ENDIF
C
C     **** NUMBER OF RIGHT, LEFT MESONS ****
      IF(ICLDI.NE.0) THEN
       WRITE(15,15)NMSR,NMSL
   15  FORMAT(/' NUM OF RT,LT MESONS:(NMSR,NMSL) ',2I6)
      ELSE
      ENDIF
C
      DO 10 I1=1,MXTRMS
       RESL0(I1)=1.D0
       RESLT(I1)=1.D0
       DO 10 I2=1,LNG
        IDEL0(I1,I2)=0
        IDELT(I1,I2)=0
   10 CONTINUE
C     **** IDEL0 CONTAINS THE DELTA FUNCTIONS FROM ****
C     **** CONTRACTING OPERATORS; IDELT HAS THE    ****
C     **** ADDITIONAL DELTAS FROM THE INTERACTIONS ****
C     **** AND MESONIC PAIRS IN STATES;            ****
C     **** RESL0, RESLT CONTAIN THE CORRESPONDING  ****
C     **** COEFFICIENT(INCLUDING SIGN) FROM        ****
C     **** CONTRACTING THEM                        ****
C
C
      DO 5 I0=1,LNG
       LNKB(I0)=0
       ICYCL(I0)=1
    5 CONTINUE
C
C
C
C
C     **** CHECK IF ELEMENT IS ZERO: SEE IF CAN MATCH ****
C     **** UP THE MOMENTA CARRIED BY CREATION AND     ****
C     **** ANNIH. OPERATORS                           ****
      DO 40 IMX=1,MOMAX+1
      DO 40 I2=1,2
      DO 40 J2=1,2
      DO 40 k2=1,nf
       MOCHK(I2,J2,IMX,k2)=0
   40 CONTINUE
C
C
      DO 20 K1=1,4
      DO 20 L1=1,LNG
       MY(K1,L1)=MX(K1,L1)
   20 CONTINUE
C
      IF(ICLDI.NE.0) THEN
       WRITE(15,35)
   35  FORMAT( ' MX: ')
       DO 60 I=1,4
        WRITE(15,55)(MX(I,J),J=1,LNG)
   55   FORMAT(25I6)
   60  CONTINUE
      ELSE
      ENDIF
C
C     ISGNC=1
      DO 90 I3=1,LNG
      DO 90 K3=1,2
C
      MOCHK(K3,1,MX(3,I3)+1,mx(4,i3))=
     >  MOCHK(K3,1,MX(3,I3)+1,mx(4,i3))+MX(K3,I3)
C
      MOCHK(K3,2,MX(3,LNG+1-I3)+1,mx(4,lng+1-i3))=
     >  MOCHK(K3,2,MX(3,LNG+1-I3)+1,mx(4,lng+1-i3))+
     >  MX(K3,LNG+1-I3)
C
c     **** check if too many annih or create ops   ****
c     **** moving either from rt to lt or lt to rt ****
c     IF((MOCHK(K3,1,MX(3,I3)+1).LT.0).OR.
c    > (MOCHK(K3,1,MX(3,I3)+1).GT.N).OR.
c    > (MOCHK(K3,2,MX(3,LNG+1-I3)+1).GT.0).OR.
c    > (MOCHK(K3,2,MX(3,LNG+1-I3)+1).LT.-N)) THEN
      IF((MOCHK(K3,1,MX(3,I3)+1,mx(4,i3)).LT.0).OR.
     > (MOCHK(K3,1,MX(3,I3)+1,mx(4,i3)).GT.N).OR.
     > (MOCHK(K3,2,MX(3,LNG+1-I3)+1,mx(4,i3)).GT.0).OR.
     > (MOCHK(K3,2,MX(3,LNG+1-I3)+1,mx(4,i3)).LT.-N)) THEN
C      ISGNC=0
       CLFACT=0.D0
       IF(ICLDI.NE.0) THEN
        WRITE(15,80)
   80   FORMAT( ' MATRIX=0 ')
       ELSE
       ENDIF
C      GO TO 5000
       RETURN
      ELSE
      ENDIF
   90 CONTINUE
C
C
C
C     **** CREATE MATRIX LNKB WHICH STORE INFO ON  ****
C     **** HOW THE BARYON INDICES ARE LINKED BY    ****
C     **** THE ANTISYMM TENSOR;                    ****
C     **** NOTE THAT THIS DEPENDS ON FILLING THE   ****
C     **** STATES FIRST WITH B BARYONS, THEN WITH  ****
C     **** MESONS; ALSO THAT B IS THE SAME         ****
C     **** IN RIGHT AND LEFT(ADJOINT) STATES       ****
C
      IF(B.NE.0) THEN
       DO 1500 NBS=1,B
C
        DO 1510 INB=N*(NBS-1)+1,N*NBS-1
         LNKB(INB)=INB+1
 1510   CONTINUE
        LNKB(N*NBS)=N*(NBS-1)+1
C
        DO 1520 JNB=LNG-N*(NBS-1),LNG-N*NBS+2 ,-1
         LNKB(JNB)=JNB-1
 1520   CONTINUE
        LNKB(LNG-N*NBS+1)=LNG-N*(NBS-1)
C
 1500  CONTINUE
C
C
       IF(((-1)**N).GT.0) THEN
C      **** ICYCL WILL TAKE INTO ACCOUNT THAT       ****
C      **** FOR EVEN N, EPS CHANGES SIGN FOR CYCLIC ****
C      **** PERMUTATIONS; NEED WHEN CONVERTING      ****
C      **** EPSILON PAIRS INTO PRODUCTS OF DELTAS   ****
        DO 1503 NBS=1,B
C
         II=1
         DO 1512 INB=N*(NBS-1)+1,N*NBS
          ICYCL(INB)=II
          II=II*(-1)
 1512    CONTINUE
C
         JJ=1
         DO 1522 JNB=LNG-N*(NBS-1),LNG-N*NBS+1 ,-1
          ICYCL(JNB)=JJ
          JJ=JJ*(-1)
 1522    CONTINUE
C
 1503   CONTINUE
       ELSE
       ENDIF
C
      ELSE
      ENDIF
C
C
C
C     **** CONSTRUCT THE INITIAL SET OF DELTAS ****
C     **** BY CONTRACTING THE ANNIH. OPERATORS ****
C     **** WITH THE NEAREST APPROPRIATE CREAT. ****
C     **** OPERATOR; KEEP TRACK OF SIGN; STORE ****
C     **** IN IDEL0(1, )                       ****
C
      IPOW=0
      DO 400 JA=2,LNG
      IF ((MY(1,JA)+MY(2,JA)).NE.-1) GO TO 400
      DO 300 KA=JA-1,1,-1
c     IF ((MY(1,KA).EQ.-MY(1,JA)).AND.(MY(2,KA).EQ.-MY(2,JA)).AND.
c    >    (MY(3,KA).EQ.MY(3,JA))) THEN
      IF ((MY(1,KA).EQ.-MY(1,JA)).AND.(MY(2,KA).EQ.-MY(2,JA)).AND.
     >    (MY(3,KA).EQ.MY(3,JA)).and.(my(4,ka).eq.my(4,ja))) THEN
       IDEL0(1,JA)=KA
       IDEL0(1,KA)=JA
       MY(1,KA)=0
       MY(2,KA)=0
       MY(1,JA)=0
       MY(2,JA)=0
       GO TO 400
      ELSE
       IF((MY(1,KA).NE.0).OR.(MY(2,KA).NE.0)) THEN
        IPOW=IPOW+1
       ELSE
       ENDIF
      END IF
  300 CONTINUE
  400 CONTINUE
C
      ISGNC=(-1)**IPOW
      RESL0(1)=DFLOAT(ISGNC)
      IF(ICLDI.NE.0) THEN
       WRITE(15,450)ISGNC
  450  FORMAT( ' ISGNC= ',I5)
       WRITE(15,475)
  475  FORMAT( ' RESL0(1),IDEL0(1, ): ')
       WRITE(15,550)RESL0(1),(IDEL0(1,L),L=1,LNG)
  550  FORMAT(e24.16,25I5)
      ELSE
      ENDIF
C
C
C     **** CONSTRUCT THE SET OF PERMUTATIONS WHICH ****
C     **** WILL GENERATE THE COMPLETE SET OF DELTA ****
C     **** FUNCTIONS POSSIBLE CONSISTENT WITH THE  ****
C     **** MOMENTA CARRIED BY CR. AND ANN. OPS;    ****
C     **** STORE THE SET OF PERMUTATIONS IN NPERMS ****
      LBCR=0
      LDCR=0
      NTRMS=1
      LPRM=0
      LPO=1
      DO 625 IA=1,LNG/2
      DO 615 IB=1,2
      NCR(IB,IA)=0
  615 CONTINUE
      DO 625 IC=1,LNG
      NPERM(IA,IC)=0
  625 CONTINUE
      DO 600 I=1,LNG
      IF(MX(2,I).EQ.0) THEN
C     **** B OPERATOR ****
        IF(MX(1,I).EQ.-1) THEN
          MP=1
          LPRM=LPRM+1
          LPCHK=0
          DO 650 J=LBCR,1,-1
           IF(LPCHK.EQ.0) THEN
c           IF(MX(3,NCR(1,J)).EQ.MX(3,I)) THEN
            IF( (MX(3,NCR(1,J)).EQ.MX(3,I)).and.
     >          (MX(4,NCR(1,J)).EQ.MX(4,I)) ) THEN
             JFSTB=J
             IF(J.NE.1) THEN
              DO 660 L=J-1,1,-1
c              IF(MX(3,NCR(1,L)).EQ.MX(3,NCR(1,J))) THEN
               IF( (MX(3,NCR(1,L)).EQ.MX(3,NCR(1,J))).and.
     >             (mx(4,ncr(1,l)).eq.mx(4,ncr(1,j))) ) then 
                NPERM(LPRM,1)=NPERM(LPRM,1)+1
                MP=2*NPERM(LPRM,1)
                NPERM(LPRM,MP)=NCR(1,J)
                NPERM(LPRM,MP+1)=NCR(1,L)
                LPCHK=1
               ELSE
               ENDIF
  660         CONTINUE
             ELSE
             ENDIF
            ELSE
            ENDIF
           ELSE
           ENDIF
  650     CONTINUE
          IF(LPCHK.EQ.0) THEN
C          **** NO MOMENTUM MATCH ****
           LPRM=LPRM-1
          ELSE
          ENDIF
C         **** REMOVE FROM NCR LOC OF B-ADJ W/ ****
C         **** SAME MOM AS B AT LOC I          ****
          DO 800 M=JFSTB,LBCR-1
           NCR(1,M)=NCR(1,M+1)
  800     CONTINUE
          NCR(1,LBCR)=0
          LBCR=LBCR-1
        ELSE
         LBCR=LBCR+1
         NCR(1,LBCR)=I
        ENDIF
      ELSE
C     **** D OPERATOR ****
        IF(MX(2,I).EQ.-1) THEN
          MP=1
          LPRM=LPRM+1
          LPCHK=0
          DO 680 J=LDCR,1,-1
           IF(LPCHK.EQ.0) THEN
c           IF(MX(3,NCR(2,J)).EQ.MX(3,I)) THEN
            IF( (MX(3,NCR(2,J)).EQ.MX(3,I)).and.
     >          (MX(4,NCR(2,J)).EQ.MX(4,I)) ) THEN
             JFSTD=J
             IF(J.NE.1) THEN
              DO 690 L=J-1,1,-1
c              IF(MX(3,NCR(2,L)).EQ.MX(3,NCR(2,J))) THEN
               IF( (MX(3,NCR(2,L)).EQ.MX(3,NCR(2,J))).and.
     >             (mx(4,ncr(2,l)).eq.mx(4,ncr(2,j))) ) then 
                NPERM(LPRM,1)=NPERM(LPRM,1)+1
                MP=2*NPERM(LPRM,1)
                NPERM(LPRM,MP)=NCR(2,J)
                NPERM(LPRM,MP+1)=NCR(2,L)
                LPCHK=1
               ELSE
               ENDIF
  690         CONTINUE
             ELSE
             ENDIF
            ELSE
            ENDIF
           ELSE
           ENDIF
  680     CONTINUE
          IF(LPCHK.EQ.0) THEN
C          **** NO MOMENTUM MATCH ****
           LPRM=LPRM-1
          ELSE
          ENDIF
C         **** REMOVE FROM NCR LOC OF D-ADJ W/ ****
C         **** SAME MOM AS D AT LOC I          ****
          DO 850 M=JFSTD,LDCR-1
           NCR(2,M)=NCR(2,M+1)
  850     CONTINUE
          NCR(2,LDCR)=0
          LDCR=LDCR-1
        ELSE
         LDCR=LDCR+1
         NCR(2,LDCR)=I
        ENDIF
       ENDIF
  600 CONTINUE
C
C
      IF(ICLDI.NE.0) THEN
       WRITE(15,1019)
 1019  FORMAT(' NPERM: ')
       DO 1015 IPM=1,LNG/2
        WRITE(15,1010)(NPERM(IPM,M),M=1,LNG)
 1010   FORMAT(25I6)
 1015  CONTINUE
      ELSE
      ENDIF
C
      IF(LPRM.GT.0) THEN
        DO 750 ID=LPRM,1,-1
         NSBAR=1
         DO 775 JD=1,NPERM(ID,1)
          KD=2*JD
          IP1=NPERM(ID,KD)
          IP2=NPERM(ID,KD+1)
C
          ISMB=0
          IF(B.GT.0) THEN
C          **** SEE IF PERMD INDICES ARE IN SAME BARYON;  ****
C          **** IF SO, SIMPLY MULT EXISTING TERM IN IDELT ****
C          **** BY APPROPR. FACTOR                        ****
           NXTB=LNKB(IP1)
           IF(NXTB.NE.0) THEN
C          **** IP1 NOT IN A MESON ****
            DO 1122 JBS=1,N-1
             IF(IP2.EQ.NXTB) THEN
C             **** IN SAME BARYON ****
              ISMB=1
             ELSE
             ENDIF
             NXTB=LNKB(NXTB)
 1122       CONTINUE
           ELSE
           ENDIF
          ELSE
          ENDIF
C
          IF(ISMB.EQ.1) THEN
           NSBAR=NSBAR+1
C          **** NUMBER OF PERMS IN SAME BARYON ****
           IF(ICLDI.NE.0) THEN
            WRITE(15,795)IP1,IP2
  795       FORMAT('  ',2I5,' ARE IN THE SAME BARYON')
            WRITE(15,797)NSBAR
  797       FORMAT(' NSBAR:  ',I5)
           ELSE
           ENDIF
          ELSE
C          **** CREATE NEW TERMS IN IDELT BY PERMUTING ****
           IF(ICLDI.NE.0) THEN
            WRITE(15,780)NTRMS
  780       FORMAT(' NTRMS(BEFORE CALL NWTERM): ',I5)
           ELSE
           ENDIF
           CALL NWTERM(IP1,IP2)
C          CALL NWTERM(NPERM(ID,KD),NPERM(ID,KD+1))
          ENDIF
C
  775    CONTINUE
         IF(NSBAR.GT.1) THEN
          DO 790 IMLT=1,NTRMS
           RESL0(IMLT)=RESL0(IMLT)*DFLOAT(NSBAR)
           IF(ICLDI.NE.0) THEN
            WRITE(15,796)IMLT,RESL0(IMLT)
  796       FORMAT(' RESL0(',I6,')= ',e24.16)
           ELSE
           ENDIF
  790     CONTINUE
         ELSE
         ENDIF
         NTRMS=NTRMS*(NPERM(ID,1)+2-NSBAR)
C        NTRMS=NTRMS*(NPERM(ID,1)+1)
C        WRITE(15,782)NTRMS
C 782    FORMAT(' NTRMS(AFTER CALL NWTERM): ',I5)
         IF(NTRMS.GT.MXTRMS) THEN
          WRITE(15,785)
  785     FORMAT(/' NTRMS EXCEEDED MXTRMS IN FN CLFACT  ')
          STOP
         ELSE
         ENDIF
  750   CONTINUE
      ELSE
      ENDIF
C
C
      IF(ICLDI.NE.0) THEN
C      **** PRINT DELTA FUNCTIONS GENERATED ****
       WRITE(15,1260)
 1260  FORMAT(//' RESL0, IDEL0 : ')
       WRITE(15,1230)NTRMS
 1230  FORMAT(/' THERE ARE',I6,' TERMS'/)
       DO 1200 I=1,NTRMS
        WRITE(15,1250)RESL0(I),(IDEL0(I,J),J=1,LNG)
 1250   FORMAT(e24.16,25I5)
 1200  CONTINUE
      ELSE
      ENDIF
C
      NTRMS0=NTRMS
C     **** SAVE NUM TERMS IN IDEL0 ****
C
C
C
C
C     **** SUM OVER TWO PARTS OF COLOR FACTOR: ****
      IF(NOPS.EQ.4)THEN
       NCTF=2
C      **** INCLUDE BOTH COLOR INT. TERMS ****
      ELSE
       NCTF=1
      ENDIF
      CLFACT=0.D0
      DO 1410 NCT=1,NCTF
       IF(ICLDI.NE.0) THEN
        WRITE(15,1005)NCT
 1005   FORMAT(////' COLOR TERM(NCT): ',I6)
       ELSE
       ENDIF
C
C
C     ***** READ IDEL0 INTO IDELT SUCH THAT      ****
C     ***** UPPER COLOR INDICES CORRESPOND TO    ****
C     ***** LOCATION I IN IDELT(I); LOWER INDEX  ****
C     ***** IS THE NUMBER 'IDELT(I)' AT THAT     ****
C     ***** LOCATION; DROP OTHER TERMS           ****
      DO 1100 J=1,LNG
       IF((MX(1,J)-MX(2,J)).EQ.1) THEN
C      **** J IS AN UPPER INDEX ****
        DO 1110 ICL=1,NTRMS
         IDELT(ICL,J)=IDEL0(ICL,J)
 1110   CONTINUE
       ELSE
       ENDIF
 1100 CONTINUE
C
      DO 1120 M1=1,NTRMS
       RESLT(M1)=RESL0(M1)
 1120 CONTINUE
C
C
      IF(NOPS.EQ.4) THEN
C      **** FILL IN DELTAS FROM COLOR VERTEX (IC1,...,IC4) ****
       IF(NCT.EQ.1) THEN
        DO 1400 JC=1,NTRMS
         IDELT(JC,LRT+IC1)=LRT+IC2
         IDELT(JC,LRT+IC3)=LRT+IC4
         RESLT(JC)=.5D0*RESLT(JC)
 1400   CONTINUE
C       **** U(N) PART OF INTERACTION ****
C
       ELSE
        DO 1405 JC=1,NTRMS
         IDELT(JC,LRT+IC1)=LRT+IC4
         IDELT(JC,LRT+IC3)=LRT+IC2
         RESLT(JC)=(-.5D0/DFLOAT(N))*RESLT(JC)
 1405   CONTINUE
C       **** U(1) PART OF INTERACTION ****
       ENDIF
C
      ELSE
       IF(NOPS.EQ.2) THEN
C       **** DELTA FN IN COLOR; LOC OF INT INDICES IN IC1,2 ****
        DO 1408 JC=1,NTRMS
         IDELT(JC,LRT+IC1)=LRT+IC2
 1408   CONTINUE
C
       ELSE
C       **** NO INTERACTION TERM ****
       ENDIF
      ENDIF

C
C
      IF(ICLDI.NE.0) THEN
C     **** PRINT DELTA FUNCTIONS GENERATED ****
C     **** INCLUDING THOSE DUE TO INTERAC. ****
       WRITE(15,960)
  960  FORMAT(//' RESLT, IDELT :(AFTER ADD INTER. DELTAS) ')
       WRITE(15,930)NTRMS
  930  FORMAT(/' THERE ARE',I6,' TERMS'/)
       DO 900 I=1,NTRMS
        WRITE(15,950)RESLT(I),(IDELT(I,J),J=1,LNG)
  950   FORMAT(e24.16,25I5)
  900  CONTINUE
      ELSE
      ENDIF
C
C
C     **** FILL IN DELTAS FROM MESONS;             ****
C     **** NOTE THAT THIS DEPENDS ON FILLING THE   ****
C     **** STATES FIRST WITH B BARYONS, THEN WITH  ****
C     **** MESONS SUCH THAT IN EACH MESON THE      ****
C     **** QUARK PRECEDES THE ANTIQUARK            ****
C
      IF(NMSR.NE.0) THEN
       DO 1300 IL=1,NTRMS
        DO 1310 IM=N*B+2,LRT,2
         IDELT(IL,IM)=IM-1
 1310   CONTINUE
 1300  CONTINUE
      ELSE
      ENDIF
C
      IF(NMSL.NE.0) THEN
       DO 1330 IL=1,NTRMS
        DO 1320 IN=LNG-N*B,LNG-LLT+2,-2
         IDELT(IL,IN)=IN-1
 1320   CONTINUE
 1330  CONTINUE
      ELSE
      ENDIF
C
C
C
C
      IF(ICLDI.NE.0) THEN
C     **** PRINT DELTA FUNCTIONS GENERATED ****
C     **** INCLUDING THOSE DUE TO INTERAC. ****
C     **** AND MESON DELTAS                ****
       WRITE(15,1760)
 1760  FORMAT(//' RESLT, IDELT :(AFTER ADD MES DELTAS) ')
       WRITE(15,1730)NTRMS
 1730  FORMAT(/' THERE ARE',I6,' TERMS'/)
       DO 1700 I=1,NTRMS
        WRITE(15,1750)RESLT(I),(IDELT(I,J),J=1,LNG)
 1750   FORMAT(e24.16,25I5)
 1700  CONTINUE
      ELSE
      ENDIF
C
C
C
C     **** CALL CLSUM TO CONTRACT DELTAS AND EPSILONS  ****
C     **** FROM IDELT AND LNKB AND RETURN COLOR FACTOR ****
C     **** IN CLRTOT                                   ****
C
      CALL CLSUM(CLRTOT)
      IF(ICLDI.NE.0) THEN
       WRITE(15,1610)CLRTOT
 1610  FORMAT(//' INTERIM COLOR FACTOR (CLRTOT):  ',e24.16)
      ELSE
      ENDIF
      CLFACT=CLFACT+CLRTOT
C
C
C     **** RE-SET NTRMS TO NUMBER OF TERMS IN IDELT0 ****
C     **** NEEDED FOR SECOND COLOR TERM              ****
      NTRMS=NTRMS0
C
C     **** END OF SUM OVER TWO PARTS OF COLOR FACTOR ****
 1410 CONTINUE
C
C
 5000 RETURN
      END
C
C
C
      SUBROUTINE NWTERM(LCR,JPR)
C     **** GIVEN INDICES LCR AND JPR,               ****
C     **** PERMUTE THESE IN ALL PREVIOUS SETS OF    ****
C     **** DELTAS AND ADD NEW TERMS TO IDEL0 AND    ****
C     **** APPROPRIATE SIGN TO RESL0; CURRENT NUM   ****
C     **** OF TERMS IN NTRMS                        ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/GL1/LPO
      COMMON/GLO/RESL0,IDEL0
C     COMMON/CLR/RESLT,IDELT
      COMMON/LN/LNG,NTRMS
      COMMON/CLFDIAG/ICLDI
      INTEGER IDEL0(12552,25)
      REAL*8 RESL0(12552)
C
C     WRITE(15,815)
C 815 FORMAT( ' CHECK ')
      IF(ICLDI.NE.0) THEN
       WRITE(15,817)LCR,JPR
  817  FORMAT(' LCR, JPR: ',2I6)
       WRITE(15,822)NTRMS
  822  FORMAT(' NTRMS= ',I5)
      ELSE
      ENDIF
      DO 800 J=1,NTRMS
       LPO=LPO+1
       DO 850 L=1,LNG
        IDEL0(LPO,L)=IDEL0(J,L)
C       WRITE(15,833)LPO,L,IDEL0(LPO,L)
C 833   FORMAT( ' IDEL0(',I3,';',I3,')=',I3)
  850  CONTINUE
      RESL0(LPO)=(-1.D0)*RESL0(J)
C     WRITE(15,860)J,LPO,RESL0(J),RESL0(LPO)
C 860 FORMAT(' RESL0(',I4,'), RESL0(',I4,'): ',2e24.16)
      IDEL0(LPO,LCR)=IDEL0(J,JPR)
      IDEL0(LPO,JPR)=IDEL0(J,LCR)
      IDEL0(LPO,IDEL0(J,JPR))=IDEL0(J,IDEL0(J,LCR))
      IDEL0(LPO,IDEL0(J,LCR))=IDEL0(J,IDEL0(J,JPR))
  800 CONTINUE
      RETURN
      END
C
C
C
C
      SUBROUTINE CLSUM(CLRTOT)
C     **** GIVEN A SET OF DELTAS AND EPSILONS AT        ****
C     **** EACH IDELT(I, ), CONTRACT THE PRODUCT        ****
C     **** OF DELTAS, CONVERT CONTRACTED EPS. INTO      ****
C     **** MORE PRODS. OF DELTAS, CONTRACT, ETC. UNTIL  ****
C     **** REDUCED TO COEFFICIENT, STORED IN RESLT(I);  ****
C     **** SUM THESE COEFFS. AND RETURN SUM IN CLRTOT   ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/CLR/RESLT,IDELT
      COMMON/LN/LNG,NTRMS
      COMMON/IDC/IDCHK(12552),IDCT
C     COMMON/BAR/LNKB(25)
      COMMON/BAR/LNKB(25),ICYCL(25)
      COMMON/CLFDIAG/ICLDI
      INTEGER IDELT(12552,25)
      REAL*8 RESLT(12552)
      INTEGER N,NF,B,K
C
      MXTRMS=12552
      MXLNG=25
C
      IF((B.NE.0).AND.(N.GT.2)) THEN
       CALL EPSB(N,1)
      ELSE
      ENDIF
C
      IDCT=1
      DO 10 I1=1,MXTRMS
       IDCHK(I1)=1
   10 CONTINUE
C
      IF(ICLDI.NE.0) THEN
       IF(B.NE.0) THEN
        WRITE(15,70)
   70   FORMAT(//' BARYON INDICES (LNKB): ')
        WRITE(15,80)(LNKB(J),J=1,LNG)
        WRITE(15,80)(ICYCL(J),J=1,LNG)
   80   FORMAT(/20I5)
       ELSE
       ENDIF
      ELSE
      ENDIF
C
C
      IDCT=1
C     MXITER=200
      MXITER=2*B+1
      DO 30 ITER=1,MXITER
       IF(IDCT.NE.0) THEN
C      **** THERE ARE STILL TERMS LEFT IN IDELT ****
C
        CALL CNTRCT
C       **** CONTRACT DELTAS ****
C
        IF(B.NE.0) THEN
C       **** WILL NEED TO CONVERT EPSILONS TO DELTAS ****
         CALL BREDCE
        ELSE
C       **** WILL ONLY NEED ONE PASS (TO CONTRACT DELTAS)****
         IDCT=0
        ENDIF
C
        IF(ICLDI.NE.0) THEN
         WRITE(15,20)NTRMS
   20    FORMAT(///' NUMBER OF TERMS (NTRMS): ',I6)
        ELSE
        ENDIF
C
       ELSE
       ENDIF
   30 CONTINUE
C
        CLRTOT=0.D0
        DO 50 ISUM=1,NTRMS
         CLRTOT=CLRTOT + RESLT(ISUM)
   50   CONTINUE
        IF(ICLDI.NE.0) THEN
         WRITE(15,60)CLRTOT
   60    FORMAT(' TOTAL COLOR SUM (CLRTOT): ',e24.16)
        ELSE
        ENDIF
C
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE CNTRCT
C     **** CONTRACT DELTA FUNCTIONS ****
C     **** MULTIPLY COEFF. IN RESLT ****
C     **** BY N WHEN INDICES CONTR. ****
C     **** INTO A CLOSED COL. LOOP  ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/CLR/RESLT(12552),IDELT(12552,25)
      COMMON/LN/LNG,NTRMS
      COMMON/IDC/IDCHK(12552),IDCT
      COMMON/CLFDIAG/ICLDI
      INTEGER N,NF,B,K
C
      IF(ICLDI.NE.0) THEN
       WRITE(15,105)
  105  FORMAT(///' SUB CNTRCT: ')
C
       WRITE(15,195)
  195  FORMAT(/' BEFORE CONTRACTION: ')
       DO 200 IO1=1,NTRMS
        WRITE(15,210)RESLT(IO1),(IDELT(IO1,L),L=1,LNG)
  210   FORMAT(e24.16,20I5)
  200  CONTINUE
      ELSE
      ENDIF
C
      DO 30 NT=1,NTRMS
       DO 40 LT=1,LNG
        ICURR=LT
   50   INXT=IDELT(NT,ICURR)
        IF(INXT.NE.0) THEN
         INXT2=IDELT(NT,INXT)
         IF(INXT2.NE.0) THEN
          IDELT(NT,ICURR)=INXT2
          IDELT(NT,INXT)=0
         ELSE
          GO TO 40
         ENDIF
        ELSE
         GO TO 40
        ENDIF
       IF (INXT2.EQ.ICURR) THEN
        IDELT(NT,ICURR)=0
        RESLT(NT)=DFLOAT(N)*RESLT(NT)
       ELSE
       ENDIF
C      ICURR=INXT2
       GO TO 50
   40  CONTINUE
   30 CONTINUE
C
      IF(ICLDI.NE.0) THEN
       WRITE(15,115)
  115  FORMAT(/' AFTER CONTRACTION: ')
       DO 100 IO1=1,NTRMS
        WRITE(15,110)RESLT(IO1),(IDELT(IO1,L),L=1,LNG)
  110   FORMAT(e24.16,20I5)
  100  CONTINUE
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
C
      SUBROUTINE BREDCE
C     **** REDUCE A PAIR OF EPSILONS WITH ONE ****
C     **** CONTRACTED INDEX TO AN ANTISYM-    ****
C     **** METRIZED SUM OF PRODUCTS OF DELTAS;****
C     **** ADD THESE NEW TERMS TO IDELT; THESE****
C     **** WILL THEN BE CONTRACTED IN SUBR.   ****
C     **** CNTRCT.                            ****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/DLPM/IBRPM(9523,25),NPRMS
      COMMON/CLR/RESLT(12552),IDELT(12552,25)
      COMMON/LN/LNG,NTRMS
      COMMON/IDC/IDCHK(12552),IDCT
C     COMMON/BAR/LNKB(25)
      COMMON/BAR/LNKB(25),ICYCL(25)
      COMMON/CLFDIAG/ICLDI
      INTEGER N,NF,B,K
      INTEGER NB1(25),NB2(25)
C
      MXTRMS=12552
C     **** MAX NUMBER OF TERMS IN IDELT ****
C
      NTPR=NTRMS
      IDCT=0
      DO 10 NT=1,NTRMS
       IF(IDCHK(NT).NE.0) THEN
C      **** THERE ARE STILL TERMS LEFT IDELT(NT) ****
        IDCHK(NT)=0
        DO 20 LT=1,LNG
         IF (IDELT(NT,LT).NE.0) THEN
          IDCHK(NT)=1
          IDCT=1
          INB1=LT
          INB2=IDELT(NT,LT)
C          WRITE(15,24)NT,RESLT(NT)
C  24      FORMAT(' NT,RESLT(NT): ',I4,2X,D10.4)
          RESLT(NT)=RESLT(NT)*ICYCL(INB1)*ICYCL(INB2)
C          WRITE(15,27)NT,INB1,INB2,ICYCL(INB1),ICYCL(INB2)
C  27      FORMAT(' NT, INB1,2; ICYCL(INB1,2): ',5I4)
C          WRITE(15,24)NT,RESLT(NT)
C         **** ACCOUNTS FOR CHANGE IN SIGN UNDER  ****
C         **** CYCLIC PERMS OF EPS (FOR EVEN N)   ****
          DO 30 I1=1,N-1
           NB1(I1)=LNKB(INB1)
           NB2(I1)=LNKB(INB2)
           INB1=LNKB(INB1)
           INB2=LNKB(INB2)
   30     CONTINUE
          IDELT(NT,LT)=0
          DO 40 J1=1,N-1
           IDELT(NT,NB2(J1))=NB1(J1)
   40     CONTINUE
          DO 50 J2=2,NPRMS
           NTPR=NTPR+1
           IF(NTPR.GT.MXTRMS) THEN
            WRITE(15,213)NTPR
  213       FORMAT(/' NTPR, ',I8,', EXCEEDED MXTRMS IN SUB BREDCE')
            STOP
           ELSE
           ENDIF
           RESLT(NTPR)=RESLT(NT)*DFLOAT(IBRPM(J2,N))
C           WRITE(15,42)NTPR,RESLT(NTPR)
C  42       FORMAT(' NTPR,RESLT(NTPR): ',I4,2X,D10.4)
           DO 60 J3=1,LNG
            IDELT(NTPR,J3)=IDELT(NT,J3)
   60      CONTINUE
           DO 70 J4=1,N-1
            IDELT(NTPR,NB2(IBRPM(J2,J4)))=NB1(J4)
   70      CONTINUE
   50     CONTINUE
          GO TO 10
         ELSE
         ENDIF
C
   20   CONTINUE
       ELSE
       ENDIF
   10 CONTINUE
C
      NTRMS=NTPR
C
      IF(NTRMS.GT.MXTRMS) THEN
       WRITE(15,211)
  211  FORMAT(/' NTRMS EXCEEDED MXTRMS IN SUB BREDCE ')
       STOP
      ELSE
      ENDIF
C
      IF(ICLDI.NE.0) THEN
       WRITE(15,105)
  105  FORMAT(///' SUB BREDCE: ')
       DO 100 IO1=1,NTRMS
        WRITE(15,110)RESLT(IO1),(IDELT(IO1,L),L=1,LNG)
  110   FORMAT(e24.16,20I5)
  100  CONTINUE
      ELSE
      ENDIF
C
      RETURN
      END
C
C
C
C
      SUBROUTINE EPSB(N,B)
C     **** GENERATE PERMS,INCL SIGN, OF COLOR  ****
C     **** INDEX FOR BARYONS,LENGTH N-1        ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/DLPM/IBRPM(9523,25),NPRMS
      COMMON/CLFDIAG/ICLDI
      INTEGER NMIN(25),NMAX(25),IPRM(25)
C     **** DIMS ARE MAX N-1 AND MAX N-1 FACT. ****
C
      INTEGER B
C
      MAXEP=9523
C     **** DIM OF IBRPM ****
C
      INDX=1
      NPRMS=1
      IF(ICLDI.NE.0) THEN
       WRITE(15,10)N
   10  FORMAT(' N= ',I4)
      ELSE
      ENDIF
C
      DO 40 II=1,25
       NMIN(II)=1
       NMAX(II)=N-1
      DO 40 JJ=1,MAXEP
       IBRPM(JJ,II)=0
   40 CONTINUE
C
C
      IF(B.NE.0) THEN
C
       IF (NFACT(N-1).GT.MAXEP) THEN
        WRITE(15,200)
  200   FORMAT(' N TOO LARGE FOR DIM OF IBRPM IN SUB EPSB ')
        STOP
       ELSE
       ENDIF
C
       IDBL=0
C      **** WANT GPRBIG TO DOUBLE COUNT ****
       CALL GPRBIG(NMIN,NMAX,N-1,INDX,NPRMS,IBRPM,IDBL)
C      **** ASSIGN SIGN TO EACH PERM  ****
       DO 100 IS1=1,NPRMS
        DO 110 IS2=1,N-1
         IPRM(IS2)=IBRPM(IS1,IS2)
  110   CONTINUE
        CALL PSIGN(N-1,IPRM,ISGNPM)
        IBRPM(IS1,N)=ISGNPM
  100  CONTINUE
C
      ELSE
      ENDIF
C
      IF(ICLDI.NE.0) THEN
C     **** WRITE IBRPM ****
       WRITE(15,50)NPRMS
   50  FORMAT(' NPRMS  = ',I6)
C
       WRITE(15,30)
   30  FORMAT(' IBRPM:')
       DO 60 J=1,NPRMS
        WRITE(15,20) (IBRPM(J,I),I=1,N)
   20   FORMAT(25I6)
   60  CONTINUE
      ELSE
      ENDIF
C
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE WFST
C     **** PRINT STATE INFO TO BE USED IN      ****
C     **** WF FORTRAN (OR WFBIG FORTRAN)       ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/STATE/MSTATE(27608,25),MSTINF(6902,8)
      COMMON/NST/NUMSTA
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
c     INTEGER NSTATB(6902,30),NSTATD(6902,30)
      INTEGER NSTATB(6902,03,15),NSTATD(6902,03,15)
      INTEGER locb(03),locd(03),lb(03),ld(03) 
      INTEGER N,NF,B,K
C
      MXLN=15
C     **** MAX LENGTH FOR STATE DATA ****
C
      DO 5 I=1,NUMSTA
       DO 7 J=1,MXLN
       DO 7 jf=1,nf
        NSTATB(I,jf,J)=0
        NSTATD(I,jf,J)=0
    7  CONTINUE
    5 CONTINUE
C
C
      DO 10 NS=1,NUMSTA
C
       LOC = MSTINF( NS, 1)
       LNG = MSTINF( NS, 2)
       NMES = MSTINF( NS, 3)
       NBAR = MSTINF( NS, 4)
       NB=N*NBAR + NMES
       ND=NMES
C
c      NSTATB(NS,1)=NB
c      NSTATD(NS,1)=-ND
C
       do 18 jfl=1,nf
        LOCB(jfl)=0
        LOCD(jfl)=0
   18  continue
       DO 20 LS=1,LNG
        ifl=mstate(loc+3,ls)
c       ****flavor ****
        IF(MSTATE(LOC,LS).EQ.1) THEN
C        *** QUARK ****
         LOCB(ifl)=LOCB(ifl)+1
         NSTATB(NS,ifl,LOCB(ifl)+1)=MSTATE(LOC+2,LS)
         NSTATB(NS,ifl,1)=locb(ifl)          
c        ****momentum ****
        ELSE
         IF(MSTATE(LOC+1,LS).EQ.1) THEN
C         *** ANTIQUARK ****
          LOCD(ifl)=LOCD(ifl)+1
          NSTATD(NS,ifl,LOCD(ifl)+1)=MSTATE(LOC+2,LS)
          NSTATD(NS,ifl,1)=-locd(ifl)         
c         ****momentum ****
         ELSE
         ENDIF
        ENDIF
   20  CONTINUE
C
   10 CONTINUE
C
C
      WRITE(15,30)
   30 FORMAT(' STATE INFO: QUARKS ')
      DO 50 JQ=1,NUMSTA
c      LB=NSTATB(JQ,1)
       do 52 jqf=1,nf
        lb(jqf)=nstatb(jq,jqf,1)
   52  continue
c      WRITE(15,60)(NSTATB(JQ,L),L=1,LB+1)
       WRITE(15,60)((NSTATB(JQ,kf,L),L=1,lb(kf)+1),kf=1,nf)
   60  FORMAT(25I6)
   50 CONTINUE
C
C
      WRITE(15,40)
   40 FORMAT(' STATE INFO: ANTIQUARKS ')
      DO 70 JA=1,NUMSTA
c      LD=ABS(NSTATD(JA,1))
       do 62 jqf=1,nf
        ld(jqf)=abs(nstatd(ja,jqf,1))
   62  continue
c      WRITE(15,80)(NSTATD(JA,M),M=1,LD+1)
       WRITE(15,80)((NSTATD(JA,kf,M),M=1,ld(kf)+1),kf=1,nf)
   80  FORMAT(25I6)
   70 CONTINUE
C
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE PRZ(Z)
C     **** PRINT INFO ABOUT MATRIX WHICH ****
C     **** ORTHONORMALIZES STATES        ****
C     **** FOR USE IN WF OR WFBIG FORTRAN****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NST/NUMSTA
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      REAL*8 Z(6902,6902)
      INTEGER IZCNT(6902,6902)
      INTEGER N,NF,B,K
C
      MXWD=6902
      EPS=1.D-12
C
      DO 5 JZ=1,NUMSTA
       DO 7 KZ=1,MXWD
        IZCNT(JZ,KZ)=0
    7  CONTINUE
    5 CONTINUE
C
      DO 10 J=1,NUMSTA
       NSX=0
       DO 20 I=1,NUMSTA
        IF(DABS(Z(I,J)).GT.EPS) THEN
         NSX=NSX+1
         IZCNT(J,NSX+1)=I
        ELSE
        ENDIF
        IZCNT(J,1)=NSX
   20  CONTINUE
   10 CONTINUE
C
      WRITE(15,30)
   30 FORMAT(' Z INFO ')
      DO 40 LP=1,NUMSTA
       ISX=IZCNT(LP,1)
       WRITE(15,50)ISX
   50  FORMAT(I7)
       WRITE(15,60)(IZCNT(LP,MM),MM=2,ISX+1)
   60  FORMAT(15I7)
   40 CONTINUE
C
      RETURN
      END
C
C
C
C
C
      SUBROUTINE PAVILL(JSTATE,LNGSTA,LCTCHK)
C     **** CHECK IF STATE SATIFIES PAULI-****
C     **** VILLARS CUT-OFF CONDITION;    ****
C     **** THEN LCTCHK =1; ELSE ZERO     ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      COMMON/CTOF/CUTOFF
      common/masses/rmq(03), iflv(03)
      INTEGER JSTATE(4,25)
      INTEGER N,NF,B,K
      real*8 rm
C
      EPS=1.D-6
C
      IF(CUTOFF.LE.0) THEN
C      **** NO CUTOFF CONDITION ****
       LCTCHK=1
      ELSE
       ELC=0
       DO 10 IC=1,LNGSTA
        X=DFLOAT(JSTATE(3,IC))/DFLOAT(K)
        IF(X.LT.EPS) THEN
         LCTCHK=0
         RETURN
        ELSE
C        ELC=ELC+ 1.D0/X
         rm = rmq(jstate(4,ic))
         ELC=ELC+ rm*rm/X
        ENDIF
   10  CONTINUE
       IF(ELC.LT.CUTOFF) THEN
        LCTCHK=1
       ELSE
        LCTCHK=0
       ENDIF
      ENDIF
C
      RETURN
      END
c
C
C
C
      SUBROUTINE flvchk(JSTATE,LNGSTA,ifchk)
C     **** check if state has correct flavor qm numbers; ****
C     **** if so, ifchk =1; else 0                       ****
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PARAM/RLAMB,N,NF,B,K,MASS
      common/masses/rmq(03), iflv(03)
      INTEGER JSTATE(4,25)
      INTEGER N,NF,B,K
      integer fct(03)
C
      ifchk=1
      do 5 id=1,nf
       fct(id) = 0
    5 continue
c
      DO 10 IC=1,LNGSTA
       mf = jstate(4,ic)
       if (jstate(1,ic).eq.1) then
        fct(mf) = fct(mf) + 1
       else
        fct(mf) = fct(mf) - 1
       endif
   10 CONTINUE
c
      do 20 id=1,nf
       if ( fct(id).ne.iflv(id) ) ifchk = 0
c      WRITE(*,15) id, fct(id), iflv(id)
c  15  FORMAT('i: ',i3,'  fct(i): ', i3, '  iflv(i)', i3)
   20 continue
C
      RETURN
      END
c
c
c
      SUBROUTINE TRR8(A,N,NP,D,E)
c     **** This is num recipe prog. tred2 with double prec. ****
      IMPLICIT REAL*8 (A-H,O-Z)
c     DIMENSION A(NP,NP),D(NP),E(NP)
      REAL*8    A(NP,NP),D(NP),E(NP)
      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.d0
          SCALE=0.d0
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.d0)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2.d0
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.d0
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.d0
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.d0
      E(1)=0.d0
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.d0)THEN
          DO 21 J=1,L
            G=0.d0
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.d0
            A(J,I)=0.d0
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE TQR8(D,E,N,NP,Z)
c     **** This is num. recipe program tqli with double prec. ****
      IMPLICIT REAL*8 (A-H,O-Z)
c     DIMENSION D(NP),E(NP),Z(NP,NP)
      REAL*8    D(NP),E(NP),Z(NP,NP)
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.d0
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
c           IF(ITER.EQ.30)PAUSE 'too many iterations'
            IF(ITER.EQ.1000)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.d0*E(L))
            R=SQRT(G**2.d0+1.d0)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.d0
            C=1.d0
            P=0.d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2.d0+1.d0)
                E(I+1)=F*R
                S=1.d0/R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2.d0+1.d0)
                E(I+1)=G*R
                C=1.d0/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.d0
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
c
c
c
      SUBROUTINE ESRTR8(D,V,N,NP)
c     **** eigenvalue sorter ****
c     **** modified for double precision ****
      implicit real*8 (a-h,o-z)
      real*8 D(NP),V(NP,NP)
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
c         IF(D(J).GE.P)THEN
          IF(D(J).lt.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END
