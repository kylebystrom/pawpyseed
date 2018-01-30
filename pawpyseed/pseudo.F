#include "symbol.inc"

MODULE PSEUDO_struct
  USE prec
  USE radial
  !
  ! pseudo potential discription include file
  ! only included if MODULES are not supported
  !only P
  INTEGER, PARAMETER:: NEKERR=100
  INTEGER, PARAMETER:: NPSPTS=1000
  INTEGER, PARAMETER:: NPSNL=100,NPSRNL=100
  !
  ! structure to support radial grid
  !
  TYPE potcar
     REAL(q) ZVALF,POMASS,RWIGS  ! valence, mass, wigner seitz radius
     REAL(q) ZVALF_ORIG          ! original valence
     ! ZVAL might be reset when core level shifts are calculated
     REAL(q) EATOM               ! atomic energy
     REAL(q) EATOM_CORRECTED     ! atomic energy with convergence correction
     REAL(q) EGGA(4)             ! atomic energy for GGA
     REAL(q) ENMAXA,ENMINA,EAUG  ! energy cutoffs
     REAL(q) QOPT1,QOPT2         ! real space optimized for this cutoffs
     REAL(q) EKECUT(NEKERR),EKEERR(NEKERR) ! description of cutofferrors
     INTEGER LEXCH             ! exchange type
     REAL(q) PSGMAX              ! maximal G for local potential
     REAL(q) PSMAXN              ! maximal G for non local potential
     REAL(q) PSRMAX              ! maximal r for non local contrib.  (in fact rmax=PSRMAX/NPSNL*(NPSNL-1))
     REAL(q) PSDMAX              ! maximal r for augmentation charge (in fact rmax=PSDMAX/NPSNL*(NPSNL-1))
     REAL(q) RDEP                ! outermost grid point on radial grid used in LDA+U
                                 ! usually equivalent to largest matching radius in PSCTR
     REAL(q) RCUTRHO             ! cutoff radius for pseudo charge dens.
     REAL(q) RCUTATO             ! cutoff radius for atomic charge dens.
     REAL(q),POINTER :: RHOSPL(:,:) => NULL() ! pseudopot. charge dens. in real space spline coeff.
     REAL(q),POINTER :: ATOSPL(:,:) => NULL() ! atomic charge dens. in real space spline coeff.
     REAL(q),POINTER :: USESPL(:,:) => NULL() ! atomic charge dens. in real space spline coeff.
     REAL(q) USEZ
     REAL(q) USECUT
     REAL(q),POINTER :: PSP(:,:) ! local pseudopotential in rec. space
     REAL(q),POINTER :: PSPCOR(:)! partial core information in rec. space
     REAL(q),POINTER :: PSPTAU(:)! partial kinetic energy density  information in rec. space
     REAL(q),POINTER :: PSPTAUVAL(:)! kinetic energy density of valence electrons information in rec. space
     REAL(q),POINTER :: PSPRHO(:)! atomic pseudo charge density in rec. space
     REAL(q),POINTER :: PSPNL(:,:)  ! non local proj. rec. space
     REAL(q),POINTER :: PSPNL_SPLINE(:,:,:) => NULL() ! non local proj. rec. space spline fit
     REAL(q),POINTER :: PSPRNL(:,:,:) ! non local proj. real space
     REAL(q),POINTER :: DION(:,:)! non local strength
     REAL(q),POINTER :: QION(:,:)! spherical augmentation charge
     REAL(q),POINTER :: QTOT(:,:)! total charge in each channel
     REAL(q),POINTER :: QPAW(:,:,:)! integrated augmentation charge 
             ! stores essentially WAE*WAE-WPS*WPS r^l
     COMPLEX(q),POINTER :: JPAW(:,:,:,:)! 
             ! stores essentially <WAE| j(r) Y_lm(r) | WAE> -<WPS| j(r) Y_lm(r) | WPS>
     REAL(q),POINTER :: QPAW_FOCK(:,:,:,:)
             ! similar to QPAW but stores coefficients that determine how much of 
             ! AUG_FOCK is added
     REAL(q),POINTER :: QATO(:,:)! initial occupancies (in atom)
     REAL(q),POINTER :: QDEP(:,:,:) ! L-dependent augmentation charges on regular grid
     REAL(q),POINTER :: QDEP_FOCK(:,:,:) ! L-dependent charge with zero moment on regular grid
                                 ! for technical reasons this also includes all elements in QDEP
                                 ! see comments in fast_aug.F
     REAL(q),POINTER :: NABLA(:,:,:)! atomic augmentation matrix elements of nabla operator
     REAL(q) PSCORE              ! is equal to V(q) + 4 Pi / q^2
     REAL(q) ESELF               ! self energy of pseudized ion (usually not used)
     INTEGER,POINTER :: LPS(:)   ! L quantum number for each PP
     INTEGER,POINTER :: NLPRO(:) ! unused
     INTEGER,POINTER :: NDEP(:,:)! number of augmentation channels per ll'
     REAL(q),POINTER :: E(:)     ! linearization energies
     REAL(q),POINTER :: E_ORIG(:)! linearization energies as stored on POTCAR
     ! quantities defined on radial grid
     TYPE (rgrid)    :: R        ! radial grid
     REAL(q),POINTER :: POTAE(:) ! frozen core AE potential on r-grid for atomic reference configuration (valence only) as read from file (must not be updated to get correct core states)
     REAL(q),POINTER :: POTAE_XCUPDATED(:)! as above, might be updated if xc potential changes
     REAL(q),POINTER :: POTPS(:) ! local PP on r-grid for atomic reference configuration (valence only)
     REAL(q),POINTER :: POTPSC(:)! local PP on r-grid (core only), V_H[\tilde{n}_Zc]
     REAL(q),POINTER :: KLIC(:)  ! averaged local exchange potential (KLI) on r-grid
     REAL(q),POINTER :: RHOAE(:) ! frozen core charge rho(r)r^2 on r-grid
     REAL(q),POINTER :: TAUAE(:) ! kinetic energy density of core
     REAL(q),POINTER :: RHOPS(:) ! frozen pseudo partial core charge rho(r)r^2 on r-grid
     REAL(q),POINTER :: TAUPS(:) ! kinetic energy density of partial core
     REAL(q),POINTER :: WAE(:,:) ! ae valence wavefunction on r-grid
     REAL(q),POINTER :: WPS(:,:) ! pseudo valence wavefunction on r-grid
     REAL(q),POINTER :: AUG(:,:) ! L-dependent augmentation charge on r-grid
     REAL(q),POINTER :: AUG_SOFT(:,:)  ! L-dependent augmentation charge on r-grid (soft version)
     REAL(q),POINTER :: AUG_FOCK(:,:,:)! L-dependent charge with zero moment on r-grid
     REAL(q)         :: DEXCCORE ! exchange correlation energy of frozen core
     REAL(q)         :: DEXCCOREM! exchange correlation energy of frozen core (MetaGGA)
     ! relaxed core stuff
     REAL(q),POINTER :: C(:,:)     ! partial wave expansion coefficients
     REAL(q),POINTER :: BETA(:,:)  ! projectors on radial grid
     REAL(q),POINTER :: RHOAE00(:) ! spherical component of valence AE charge density
     REAL(q),POINTER :: RHOPS00(:) ! spherical component of valence PS charge density
     REAL(q),POINTER :: RHOPSPW(:) ! spherical component of valence PS charge density from PW grid
     REAL(q),POINTER :: V00_AE(:)  ! spherical component of AE potential
     REAL(q),POINTER :: V00_PS(:)  ! spherical component of PS potential
     REAL(q),POINTER :: WKINAE(:,:)! T | \psi_i >
     REAL(q),POINTER :: WKINPS(:,:)! T | \tilde{\psi}_i >
     REAL(q),POINTER :: DIJ(:,:)   ! < \psi_i | T+V |\psi_j > - < \tilde{\psi}_i | T+\tilde{V} |\tilde{\psi}_j >
     REAL(q),POINTER :: QIJ(:,:)   ! < \psi_i | \psi_j > - < \tilde{\psi}_i | \tilde{\psi}_j >
     REAL(q),POINTER :: POTAEC(:)  ! AE core potential 
     REAL(q)         :: AVERAGEPOT(3) ! average local potential (AE,PS,PW)
     REAL(q)         :: VPSRMAX    ! local PP at boundary
     REAL(q),POINTER :: CLEV(:)    ! core state eigenenergies
     REAL(q)         :: ECORE(2)   ! core contributions to total energy
     REAL(q)         :: VCA        ! weight of this potential (can be overwritten by INCAR VCA tag)
     ! atomic stuff
     REAL(q),POINTER :: ATOMIC_J(:)   ! relativistic quantum number
     REAL(q),POINTER :: ATOMIC_E(:)   ! eigenenergy
     REAL(q),POINTER :: ATOMIC_OCC(:) ! occupation number
     INTEGER,POINTER :: ATOMIC_N(:)   ! main quantum number
     INTEGER,POINTER :: ATOMIC_L(:)   ! angular momentum
     ! most integers are on the end (alignment ?)
     INTEGER      LDIM           ! (leading) dimension for l channels (>=LMAX)
     INTEGER      LMAX           ! total number of l-channels for non local PP
     INTEGER      LMDIM          ! (leading) dimension for lm channels (>=LMMAX)
     INTEGER      LMMAX          ! total number nlm-channels for non local PP
     INTEGER      LDIM2          ! dimension for augmentation arryas
     LOGICAL      LREAL          ! real space optimized ?
     LOGICAL      LUNSCR         ! partial core has been unscreened
     INTEGER      LMAX_CALC      ! maximum L for onsite terms in PAW
     ! might be overwritten by LMAXPAW line in the INCAR file
     REAL(q)      ZCORE          ! charge of core
     CHARACTER*40 SZNAMP         ! header
     CHARACTER*2  ELEMENT        ! Name of element
  END TYPE potcar
END MODULE PSEUDO_struct

MODULE PSEUDO
      USE PSEUDO_struct
!**********************************************************************
!  interface to a function that returns a pointer
!  to the current pseudopotential, the function itself
!  can be found in core_rel.F 
!**********************************************************************
      INTERFACE
      FUNCTION PP_POINTER(P, NI, NT)
        USE pseudo_struct
        IMPLICIT NONE
        TYPE (potcar),TARGET :: P(:)
        INTEGER :: NI      ! ion index
        INTEGER :: NT      ! type index
        TYPE (potcar), POINTER :: PP_POINTER
      END FUNCTION PP_POINTER
      END INTERFACE      

      CONTAINS

!**************** SUBROUTINE RD_PSEUDO *********************************
! RCS:  $Id: pseudo.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
!  reads in all pseudopotential from the POTCAR file
!
!  check if LDIM and LMDIM is sufficient
!  if not tell the user to increase those numbers
!
!***********************************************************************

      SUBROUTINE RD_PSEUDO(INFO,P, &
     &           NTYP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           POMASS,RWIGS,TYPE,VCA, &
     &           IU0,IU6,NWRITE,LPAW)
      USE base
      USE ini
      USE main_mpi
      IMPLICIT NONE

      TYPE(INFO_STRUCT) INFO
      INTEGER  NTYP,NTYPD      ! number of types (acutal/dimension)
      INTEGER LDIM,LDIM2,LMDIM ! see pseudo.inc
      TYPE (potcar) P(NTYPD)   ! PP information
      LOGICAL LPAW
      INTEGER IU6,IU0          ! where I/O goes
! wigner seitz radius and mass found on INCAR
      REAL(q)    RWIGS(NTYPD),POMASS(NTYPD),VCA(NTYPD)
      REAL(q)    BOUNDARY
      CHARACTER(LEN=2) :: TYPE(NTYPD)
! dynamical work space
      INTEGER,PARAMETER :: ISDIM=100
      CHARACTER (80) STRING(ISDIM)
      CHARACTER (40) FORM
! temporary varibales
      INTEGER IDUMM1,IDUMM2,IDUMM3,L2,NREAD,NL,IC,L,LP,I,J,K,CHANNELS,LMAX
      LOGICAL LDUM
      CHARACTER (1)  CSEL
      INTEGER IREAD,IWRITE,NWRITE,NMAX
      LOGICAL LPOTPSC
      INTEGER IERR
#ifdef debug
      REAL(q) :: QR, F, DUMMY
#endif

      OPEN(UNIT=10,FILE=DIR_APP(1:DIR_LEN)//'POTCAR',STATUS='OLD',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=10,FILE='POTCAR',STATUS='OLD')
      ENDIF

      LPAW = .FALSE.
      INFO%LOVERL=.FALSE.
      INFO%LCORE =.FALSE.
!-----------------------------------------------------------------------
! loop over all pseudpotentials on POTCAR file
!-----------------------------------------------------------------------
      NTYP=1

      forever: DO NTYP=1,NTYPD
      P(NTYP)%LDIM =LDIM
      P(NTYP)%LMDIM=LMDIM
      P(NTYP)%LDIM2=LDIM2
      P(NTYP)%ESELF=0
      NULLIFY(P(NTYP)%QPAW)
      NULLIFY(P(NTYP)%JPAW)
      NULLIFY(P(NTYP)%QPAW_FOCK)
      NULLIFY(P(NTYP)%QTOT)
      NULLIFY(P(NTYP)%QATO)
      NULLIFY(P(NTYP)%NABLA)
      NULLIFY(P(NTYP)%ATOMIC_N)
      
      READ(10,'(A40)',END=100,ERR=100) P(NTYP)%SZNAMP

      IF (IU6>=0) &
      WRITE(IU6,*)'POTCAR:  ',P(NTYP)%SZNAMP

      READ(10,*) P(NTYP)%ZVALF
      P(NTYP)%ZVALF_ORIG=P(NTYP)%ZVALF
!-----------------------------------------------------------------------
!  if there is PSCTR header read in and parse information from this
!  header
!-----------------------------------------------------------------------
      READ(10,'(1X,A1)') CSEL
      IF (CSEL /= 'p') THEN
        IF (IU6>=0) &
        WRITE(IU6,*)'this version requires full pseudpotential ', &
     &              ' generation information'
        IF (IU0>=0) &
        WRITE(IU0,*)'this version requires full pseudpotential ', &
     &              ' generation information'
        STOP
      ENDIF
! read contents to STRING
        CALL RDPARA(10,ISDIM,STRING,IREAD,IWRITE)

        IF (NWRITE>=2.AND. IU6>=0 ) &
           WRITE(IU6,'(A80)') (STRING(I),I=1,IWRITE)

! parse from STRING necessary information
        CALL RDPARS(ISDIM,STRING,IREAD,P(NTYP),IU6)
        IF (TYPE(NTYP)=='  ') THEN
           TYPE(NTYP)=P(NTYP)%ELEMENT
        ENDIF
        IF (TYPE(NTYP)/=P(NTYP)%ELEMENT.AND.NWRITE>=0.AND. IU0>=0) THEN
          WRITE(IU0,*)'WARNING: type information on POSCAR and POTCAR are incompatible'
          WRITE(IU0,*)'POTCAR overwrites the type information in POSCAR'
          WRITE(IU0,'(A,I3,A,2A3)')' typ ',NTYP,' type information: ',TYPE(NTYP),P(NTYP)%ELEMENT
        ENDIF
! set and check POMASS
        IF (POMASS(NTYP)<0) POMASS(NTYP)=P(NTYP)%POMASS
        IF (ABS(POMASS(NTYP)-P(NTYP)%POMASS)>1E-3_q) THEN
        IF (NWRITE>=0.AND. IU0>=0) THEN
          WRITE(IU0,*)'WARNING: mass on POTCAR and INCAR are incompatible'
          WRITE(IU0,*)' typ',NTYP,' Mass',POMASS(NTYP),P(NTYP)%POMASS
        ENDIF
        ENDIF
! set VCA tag
        IF (VCA(NTYP)<0) VCA(NTYP)=P(NTYP)%VCA
! set RWIGS
        IF (RWIGS(NTYP)==0) RWIGS(NTYP)=P(NTYP)%RWIGS

!-----------------------------------------------------------------------
! local potential, gradient corrections and type of gradient corrections
! partial core and atomic charge density
!-----------------------------------------------------------------------
      READ(10,'(1X,A1)') CSEL
      READ(10,*) P(NTYP)%PSGMAX
      ALLOCATE(P(NTYP)%PSP(NPSPTS,5))

      READ(10,*) (P(NTYP)%PSP (I,2),I=1,NPSPTS)
      DO I=1,NPSPTS
          P(NTYP)%PSP(I,1)=(P(NTYP)%PSGMAX/NPSPTS)*(I-1)
      ENDDO
! PSCORE is equal to V(q) + 4 Pi / q^2
      P(NTYP)%PSCORE=P(NTYP)%PSP(1,2)

      CALL SPLCOF(P(NTYP)%PSP(1,1) ,NPSPTS,NPSPTS,0._q)
      IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' local pseudopotential read in'

!  gradient corrections this is now overwritten by PSCTR header
      READ(10,'(1X,A1)') CSEL
      IF (CSEL=='g') THEN
         READ(10,*)
         READ(10,'(1X,A1)') CSEL
      ENDIF
!
!  partial core and atomic chargedensity
!
      IF (CSEL=='c') THEN
        INFO%LCORE=.TRUE.
        ALLOCATE(P(NTYP)%PSPCOR(NPSPTS))
        IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' partial core-charges read in'

        READ(10,*) (P(NTYP)%PSPCOR (I),I=1,NPSPTS)
        READ(10,'(1X,A1)') CSEL
      ELSE
         NULLIFY(P(NTYP)%PSPCOR)
      ENDIF
      IF (CSEL=='k') THEN
        ALLOCATE(P(NTYP)%PSPTAU(NPSPTS))
        IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' partial kinetic energy density read in'

        READ(10,*) (P(NTYP)%PSPTAU (I),I=1,NPSPTS)
        READ(10,'(1X,A1)') CSEL
      ELSE
         NULLIFY(P(NTYP)%PSPTAU)
      ENDIF
      IF (CSEL=='K') THEN
        ALLOCATE(P(NTYP)%PSPTAUVAL(NPSPTS))
        IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' kinetic energy density of atom read in'

        READ(10,*) (P(NTYP)%PSPTAUVAL(I),I=1,NPSPTS)
        READ(10,'(1X,A1)') CSEL
      ELSE
         NULLIFY(P(NTYP)%PSPTAUVAL)
      ENDIF

      ALLOCATE(P(NTYP)%PSPRHO(NPSPTS))
      READ(10,*) (P(NTYP)%PSPRHO (I),I=1,NPSPTS)
      IF (NWRITE>=0 .AND. IU6>=0) WRITE(IU6,*)' atomic valenz-charges read in'

      CHANNELS=0
      P(NTYP)%LMMAX =0
      P(NTYP)%LMAX  =0
      P(NTYP)%PSMAXN=0
!-----------------------------------------------------------------------
! depletion charges to 0
!-----------------------------------------------------------------------
      ALLOCATE(P(NTYP)%DION(LDIM,LDIM),P(NTYP)%QION(LDIM,LDIM), &
     &         P(NTYP)%LPS(LDIM),P(NTYP)%NLPRO(LDIM), &
     &         P(NTYP)%PSPNL(0:NPSNL,LDIM),P(NTYP)%PSPRNL(NPSRNL,5,LDIM),P(NTYP)%E(LDIM))
      IF (INFO%NLSPLINE) THEN
         ALLOCATE(P(NTYP)%PSPNL_SPLINE(0:NPSNL,5,LDIM))
      ELSE
         NULLIFY (P(NTYP)%PSPNL_SPLINE)
      ENDIF

! Intel efc compiler workaround (at least version 6.X)
!      P(NTYP)%DION=0
!      P(NTYP)%QION=0
      DO I=1,LDIM
         DO J=1,LDIM
            P(NTYP)%DION(I,J)=0
            P(NTYP)%QION(I,J)=0
         ENDDO
      ENDDO
      P(NTYP)%LPS=0
      P(NTYP)%NLPRO=0

      P(NTYP)%PSDMAX=0
      NULLIFY(P(NTYP)%QDEP,P(NTYP)%AUG,P(NTYP)%AUG_SOFT)
      NULLIFY(P(NTYP)%NDEP)
!-----------------------------------------------------------------------
!  read reciprocal projection operators
!-----------------------------------------------------------------------
      READ(10,*,ERR=620,END=600)  P(NTYP)%PSMAXN,LDUM
  620 IF (P(NTYP)%PSMAXN==0) GOTO 600

  610 READ(10,'(1X,A1)',ERR=600,END=600) CSEL
      IF (CSEL== 'D'  .OR. CSEL== 'A' .OR. CSEL=='P' ) GOTO 650
      IF (CSEL== 'E' ) GOTO 600

      READ(10,*)  P(NTYP)%LPS(CHANNELS+1),P(NTYP)%NLPRO(CHANNELS+1), &
     &            P(NTYP)%PSRMAX

      IC= P(NTYP)%NLPRO(CHANNELS+1)
      IF (CHANNELS+IC>LDIM) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'RD_PSEUDO: internal ERROR: increase LDIM'
         STOP
      ENDIF

      DO 615 L=CHANNELS+1,CHANNELS+IC
        P(NTYP)%LPS(L)   =P(NTYP)%LPS   (CHANNELS+1)
        P(NTYP)%NLPRO(L) =P(NTYP)%NLPRO (CHANNELS+1)
  615 CONTINUE

!-----Multipliers DION
      READ(10,*) &
     & ((P(NTYP)%DION(L,LP),L=CHANNELS+1,CHANNELS+IC),LP=CHANNELS+1,CHANNELS+IC)

      P(NTYP)%LMMAX=P(NTYP)%LMMAX+(2*P(NTYP)%LPS(CHANNELS+1)+1)*IC
      IF (P(NTYP)%LMMAX>LMDIM) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'MAIN: ERROR: increase LMDIM to ',P(NTYP)%LMMAX
         STOP
      ENDIF

      DO 630 NL=1,IC
!-----reciprocal projection operator
        READ(10,*)
        READ(10,*)(P(NTYP)%PSPNL  (I,CHANNELS+NL),I=1,NPSNL)

        IF (NWRITE>=0.AND.IU6>=0) &
             WRITE(IU6,*)' non local Contribution for L=', &
                       P(NTYP)%LPS(CHANNELS+1),' read in'
        IF (MOD(P(NTYP)%LPS(CHANNELS+1),2)==0) THEN
           P(NTYP)%PSPNL(0,CHANNELS+NL) =  P(NTYP)%PSPNL(2,CHANNELS+NL)
        ELSE
           P(NTYP)%PSPNL(0,CHANNELS+NL) = -P(NTYP)%PSPNL(2,CHANNELS+NL)
        ENDIF
        IF (INFO%NLSPLINE) THEN
           P(NTYP)%PSPNL_SPLINE(:,2,CHANNELS+NL)=P(NTYP)%PSPNL(:,CHANNELS+NL)
           DO I=0,NPSNL
              P(NTYP)%PSPNL_SPLINE (I,1,CHANNELS+NL)=(P(NTYP)%PSMAXN/NPSNL)*(I-1)
           ENDDO
           ! possibly this works, but one would need to inspect visually the behaviour at small q
          ! BOUNDARY=(P(NTYP)%PSPNL (2,2,CHANNELS+NL)-P(NTYP)%PSPNL (1,2,CHANNELS+NL))/(P(NTYP)%PSMAXN/NPSNL)
          ! IF (P(NTYP)%LPS(CHANNELS+1)/=1) THEN
          !    BOUNDARY=0
          ! ENDIF
           BOUNDARY=1E30
           CALL SPLCOF(P(NTYP)%PSPNL_SPLINE (0,1,CHANNELS+NL) ,NPSNL+1,NPSNL+1,BOUNDARY)
        ENDIF
!-----real space projection operator
        READ(10,*)
        READ(10,*)(P(NTYP)%PSPRNL (I,2,CHANNELS+NL),I=1,NPSRNL)

        IF (NWRITE>=0.AND.IU6>0) &
               WRITE(IU6,*)'   real space projection operators read in'

        DO I=1,NPSRNL
           P(NTYP)%PSPRNL (I,1,CHANNELS+NL)=(P(NTYP)%PSRMAX/NPSRNL)*(I-1)
        ENDDO

        BOUNDARY=(P(NTYP)%PSPRNL (2,2,CHANNELS+NL)-P(NTYP)%PSPRNL (1,2,CHANNELS+NL))/(P(NTYP)%PSRMAX/NPSRNL)
        IF (P(NTYP)%LPS(CHANNELS+1)/=1) THEN
           BOUNDARY=0
        ENDIF
        CALL SPLCOF(P(NTYP)%PSPRNL (1,1,CHANNELS+NL) ,NPSRNL,NPSRNL,BOUNDARY)
#ifdef debug
        DO K=1,(NPSRNL-1)*10+1
           QR=(P(NTYP)%PSRMAX/NPSRNL)/10*(K-1)
           CALL SPLVAL(QR,F,DUMMY,P(NTYP)%PSPRNL(1,1,CHANNELS+NL),NPSRNL,NPSRNL)
           WRITE(78,'(2F14.7)') QR,F
        ENDDO
        WRITE(78,*)
#endif


  630 CONTINUE

      CHANNELS=CHANNELS+IC
      P(NTYP)%LMAX=CHANNELS
      GOTO 610

!-----------------------------------------------------------------------
!  read depletion charges
!  old dataset start with 'D' new one with tag 'A'
!-----------------------------------------------------------------------
  650 CONTINUE
   aug: IF ( CSEL /= 'P') THEN

      INFO%LOVERL=.TRUE.
      ALLOCATE(P(NTYP)%QDEP(NPSRNL,5,LDIM2),P(NTYP)%NDEP(LDIM,LDIM))

      IF (NWRITE>=0.AND.IU6>=0) &
           WRITE(IU6,*)'   augmentation charges read in'

      L2=1
      DO L =1 ,CHANNELS
      DO LP=L ,CHANNELS
      NREAD=0

  661   CONTINUE

        IF (L2> LDIM2) THEN
           IF (IU0>=0) &
           WRITE(IU0,*)'RD_PSEUDO: internal ERROR increase LDIM2'
           STOP
        ENDIF

        IF (CSEL=='D') THEN
         READ(10,*) IDUMM1,IDUMM2,P(NTYP)%QION(L,LP),P(NTYP)%PSDMAX
         P(NTYP)%NDEP(L,LP)=-1
        ELSE
         READ(10,*) IDUMM1,IDUMM2,IDUMM3,P(NTYP)%NDEP(L,LP), &
     &              P(NTYP)%QION(L,LP),P(NTYP)%PSDMAX
         IF (P(NTYP)%LPS(L)/=IDUMM1 .OR. P(NTYP)%LPS(LP)/=IDUMM2) THEN
            IF (IU0>=0) &
            WRITE(IU0,*)'ERROR: on reading POTCAR augmentation charges', &
            '       are wrong'
            STOP
         ENDIF
        ENDIF
        P(NTYP)%QION(LP,L) = P(NTYP)%QION(L,LP)

        IF (P(NTYP)%NDEP(L,LP)==0) GOTO 667

        READ(10,*) (P(NTYP)%QDEP(I,2,L2),I=1,NPSRNL)

        DO 665 I=1,NPSRNL
         P(NTYP)%QDEP(I,1,L2)=(I-1)*P(NTYP)%PSDMAX/NPSRNL
  665   CONTINUE
        ! changed  12Dec96 (gK)
        P(NTYP)%QDEP(NPSRNL,2,L2)=0 ! fix last value to zero

        CALL SPLCOF(P(NTYP)%QDEP(1,1,L2),NPSRNL,NPSRNL,0._q)

        L2   =L2   +1
        NREAD=NREAD+1

  667   READ(10,*,ERR=600,END=600)
        IF (NREAD<ABS(P(NTYP)%NDEP(L,LP))) GOTO 661

      ENDDO
      ENDDO
      ELSE aug
!-----------------------------------------------------------------------
!  read paw data sets
!-----------------------------------------------------------------------
         INFO%LOVERL=.TRUE.
         LPAW=.TRUE.
         READ(10,*) NMAX,P(NTYP)%PSDMAX

         LMAX=0
         DO I=1,CHANNELS
            LMAX=MAX( P(NTYP)%LPS(I),LMAX )
         ENDDO
         LMAX=LMAX*2
         P(NTYP)%LMAX_CALC=LMAX

         READ(10,'(A)') FORM            ! format of remaining entities
         ALLOCATE( P(NTYP)%QPAW(CHANNELS,CHANNELS,0:LMAX),P(NTYP)%QATO(CHANNELS,CHANNELS), &
                   P(NTYP)%R%R(NMAX), P(NTYP)%POTAE(NMAX),P(NTYP)%POTAE_XCUPDATED(NMAX), P(NTYP)%POTPS(NMAX),  &
                   P(NTYP)%POTPSC(NMAX),P(NTYP)%TAUAE(NMAX),P(NTYP)%TAUPS(NMAX), &
                   P(NTYP)%RHOAE(NMAX), P(NTYP)%RHOPS(NMAX),P(NTYP)%QTOT(CHANNELS,CHANNELS), &
                   P(NTYP)%WAE(NMAX,CHANNELS), P(NTYP)%WPS(NMAX,CHANNELS), &
                   P(NTYP)%NABLA(3,LMDIM,LMDIM))

! Intel efc compiler workaround (at least version 6.X)
!         P(NTYP)%NABLA=0
         DO I=1,LMDIM
            DO J=1,LMDIM
               P(NTYP)%NABLA(1,I,J)=0
               P(NTYP)%NABLA(2,I,J)=0
               P(NTYP)%NABLA(3,I,J)=0
            ENDDO
         ENDDO

         READ(10,*)
! Intel efc compiler workaround (at least version 6.X)
!         P(NTYP)%QPAW=0
         DO I=1,CHANNELS
            DO J=1,CHANNELS
               DO K=0,LMAX
                  P(NTYP)%QPAW(I,J,K)=0
               ENDDO
            ENDDO
         ENDDO

         READ(10,FORM) P(NTYP)%QPAW(:,:,0)

         READ(10,'(1X,A1)') CSEL
         IF (CSEL=='t') THEN
            READ(10,*) P(NTYP)%QTOT
            READ(10,*)
         ELSE
            DO I=1,CHANNELS
               DO J=1,CHANNELS
                  P(NTYP)%QTOT(I,J)=0
               ENDDO
               P(NTYP)%QTOT(I,I)=1
            ENDDO
         ENDIF
         READ(10,FORM) P(NTYP)%QATO

         CALL READGRD(10,FORM,P(NTYP)%R%R,'g',IU0)
         CALL READGRD(10,FORM,P(NTYP)%POTAE,'a',IU0)
         P(NTYP)%POTAE_XCUPDATED=P(NTYP)%POTAE
         CALL READGRD(10,FORM,P(NTYP)%RHOAE,'c',IU0)
         READ(10,'(1X,A1)') CSEL
         IF (CSEL=='k') THEN
            CALL READGRD_NO_SELECTOR(10,FORM,P(NTYP)%TAUAE,'k',CSEL,IU0)
            READ(10,'(1X,A1)') CSEL
         ELSE
            P(NTYP)%TAUAE=0
         ENDIF
         IF (CSEL=='m') THEN
            CALL READGRD_NO_SELECTOR(10,FORM,P(NTYP)%TAUPS,'m',CSEL,IU0)
            READ(10,'(1X,A1)') CSEL
         ELSE
            P(NTYP)%TAUPS=0
         ENDIF
         LPOTPSC=.FALSE.
         IF (CSEL=='l') THEN
            CALL READGRD_NO_SELECTOR(10,FORM,P(NTYP)%POTPSC,'l',CSEL,IU0)
            LPOTPSC=.TRUE.
            READ(10,'(1X,A1)') CSEL
         ENDIF
         CALL READGRD_NO_SELECTOR(10,FORM,P(NTYP)%POTPS,'p',CSEL,IU0)
         CALL READGRD(10,FORM,P(NTYP)%RHOPS,'c',IU0)

         P(NTYP)%R%NMAX  =NMAX
         P(NTYP)%R%RSTART=P(NTYP)%R%R(1)
         P(NTYP)%R%REND  =P(NTYP)%R%R(NMAX)
         P(NTYP)%R%D     =(P(NTYP)%R%REND/P(NTYP)%R%RSTART)**(1._q/(NMAX-1))
         P(NTYP)%R%H     =LOG(P(NTYP)%R%D)
         P(NTYP)%R%RMAX=P(NTYP)%PSDMAX
 
         IF (.NOT. LPOTPSC) THEN
            CALL POTTORHO( P(NTYP)%ZVALF, NPSPTS, P(NTYP)%PSP(:,2), P(NTYP)%PSGMAX/NPSPTS, &
                 .TRUE. , NMAX, P(NTYP)%R%R ,  P(NTYP)%POTPSC )                        
         ENDIF

         DO I=1,CHANNELS
            CALL READGRD(10,FORM,P(NTYP)%WPS(:,I),'p',IU0)
            CALL READGRD(10,FORM,P(NTYP)%WAE(:,I),'a',IU0)
         ENDDO
         READ(10,*,ERR=600,END=600) CSEL
         IF (NWRITE>=0.AND.IU6>=0) &
              WRITE(IU6,*)'   PAW grid and wavefunctions read in'
      ENDIF aug

  600 CONTINUE

      IF (NWRITE>=0.AND.IU6>=0) THEN
      WRITE(IU6,*)
      WRITE(IU6,*)'  number of l-projection  operators is LMAX  =', &
     &            P(NTYP)%LMAX
      WRITE(IU6,*)'  number of lm-projection operators is LMMAX =', &
     &            P(NTYP)%LMMAX
      WRITE(IU6,*)
      ENDIF
      
      CALL RDPARS_ATOMIC_CONFIGURATION(ISDIM,STRING,IREAD,P(NTYP),IU6)

      CALL RDPARS_ENERGY(ISDIM,STRING,IREAD,P(NTYP),IU6)

      ENDDO forever

  100 CONTINUE


      NTYP=NTYP-1
      

      CLOSE(10)
      RETURN
      END SUBROUTINE

!******************* SUBROUTINE READGRD *******************************
!
! small helper routine to read a grid based entry from the POTCAR
! file
!
!**********************************************************************

      SUBROUTINE READGRD(IUNIT,FORM,A,STRING,IU0)
        IMPLICIT NONE
        INTEGER  :: IUNIT,IU0
        CHARACTER (LEN=*) :: FORM
        REAL(q)  :: A(:)
        CHARACTER (1) :: STRING,CSEL

        READ(IUNIT,'(1X,A1)') CSEL
        IF (CSEL /= STRING) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'READGRD: POTCAR file has an error: '
              WRITE(IU0,*)' expected: ',STRING,' found: ',CSEL
           ENDIF
           STOP
        ENDIF

        READ(IUNIT,FORM) A
      END SUBROUTINE READGRD

      ! selector already read, just compare CSEL and STRING

      SUBROUTINE READGRD_NO_SELECTOR(IUNIT,FORM,A,STRING,CSEL,IU0)
        IMPLICIT NONE
        INTEGER  :: IUNIT,IU0
        CHARACTER (LEN=*) :: FORM
        REAL(q)  :: A(:)
        CHARACTER (1) :: STRING,CSEL

        IF (CSEL /= STRING) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'READGRD: POTCAR file has an error: '
              WRITE(IU0,*)' expected: ',STRING,' found: ',CSEL
           ENDIF
           STOP
        ENDIF

        READ(IUNIT,FORM) A

      END SUBROUTINE READGRD_NO_SELECTOR

!******************* SUBROUTINE DEALLOC_PP ****************************
!
!  deallocate PP arrays
!  (I always try to do this in reverse order,
!   might help the memory subsystem to remerge free junks)
!**********************************************************************

      SUBROUTINE DEALLOC_PP(P,NTYPD)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (potcar) P(NTYPD)

      DO NTYP=NTYPD,1,-1
        CALL DEALLOC_PP1(P(NTYP))
      ENDDO
      RETURN

      END SUBROUTINE



      SUBROUTINE DEALLOC_PP1(P)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (potcar) P

      IF (ASSOCIATED(P%QDEP)) THEN
        DEALLOCATE(P%QDEP,P%NDEP)
      ENDIF
      DEALLOCATE(P%DION,P%QION, &
     &         P%LPS,P%NLPRO, &
     &         P%PSPNL,P%PSPRNL,P%E)

      ! Allocated in RDPARS_ATOMIC_CONFIGURATION
      IF(ASSOCIATED(P%ATOMIC_N)) THEN
         DEALLOCATE(P%ATOMIC_N,P%ATOMIC_L,P%ATOMIC_J,P%ATOMIC_OCC,P%ATOMIC_E)
      ENDIF

      DEALLOCATE(P%PSPRHO)

      IF (ASSOCIATED(P%PSPCOR)) THEN
        DEALLOCATE(P%PSPCOR)
      ENDIF
      DEALLOCATE(P%PSP)

      IF (ASSOCIATED(P%PSPNL_SPLINE)) THEN
         DEALLOCATE(P%PSPNL_SPLINE)
      ENDIF

      RETURN
      END SUBROUTINE

!******************* SUBROUTINE RDPARA ********************************
!
!  this subroutine  reads the implicit PSCTR description from the
!  POTCAR file (UNIT=IU)
!  information is stored in the string array
!   STRING(ISDIM)
!  the number of lines is returned in IREAD
!
!**********************************************************************

      SUBROUTINE RDPARA(IU,ISDIM,STRING,IREAD,IWRITE)
      USE mpimy
      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER (80) STRING(ISDIM)
      CHARACTER (6)  TAG

      I=0
  100 I=I+1
      IF (I>ISDIM) THEN
        WRITE(*,*)'ERROR: Description of pseudopotential is too long'
        STOP
      ENDIF
      READ(IU,'(A)') STRING(I)
      IF (STRING(I)(1:3)=='END') GOTO 120
      GOTO 100
  120 CONTINUE
      IREAD=I

      DO 200 I=1,IREAD
        TAG=STRING(I)(4:9)
        IF (TAG(1:5)=='Error') GOTO 210
  200 CONTINUE
  210 IWRITE=I-1

      RETURN
      END SUBROUTINE

!******************* SUBROUTINE RDPARS ********************************
!
!  this subroutine interprets the STRING array read in the
!  subroutine RDPARA
!  only a few number of items are retrieved from the STRING array
!   LREAL
!
!**********************************************************************


      SUBROUTINE RDPARS(ISDIM,STRING,IREAD,P,IU6)
      USE constant
      USE setexm
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (potcar) P
      CHARACTER (80) STRING(ISDIM)

      CHARACTER (6) TAG
      CHARACTER (80) VALUE
      LOGICAL LDUM
      EXTERNAL  LENGTH
!
!     initialize all values
!
      P%LREAL=.FALSE.
      P%LEXCH=0
      P%ZVALF=0
      P%POMASS=0
      P%RWIGS=0
      P%EATOM=0
      P%ENMAXA=0
      P%ENMINA=0
      P%QOPT1=0
      P%QOPT2=0
      P%EAUG=0
      P%EKECUT=0
      P%EKEERR=0
      P%EGGA=0
      P%DEXCCORE=0
      P%ELEMENT='  '
      P%RDEP=0
      P%VCA=1.0

      I=1
!  next line
  120 CONTINUE

      ITEM=1
!  next item in line
  110 CONTINUE

      IF (ITEM==1) THEN
        TAG=STRING(I)(4:9)
        VALUE=STRING(I)(13:80)
      ELSE
        TAG=STRING(I)(23:28)
        VALUE=STRING(I)(32:80)
      ENDIF

      EDUML=0
      EDUM=0
      IDUM=0
      LDUM=.FALSE.

      READ(VALUE,'(G10.4)',IOSTAT=IERR) EDUML
      READ(VALUE,'(G8.3)',IOSTAT=IERR) EDUM
      READ(VALUE,'(L8)',IOSTAT=IERR) LDUM
      READ(VALUE,'(I8)',IOSTAT=IERR) IDUM
!
!     set values
!
      L=LENGTH(TAG)
      IF (TAG(1:L)=='QCUT')   THEN
        P%LREAL=.TRUE.
        P%QOPT1=EDUM
        P%QOPT2=EDUM*2
      ENDIF
      IF (TAG(1:L)=='QGAM')   THEN
        P%LREAL=.TRUE.
        P%QOPT2=EDUM
      ENDIF
      IF (TAG(1:L)=='ZVAL')   THEN
         P%ZVALF  =EDUM
         P%ZVALF_ORIG=EDUM
      ENDIF
      IF (TAG(1:L)=='POMASS') P%POMASS=EDUM
      IF (TAG(1:L)=='RWIGS')  P%RWIGS =EDUM
      IF (TAG(1:L)=='RDEP')  THEN
         P%RDEP  =EDUM*AUTOA
      ENDIF
      IF (TAG(1:L)=='ENMAX')  THEN
        P%ENMAXA=EDUM
      ENDIF
      IF (TAG(1:L)=='ENMIN')  P%ENMINA=EDUM
      IF (TAG(1:L)=='EAUG')   P%EAUG  =EDUM
      IF (TAG(1:L)=='VCA')    P%VCA   =EDUM
      IF (TAG(1:L)=='EATOM')  P%EATOM =EDUML
      IF (TAG(1:L)=='DEXC')   P%DEXCCORE=EDUM
      IF (TAG(1:L)=='LEXCH')  THEN
        CALL EXTYP(VALUE(1:2),P%LEXCH)
      ENDIF
      IF (TAG(1:L)=='GGA')  THEN
        READ(STRING(I)(13:80),'(4G10.4)',IOSTAT=IERR) P%EGGA
      ENDIF
!  kinetic energy error lines
      IF (TAG(1:L)=='Error')  THEN
        NLINE=8
        I=I+1
        READ(STRING(I)(13:80),'(I8)',IOSTAT=IERR)   NDATA
        I=I+1
        READ(STRING(I)(13:80),'(2G8.3)',IOSTAT=IERR) START,STEP
        IF (NDATA/=NEKERR) GOTO 130
        DO I1=0,NDATA-1,NLINE
        I=I+1
        IM=MIN(NDATA-I1,NLINE)
        READ(STRING(I),'(8G10.4)',IOSTAT=IERR) (P%EKEERR(I1+I2),I2=1,IM)
        DO I2=1,IM
          P%EKECUT(I1+I2)=STEP**(I1+I2-1)*START
        ENDDO
        ENDDO
!  next item
      ENDIF
      IF (TAG(1:L)=='VRHFIN')  THEN
         P%ELEMENT(1:1)=STRING(I)(12:12)
         IF (STRING(I)(13:13).NE.':') P%ELEMENT(2:2)=STRING(I)(13:13)
      ENDIF
      ITEM=ITEM+1
      IF (ITEM<=2)  GOTO 110
      I=I+1
      IF (I<=IREAD) GOTO 120
      RETURN

!   error in kinetic energy lines
  130 CONTINUE
      IF (IU6>=0) &
      WRITE(*,*)'WARNING: can not read POTCAR file (kinetic energy)'

      END SUBROUTINE

!******************* SUBROUTINE RDPARS ********************************
!
!  this subroutine interprets the STRING array and
!  reads the eigenenergies at which the pseudisation was performed
!
!**********************************************************************


      SUBROUTINE RDPARS_ENERGY(ISDIM,STRING,IREAD,P,IU6)
      USE setexm
      IMPLICIT NONE

      INTEGER ISDIM
      INTEGER IREAD
      CHARACTER (80) STRING(ISDIM)
      TYPE (potcar) P
      INTEGER IU6
   ! local
      CHARACTER (6) TAG
      CHARACTER (80) VALUE
      INTEGER, PARAMETER :: LMAX=4, NMAX=4
      INTEGER I, J, ILOW, LN(0:LMAX), ITYP, LL
      REAL(q) E, ENERGY(NMAX,0:LMAX), RCUT

      DO I=1,IREAD
        TAG=STRING(I)(4:9)
        VALUE=STRING(I)(13:80)
        IF (TAG=='Descri')  EXIT
      ENDDO

      LN=0

      ILOW=I
      DO I=ILOW+2,IREAD
         TAG=STRING(I)(4:9)
         VALUE=STRING(I)
         IF (TAG=='Error') EXIT
         IF (STRING(I)(1:3)=='END') EXIT
         READ(VALUE,*) LL,E,ITYP,RCUT
         IF (LL>LMAX) THEN
            WRITE(*,*)'internal error in RDPARS_ENERGY: L is too large',LL
            STOP
         ENDIF
         LN(LL)=LN(LL)+1
         IF (LN(LL)>NMAX) THEN
            WRITE(*,*)'internal error in RDPARS_ENERGY:  too many entries for L=',LL
            STOP
         ENDIF
         ENERGY(LN(LL),LL)=E
      ENDDO

      DO I=1,P%LMAX
         LL=P%LPS(I)
         IF (LN(LL)==0) THEN
            WRITE(*,*)'internal error in RDPARS_ENERGY: not sufficient entries for L=',LL
            STOP
         ENDIF
         
         ! use bottom most entry from ENERGY
         P%E(I)=ENERGY(1,LL)
         ! remove one entry from ENERGY
         LN(LL)=LN(LL)-1
!        DO J=LN(LL),1,-1
         DO J=1,LN(LL)
            ENERGY(J,LL)=ENERGY(J+1,LL)
         ENDDO
      ENDDO

      END SUBROUTINE


!******************* SUBROUTINE RDPARS ********************************
!
!  this subroutine interprets the STRING array and
!  reads the occupation numbers and energies of the
!  atomic configuration
!
!**********************************************************************


      SUBROUTINE RDPARS_ATOMIC_CONFIGURATION(ISDIM,STRING,IREAD,P,IU6)
      USE setexm
      IMPLICIT NONE

      INTEGER ISDIM
      INTEGER IREAD
      CHARACTER (80) STRING(ISDIM)
      TYPE (potcar) P
      INTEGER IU6
   ! local
      CHARACTER (6) TAG
      CHARACTER (80) VALUE
      INTEGER, PARAMETER :: LMAX=4, NMAX=4
      INTEGER I, J, ILOW, LN(0:LMAX), ITYP, LL
      REAL(q) E, ENERGY(NMAX,0:LMAX), RCUT

       DO I=1,IREAD
         TAG=STRING(I)(4:9)
         VALUE=STRING(I)(13:80)
         IF (TAG=='Atomic')  EXIT
       ENDDO

       IF (I==IREAD+1) RETURN

       I=I+1
       VALUE=STRING(I)
       READ(VALUE,*) J

       ALLOCATE(P%ATOMIC_N(J),P%ATOMIC_L(J),P%ATOMIC_J(J),P%ATOMIC_OCC(J), &
               P%ATOMIC_E(J))

       ILOW=I
       J=1
       DO I=ILOW+2,IREAD

          TAG=STRING(I)(4:9)
          VALUE=STRING(I)
          IF (TAG=='Descri') EXIT
          READ(VALUE,*)P%ATOMIC_N(J),P%ATOMIC_L(J),P%ATOMIC_J(J),P%ATOMIC_E(J), &
               P%ATOMIC_OCC(J)
          J=J+1
       ENDDO

    END SUBROUTINE RDPARS_ATOMIC_CONFIGURATION


!******************* SUBROUTINE PSPOST ********************************
!
!  postprocessing of pseudopotential reader
!  checks
!
!**********************************************************************

      SUBROUTINE POST_PSEUDO(NTYPD,NTYP_PP,NTYP,NIONS,NITYP,VCA,P,INFO, &
     &        LREALD,ROPT,IDIOT,IU6,IU0,LMAX_CALC,L_NO_US,LSPIRAL)
      USE base
      USE ini
      USE constant
      USE setexm
      IMPLICIT NONE

      INTEGER NTYPD
      TYPE (potcar) P(NTYPD)
      TYPE (info_struct) INFO
      INTEGER NIONS         ! number of ions
      INTEGER NTYP,NTYP_PP  ! number of spec. on INCAR  / POTCAR
      INTEGER NITYP(NTYPD)  ! numer of ions of each species
      REAL(q) VCA(NTYPD)    ! weight for each atom species
      LOGICAL LREALD        ! default for LREAL from INCAR
      REAL(q)    ROPT(NTYPD)   ! cutoff for automatic real space opt
      LOGICAL L_NO_US       ! no US PP
      INTEGER LMAX_CALC     ! maximum L quantum number in PAW method
      LOGICAL LSPIRAL
! arryas for tutor call
      INTEGER IU6,IU0,IDIOT,ITUT(3),CDUM,LDUM
      REAL(q)    RTUT(3)
! temporary variables
      LOGICAL LREALT
      INTEGER NT,LEXCH_PP
      REAL(q) NVALEL
      REAL(q) QMAXL,QMAXNL,QMAXNL2,EERROR,ENAUG2,DUMMY
      REAL(q) SPLFIT(NEKERR,5)
      REAL(q) GMAX,G1,G2

      INTEGER, PARAMETER :: NFFT=32768

!-----------------------------------------------------------------------
! check number of species
!-----------------------------------------------------------------------
      IF (NTYP/=NTYP_PP) THEN
        IF (IU0>=0) &
        WRITE(IU0,*) 'ERROR: number of potentials on File POTCAR', &
     &              ' incompatible with number of species', &
     &             'INCAR :',NTYP,'POTCAR: ',NTYP_PP
        STOP
      ENDIF
!-----------------------------------------------------------------------
! all PP generated with the same XC-type
! also find out whether all PP are real space optimized
!-----------------------------------------------------------------------
      LEXCH_PP=P(1)%LEXCH
      DO NT=2,NTYP
         IF (LEXCH_PP /= P(NT)%LEXCH) THEN
         ITUT(1)=P(NT)%LEXCH
         ITUT(2)=LEXCH_PP
         ITUT(3)=NT
         CALL VTUTOR('E','DIFFERENT LDA-XC TYPES',RTUT,1, &
     &                ITUT,3,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('S','DIFFERENT LDA-XC TYPES',RTUT,1, &
     &               ITUT,3,CDUM,1,LDUM,1,IU0,IDIOT)
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!  set default value for LREAL if not given in INCAR
!-----------------------------------------------------------------------
      LREALT=P(1)%LREAL
      DO NT=2,NTYP
         LREALT=LREALT.AND.P(NT)%LREAL
      ENDDO
      IF (LREALD) INFO%LREAL=.FALSE.
!---- 'very large' cell and no real-space-projection scheme ???
      IF ((.NOT.INFO%LREAL).AND.(NIONS>16)) THEN
         IF (LREALT) THEN
            CALL VTUTOR('A','NO REAL-SPACE AND YOU COULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
            CALL VTUTOR('A','NO REAL-SPACE AND YOU COULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
         ELSE
            CALL VTUTOR('A','NO REAL-SPACE AND YOU SHOULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
            CALL VTUTOR('A','NO REAL-SPACE AND YOU SHOULD', &
     &                  RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
         ENDIF
      ENDIF
!---- 'very small' cell and still real-space projection scheme ???
      IF (INFO%LREAL.AND.(NIONS<=8)) THEN
         CALL VTUTOR('A','REAL-SPACE NOMORE RECOMMENDED', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('A','REAL-SPACE NOMORE RECOMMENDED', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF
!-----------------------------------------------------------------------
!  set energy cutoff and
!  check whether range of potentials is sufficient
!-----------------------------------------------------------------------
      IF (INFO%ENMAX==-1) THEN
      DO NT=1,NTYP
        IF (P(NT)%ENMAXA==0 &
     &   .OR.INFO%SZPREC(1:1)=='l'.AND.P(NT)%ENMINA==0 ) THEN
         IF (IU0>=0) &
         WRITE(IU0,*)'Fatal error! Could not find entry for ENMAX'// &
     &               ' on file INCAR. MUST be specified'
         STOP
        ENDIF
        IF  (INFO%SZPREC(1:1)=='l') THEN
          INFO%ENMAX=MAX(INFO%ENMAX,P(NT)%ENMINA)
        ELSE IF (INFO%SZPREC(1:1)=='h') THEN
          INFO%ENMAX=MAX(INFO%ENMAX,P(NT)%ENMAXA*1.25_q)
          IF (IU0>0) THEN
             WRITE(0,*) 'WARNING: for PREC=h ENMAX is automatically increase by 25 %'
             WRITE(0,*) '       this was not the case for versions prior to vasp.4.4'
          ENDIF
        ELSE
          INFO%ENMAX=MAX(INFO%ENMAX,P(NT)%ENMAXA)
        ENDIF
      ENDDO
      ENDIF
      IF (INFO%ENINI==-1) INFO%ENINI=INFO%ENMAX

      QMAXNL=   SQRT(INFO%ENMAX /RYTOEV)/AUTOA
      QMAXL = 2*SQRT(INFO%ENMAX /RYTOEV)/AUTOA

      DO NT=1,NTYP
      IF (ROPT(NT)/=0) THEN
        QMAXNL2    =QMAXNL*2.0
        IF  (INFO%SZPREC(1:1)=='h') THEN
           QMAXNL2    =QMAXNL*2.5
        ENDIF

        IF ( ABS(ROPT(NT)) > 0.1) THEN
          P(NT)%PSRMAX=12.5_q/ABS(QMAXNL)*ABS(ROPT(NT))**(1._q/3._q)
          IF  (INFO%SZPREC(1:1)=='h') &
          P(NT)%PSRMAX=10/ABS(QMAXNL)*ABS(ROPT(NT))**(1._q/3._q)
        ELSE
          P(NT)%PSRMAX=12.5_q/ABS(QMAXNL)*(0.5)**(1._q/3._q)
          IF  (INFO%SZPREC(1:1)=='h') &
          P(NT)%PSRMAX=10/ABS(QMAXNL)*(0.5)**(1._q/3._q)
        ENDIF

        IF (ROPT(NT)>0) THEN
           CALL OPTREAL(IU6,P(NT)%LMAX,P(NT)%LPS(1),NPSRNL, &
     &         P(NT)%PSMAXN/NPSRNL, &
     &         P(NT)%PSPRNL (1,1,1),P(NT)%PSPNL(0,1), &
     &         QMAXNL2,QMAXNL,P(NT)%PSRMAX,ROPT(NT)<0.1,ROPT(NT))
        ELSE
           CALL OPTREAL_NEW(IU6,P(NT)%LMAX,P(NT)%LPS(1),NPSRNL, &
     &         P(NT)%PSMAXN/NPSRNL, &
     &         P(NT)%PSPRNL (1,1,1),P(NT)%PSPNL(0,1), &
     &         QMAXNL2,QMAXNL,P(NT)%PSRMAX,ABS(ROPT(NT))<0.1,ABS(ROPT(NT)))
        ENDIF
        P(NT)%LREAL=.TRUE.
        P(NT)%QOPT1=QMAXNL *AUTOA
        P(NT)%QOPT2=QMAXNL2*AUTOA
      ENDIF

      IF (P(NT)%PSGMAX< QMAXL ) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WARNING: PSGMAX for local potential too small'
        IF (IU6>=0) THEN
        WRITE(IU6,*)'WARNING: PSGMAX for local potential too small'
        WRITE(IU6,*)'         PSGMAX should be >',QMAXL,' NTYP=',NT
        ENDIF
      ENDIF
      IF (P(NT)%PSMAXN<  QMAXNL ) THEN
        IF (IU0>=0) &
        WRITE(IU0,*)'WARNING: PSMAXN for non-local potential too small'
        IF (IU6>=0) THEN
        WRITE(IU6,*)'WARNING: PSMAXN for non-local potential too small'
        WRITE(IU6,*)'         PSMAXN should be >', QMAXNL,' NTYP=',NT
        ENDIF
      ENDIF
      IF (INFO%LREAL.AND. P(NT)%QOPT1/=0) THEN
        IF ((P(NT)%QOPT1**2*RYTOEV-INFO%ENMAX)/INFO%ENMAX &
     &      > 0.1_q) THEN
         RTUT(1)=P(NT)%QOPT1**2*RYTOEV
         RTUT(2)=INFO%ENMAX
         RTUT(3)=SQRT(INFO%ENMAX/RYTOEV)
         CALL VTUTOR('A','WRONG OPTIMZATION REAL-SPACE', &
     &               RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('A','WRONG OPTIMZATION REAL-SPACE', &
     &               RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
        ENDIF
      ENDIF
      ENDDO

      LREALT=P(1)%LREAL
      DO NT=2,NTYP
         LREALT=LREALT.AND.P(NT)%LREAL
      ENDDO
!---- complain about using non-optimized projectors if LREAL=.T.
      IF (INFO%LREAL.AND.(.NOT.LREALT)) THEN
         CALL VTUTOR('W','REAL-SPACE WITHOUT OPTIMIZATION', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('W','REAL-SPACE WITHOUT OPTIMIZATION', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF
!-----------------------------------------------------------------------
! switch to specific LDA/GGA if tag was found on INCAR file
!-----------------------------------------------------------------------
      IF (LEXCH==-1) LEXCH=LEXCH_PP
      IF (LEXCH/=LEXCH_PP) THEN
         IF (ISGGA()) THEN
!  add GGA corrections to EATOM
!  for some of the LDA potential this is possible
            DO NT=1,NTYP
               P(NT)%EATOM=P(NT)%EATOM-P(NT)%EGGA(MOD(LEXCH-4,4)+1)
            ENDDO
         ENDIF
         CALL VTUTOR('A','ENFORCED LDA', &
     &        RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('A','ENFORCED LDA', &
     &        RTUT,3,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF
!-----------------------------------------------------------------------
!  set number of electrons and total atomic energy
!-----------------------------------------------------------------------
      INFO%EALLAT=0
      NVALEL=0
      ENAUG2=0
      DO NT=1,NTYP
       ENAUG2=MAX(P(NT)%EAUG,ENAUG2)
       EERROR=0
       IF (IU6>=0) &
       WRITE(IU6,449) P(NT)%SZNAMP,NT,-P(NT)%EATOM

  449  FORMAT( A40,':'/ &
     & ' energy of atom ',I2,'       EATOM=',F10.4)
       IF ((P(NT)%EKEERR(1)/=0 .OR. P(NT)%EKECUT(1)/=0)) &
     & THEN
       CALL SPLCPY(P(NT)%EKECUT,P(NT)%EKEERR,SPLFIT,NEKERR,NEKERR,10E30_q)
       IF (LSPIRAL) THEN
          CALL SPLVAL(INFO%ENINI,EERROR,DUMMY,SPLFIT,NEKERR,NEKERR)
       ELSE
          CALL SPLVAL(INFO%ENMAX,EERROR,DUMMY,SPLFIT,NEKERR,NEKERR)
       ENDIF
       IF (IU6>=0) &
       WRITE(IU6,448) EERROR

  448  FORMAT( &
     & ' kinetic energy error for atom=',F10.4, &
     & ' (will be added to EATOM!!)')
       ENDIF
       ! energy including convergence corrections
       P(NT)%EATOM_CORRECTED=P(NT)%EATOM-EERROR
       INFO%EALLAT=INFO%EALLAT+NITYP(NT)*(P(NT)%EATOM_CORRECTED)*VCA(NT)
       NVALEL=NVALEL+NITYP(NT)*P(NT)%ZVALF*VCA(NT)
      ENDDO
      IF (INFO%ENAUG==-1) INFO%ENAUG=ENAUG2

      IF (INFO%NELECT==0) INFO%NELECT= NVALEL
      IF (IU6>=0) WRITE(IU6,*)
!-----------------------------------------------------------------------
!  set truncation for LMAX for PAW
!-----------------------------------------------------------------------
      IF (LMAX_CALC >= 0) THEN
         P%LMAX_CALC=MIN( LMAX_CALC, P%LMAX_CALC)
      ENDIF

      L_NO_US=.TRUE.
      DO NT=1,NTYP
         IF (ASSOCIATED( P(NT)%QDEP )) THEN
            L_NO_US=.FALSE.
         ENDIF
      ENDDO

      IF(INFO%TURBO>=1)THEN
         DO NT=1,NTYP
            GMAX=SQRT(INFO%ENMAX/HSQDTM)
            IF (INFO%SZPREC(1:1)=='s') THEN
               G1=GMAX
               G2=MIN(3.0_q*GMAX,P(NT)%PSGMAX)/2
            ELSE
               G1=2.0_q*GMAX
               G2=MIN(3.0_q*GMAX,P(NT)%PSGMAX)
            ENDIF
            ALLOCATE(P(NT)%RHOSPL(NFFT/2,5))
            CALL POTTORHO_FILTER(G1, G2, P(NT)%ZVALF, NPSPTS, P(NT)%PSP(:,2), &
                 P(NT)%PSGMAX/NPSPTS, NFFT, P(NT)%RHOSPL, P(NT)%RCUTRHO, P(NT)%ESELF, IU6)
         ENDDO
      ENDIF

      IF (INFO%NLSPLINE) THEN
         CALL VTUTOR('W','SPLINE INTERPOLATE PROJECTORS', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU6,IDIOT)
         CALL VTUTOR('W','SPLINE INTERPOLATE PROJECTORS', &
     &               RTUT,1,ITUT,1,CDUM,1,LDUM,1,IU0,IDIOT)
      ENDIF

      END SUBROUTINE

!**********************************************************************
!
! determine the maximum L quantum number (LDIMP)
!
!**********************************************************************

      SUBROUTINE LDIM_PSEUDO(LORBIT, NTYPD, P, LDIMP, LMDIMP)
      USE base
      USE  ini
      USE constant

      IMPLICIT NONE

      INTEGER LORBIT       ! how do calculate partial dos
      INTEGER NTYPD        ! how many types are there
      INTEGER LDIMP,LMDIMP ! required dimension of arrays
      TYPE (potcar) P(NTYPD)
! local
      INTEGER MAXIMUM_L,NTYP,I

      MAXIMUM_L=0

      DO NTYP=1,NTYPD
         DO I=1,P(NTYP)%LMAX
            MAXIMUM_L = MAX(MAXIMUM_L,P(NTYP)%LPS(I))
         ENDDO
      ENDDO
      
      ! if the maximum L is larger than 2 we have to set LDIMP

      IF ( MAXIMUM_L > 2 ) THEN
         LDIMP= MAXIMUM_L+1
      ELSE
         LDIMP=3
      ENDIF
      LMDIMP=LDIMP*LDIMP

      END SUBROUTINE LDIM_PSEUDO

      END MODULE
