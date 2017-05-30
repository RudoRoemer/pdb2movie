      program diagstd
      implicit none
      integer natmax, ndim
c======================================================================

c     Diagonalization of a matrix, real, symmetrical

c     Diagonalization routine: TQLI (EISPACK).
c    -simple, public domain, but slow (all eigenvectors are computed).
 
c======================================================================
 
c     INPUT matrix filename (expected from PDBMAT):
c     *********************************************
c     CERFACS -Formatted : pdbmat.sdijf
c     CERFBIN -Binary    : pdbmat.sdijb

c     Input matrix format: i, j, element-ij-non-nul.

c     OUTPUT:
c     *******
c     Eigenvector filename: pdbmat.eigenfacs

c.......................................................................

c     MEMORY LIMITS:
c     **************
c     NATMAX: Maximum number of atoms allowed.

      parameter(natmax=4000,ndim=3*natmax)
 
c.......................................................................
c     YHS-Sep-2002: Premiere version, from Diagijr v1.13 (Bordeaux).
c.......................................................................
      logical qcrois, qexist, qinterr
      integer evord(ndim), i, ii, j, jj, k, natom, nbig, nord, 
     .        nredond, ntrace, nvec, nvecout, nunit, rdunit, 
     .        unmess, unmodes
      double precision amat(ndim,ndim), ev(ndim), evsort(ndim),
     .       matrd, trace, work(ndim)
      character cformat*20, cstatus*20, eige*4, matrice*20, nomfich*20, 
     .       program*9, progrer*12, progrwn*12, version*32
c.......................................................................
      version=' Version 1.05, August 2004.'
 
      program=' Diagstd>'
      progrer='%Diagstd-Er:'
      progrwn='%Diagstd-Wn:'
 
c     Sortie standard: 
      unmess=6
      write(unmess,'(2A)') program,version

c     Eigenvector are given by increasing eigenvalues (LOWE).
c     Or by decreasing values (HIGH).
      eige='LOWE'

      nvec=3*natmax
 
c     Detection de la matrice d'entree:
c     --------------------------------
      nunit=10
      rdunit=nunit
      nunit=nunit+1
      cformat='UNFORMATTED'
      cstatus='OLD'
      nomfich='pdbmat.sdijb'
      call openam(nomfich,cformat,cstatus,
     .     rdunit,.false.,
     .     qinterr,qexist)
      if (qexist) then
          matrice='CERFBIN'
      else
      nomfich='pdbmat.sdijf'
      cformat='FORMATTED'
      call openam(nomfich,cformat,cstatus,
     .     rdunit,.false.,
     .     qinterr,qexist)
      if (qexist) then
          matrice='CERFACS'
      else
          matrice='NONE'
          write(unmess,'(2A)') progrer,
     .' Matrice des derivees secondes non trouvee.',
     .' Noms possibles: ',
     .' pdbmat.sdijf  (CERFACS) -formattee, creuse, format libre.',
     .' pdbmat.sdijb  (CERFBIN) -non formattee, creuse, format libre.'
          stop
      endif
      endif
 
      write(unmess,'(3A)') program,
     .    ' Matrix to be read from file: ',nomfich
 
c     Lecture matrice d'entree (CERFACS, CERFBIN):
c     --------------------------------------------
 
      write(unmess,'(/2A)') program,
     .     ' Lecture de la matrice -Format CERFACS- '
 
c     1) Ordre de la matrice, nombre de lignes:
c     -----------------------------------------
      k=0
      nord=0
  90  continue
      if (matrice.eq.'CERFACS') then
      read(rdunit,*,end=100) i,j
      else
      read(rdunit,end=100) i,j
      endif
      k=k+1
      if (i.le.0.or.j.le.0) then
          write(unmess,'(/2A,I9,2(A,I6))')
     .    progrer,' in ligne: ',k,' I= ',i,' J= ',j
          stop
      endif
      if (i.gt.nord) nord=i
      if (j.gt.nord) nord=j
      goto 90
 100  continue
 
      write(unmess,'(2A,I9)')
     .     program,' Matrix dimension  (Nord)  =',nord
      write(unmess,'(2A,I9)')
     .     program,' Number of non-zero elements',k
 
      natom=nord/3
      if (natom.gt.natmax.or.nord.gt.ndim) then
          write(unmess,'(2A)')
     .    progrer,' Matrix can not be read.'
          if (natom.gt.natmax) write(unmess,'(2(A,I9))')
     .   ' Natom= ',natom,' > natmax= ',natmax
          if (nord.gt.ndim) write(unmess,'(2(A,I9))')
     .   ' Nord=  ',nord,' > Ndim=  ',ndim
          stop
      endif
 
c     2) Lecture de la matrice:
c     -------------------------
      rewind(rdunit)
 
      nredond=0
      ntrace=0
      trace=0.d0
      nbig=0
      do i=1,nord
        do j=1,nord
         amat(i,j)=0.d0
        enddo
      enddo

      do jj=1,k
         if (matrice.eq.'CERFACS') then
         read(rdunit,*,err=95) i,j,matrd
         else
         read(rdunit,err=95) i,j,matrd
         endif
 
         if (dabs(matrd).gt.0.d0) then
             amat(i,j)=matrd
             amat(j,i)=matrd
             if (i.eq.j) then 
                trace=trace+matrd
                ntrace=ntrace+1
             endif
             if (matrd.gt.1E+10) then
                 nbig=nbig+1
                 if (nbig.lt.10) then
                     write(unmess,'(2A,2I12,A,G12.3)') 
     .               progrwn,' Element: ',i,j,' = ',matrd
                 else 
                     if (nbig.eq.10) write(unmess,*) '...'
                 endif
             endif
         else
             nredond=nredond+1
         endif
      enddo
      goto 105
  95  continue
      write(unmess,'(2A,I6)')
     .     progrer,' while reading ligne ',k
      write(unmess,'(2I6,F16.8)') ' i, j, matrd= ',i,j,matrd
      stop
 105  continue
 
      write(unmess,'(2A,I9)') program,
     .    ' Nb of elements found twice:',nredond
      if (nredond.gt.0) 
     .write(unmess,'(2A/)') progrwn,' Ok ?'
      write(unmess,'(2A,I9)') program,
     .    ' Nb of elements    > 1E+10 :',nbig
      if (nbig.gt.0) 
     .write(unmess,'(2A/)') progrwn,' Ok ?'
      write(unmess,'(2A,F31.7)') program,
     .    ' Matrix trace:',trace
      if (ntrace.gt.0)
     .write(unmess,'(2A,I11)') progrwn,
     .    ' Nb on non-zero elements there:',ntrace
 
c     Diagonalisation:
c     ----------------
      nomfich='pdbmat.eigenfacs'
      cformat='FORMATTED'
      cstatus='ove'
      unmodes=nunit
      nunit=nunit+1
      call openam(nomfich,cformat,cstatus,
     .     unmodes,.true.,
     .     qinterr,qexist)
 
      write(unmess,'(/2A)') program,' Diagonalization.'
 
      if (nvec.gt.nord) nvec=nord
      write(unmess,'(A,I6,A)') program,
     .      nvec,' eigenvectors are about to be computed. '
      nvecout=nvec

c     Initialisations:
      do i=1,ndim
         ev(i)=0.d0
      enddo
 
c     Eigenvalues/Matrix Diagonalization

c     The following routines (based on the original EISPACK library) 
c     perform a diagonalization of a real symmetric matrix based 
c     on the QL algorithm. 

      CALL TRED2(amat,nord,ndim,ev,work)
      CALL TQLI(ev,work,nord,ndim,amat)
 
      trace=0.d0
      do i=1,nvec
         trace=trace+ev(i)
      enddo
      write(unmess,'(/2A,F24.7)') program,
     .     ' Sum of eigenvalues =',trace
 
c     Trier par ordre croissant ou decroissant:
      qcrois=.true.
      if (eige.eq.'HIGH') qcrois=.false.

      call trier(ev,nvec,ndim,evsort,evord,qcrois)

      write(unmess,'(/2A/(5F15.7))') program,
     .    ' Eigenvalues: ',(ev(evord(i)),i=1,nvecout)
      WRITE(unmess,'(/2A/(5F15.7))') program,
     .    ' Frequencies (cm-1, '//
     .     'if the matrix is a hessien in CHARMM units):',
     .    (sqrt(dabs(ev(evord(i))))*108.591365,i=1,nvecout)
 
c     Ecriture des modes normaux au format 'CERFACS':
c     -----------------------------------------------
      do j=1,nvecout
         i=evord(j)
         write(unmodes,'(A,I5,7X,A,1PG12.4)') 
     .       ' VECTOR',j,'VALUE',ev(i)
         write(unmodes,'(1X,35(1H-))') 
         write(unmodes,'(3(1PG12.4))') 
     .        (amat(k,i),k=1,nord)
      enddo
 
      write(unmess,'(/2A)')
     .      program,' Normal end.'
 
      stop 
      end
c----------------------------------------------------------------
      SUBROUTINE openam(namfil,cformat,cstatus,unit,qverbos,
     .                  qinterr,qexist)
c
c     Ouverture d'un fichier de nom NAMFIL, sur l'unite UNIT,
c     a priori suite a une interrogation...
c
c     input:
c        namfil: nom du fichier a ouvrir. 
c        "stop", "end", "fin", "quit" : arretent le programme.
c        cstatus: mots-cles fortran... ou "OVE" pour overwrite.
c     output: 
c        qexist: flag / existence du fichier 
c        qinterr: Pas de nom pour le fichier cherche.
c
c     YHS-oct-93
c     YHS-jan-95
c I/O:
      logical qinterr, qverbos, qexist
      integer unit
      character*(*) namfil, cformat, cstatus
c Local
      character*132 ordrunix
c begin:
      if (cstatus.eq.'old') cstatus='OLD'
      if (cstatus.eq.'new') cstatus='NEW'
      if (cstatus.eq.'ove') cstatus='OVE'
      if (cstatus.eq.'unknown') cstatus='UNKNOWN'
c
      qinterr=.false.
      qexist=.false.
c
      if (namfil.eq.' ') then 
          qinterr=.true.
          write(6,'(A)') '%Openam-Err> No filename.'
          return
      endif
c
      if (namfil.eq.'stop'.or.namfil.eq.'end'                         
     &    .or.namfil.eq.'fin'.or.namfil.eq.'quit') then 
         write(*,*) 'Openam> Program is stopping on user request.'
         stop                                                                   
      endif 
c
c     Checks if filename is consistent with the opening:
c
      inquire(file=namfil,exist=qexist)
      if (.not.qexist.and.cstatus.eq.'OLD') then
          qinterr=.true.
          if (qverbos) write(6,'(A)') '%Openam-Err> File not found.'
          return
      endif
c
      if (qexist.and.cstatus.eq.'NEW') then
         write(*,'(/A)') 
     .      '%Openam-Err> This file exists:',namfil
         stop
      else if (qexist.and.cstatus.eq.'OVE') then
         ordrunix='rm '//namfil
         call system(ordrunix)
      endif
      if (cstatus.eq.'OVE') cstatus='NEW'
c                                                                   
      if (qverbos) then
         write(*,'(/A,I6,A)')
     .           ' Openam> file on opening on unit ',unit,':'
         write(*,*) namfil
      endif
      open(file=namfil,form=cformat,
     .     status=cstatus,unit=unit)                
c        
      return  
      end

      SUBROUTINE TRED2(A,N,NP,D,E)

c     Reduce the matrix to tridiagonal form.

      integer i, j, k, l, n, np
      double precision A(NP,NP), D(NP), E(NP), f, g, h, hh,
     .  scale

      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
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
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
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
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END

      SUBROUTINE TQLI(D,E,N,NP,Z)

c     Finds the eigenvalues and eigenvectors of a tridiagonal matrix:

      integer i, iter, j, k, l, m, n, np
      double precision b, c, D(NP), dd, E(NP), f, g, p, r, s,
     .       Z(NP,NP)
      
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
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
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END

      subroutine trier(y,npoint,nmax,ysort,iord,qcrois)
c
c     Tri par ordre croissant (qcrois=T) ou non.
c     YHS-Jun-2002: Premiere version (Bordeaux).
c     YHS-Sep-2002: Derniere version (Bordeaux).

      implicit none
      logical qcrois
      integer i, icur, iord(*), j, nmax, npoint
      double precision y(*), ycur, ysort(*)
      character progrer*10

      progrer='%Trier-Er>'

      if (npoint.gt.nmax) then
          write(6,'(A,I9,A,I9,A)') progrer,npoint,
     .  ' points to be sorted, i.e., more than ',nmax,' Sorry.'
          stop
      endif

      do i=1,npoint
         ysort(i)=y(i)
         iord(i)=i
      enddo

      do i=1,npoint
        do j=1,npoint
          if (qcrois) then
            if (ysort(i).lt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          else
            if (ysort(i).gt.ysort(j)) then
                ycur=ysort(i)
                icur=iord(i)
                ysort(i)=ysort(j)
                ysort(j)=ycur
                iord(i)=iord(j)
                iord(j)=icur
            endif
          endif
        enddo
      enddo

      return
      end
