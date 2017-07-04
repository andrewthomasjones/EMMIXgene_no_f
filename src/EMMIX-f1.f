      PROGRAM EMMIXlite
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     INCLUDE 'EMMIXFAC.max'
      INCLUDE 'EMMIX-f1.max'

      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG,        ! Number of components or groups
     &        NCOV       ! Covariance structure
      REAL X(MNIND,MNATT)  ! Data set/sample
      INTEGER SUB(2+MNATT*10)

       DOUBLE PRECISION TOLS(2)
       INTEGER MAXITS(2),N(MAXNG),
     &          FACT(3)        ! Control mixtures of Factor analyzers
       CHARACTER*255 INFYLE
       CHARACTER*255 OFYLE
       CHARACTER*255 MODFYLE
       CHARACTER*255 MFYLE
       DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &           XMU(MAXNG,MNATT),
     &           T(MAXNG)
       DIMENSION XM(MNATT),XVR(MNATT)

        TOLS(1)=TAUTO
	MAXITS(1)=MITAUT
	TOLS(2)=TFINAL
	MAXITS(2)=MITFIN


C    Check to see if the random number generator works (try two formats)
      CALL DETERRANDOM(IER)
	              
      CALL SETUP(NIND,NATT,NCOV,NG0,NG1,NKMEANS,NRANDS,RDSB,INFYLE,
     & OFYLE,MFYLE,MODFYLE,OPTION,START,TOLS,SUB,SCAL,FACT,NATTQ,
     & REV,IER)

      OPEN (UNIT=21,FILE=INFYLE,STATUS = 'OLD',ERR=505)
      OPEN (UNIT=22,FILE=OFYLE,STATUS = 'UNKNOWN')
      OPEN (UNIT=23,FILE=MFYLE,STATUS = 'UNKNOWN')
      OPEN (UNIT=24,FILE=MODFYLE,STATUS = 'UNKNOWN')
       
      WRITE(22,*) '--------------------------------------'
      WRITE(22,*) 'Output from analysis of ',INFYLE
      WRITE(22,*) '--------------------------------------'
      WRITE(22,*) 'Number of points=',NIND
      WRITE(22,*) 'Number of Variables=',NATT
      WRITE(22,*) 'Number of Groups=',NG0
      WRITE(22,*) '--------------------------------------'
      WRITE(22,*) 

      CALL READX(NIND,NATT,X,SUB,REV)
      IF (SCAL.EQ.1) THEN
	CALL KSTAND(NIND,NATT,X,XM,XVR)
        DO 220 KK=1,NIND
          DO 220 J=1,NATT
           X(KK,J)=(DBLE(X(KK,J))-XM(J))/SQRT(XVR(J))
220     CONTINUE
      ENDIF

      IF (START.EQ.2) THEN
        CALL READPARA(NATT,NG0,XMU,XVAR,T,SUB)
      ELSEIF (START.EQ.1) THEN
        NATT=SUB(2)
        CALL READALLOC(NIND,NATT,NG0,XMU,XVAR,T,X)
      ENDIF
       
      NATT=SUB(2)
      CALL MAIN(NIND,NATT,NCOV,X,NG0,NG1,XMU,XVAR,T,MAXITS,
     & TOLS,AIT,OPTION,START,NRANDS,NKMEANS,FACT,RDSB,NATTQ)

      CLOSE(21)
      CLOSE(22)
505   CONTINUE
      END


      SUBROUTINE READX(NIND,NATT,X,SUB,REV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
        REAL X(MNIND,MNATT)
        REAL TEMP(MNATT*10)
        INTEGER SUB(2+MNATT*10)
        IF ((SUB(1).EQ.100).AND.(SUB(2).EQ.NATT)) THEN
	 IF (REV.EQ.0)  THEN
          DO 70 I=1,NIND
 	    READ(21,*) (X(I,J),J=1,NATT)
70        CONTINUE
	 ELSEIF (REV.EQ.1) THEN
          DO 170 I=1,NIND
 	    READ(21,*) (X(J,I),J=1,NATT)
170        CONTINUE
	   SUB(2)=NIND
	   NIND=NATT
	   NATT=SUB(2)
	   IF (NIND.GT.MNIND) THEN
	    WRITE (*,*)'Adjust EMMIX-f1.max'
	    stop
           ENDIF	     
	   IF (NATT.GT.MNATT) THEN
	    WRITE (*,*)'Adjust EMMIX-f1.max'
	    stop
           ENDIF	     
	 WRITE (*,*) 'Transposed new dim= ',NIND,NATT
	 WRITE (22,*) 'Transposed new dim= ',NIND,NATT
	 ENDIF
        ELSEIF (SUB(1).EQ.100) THEN 
          DO 71 I=1,NIND
            READ(21,*) (TEMP(J),J=1,NATT)
             DO 72 J=1,SUB(2)
              X(I,J)=TEMP(SUB(2+J))
72        CONTINUE
71        CONTINUE
        ELSE
         COUNT=0
         DO 100 I=1,NIND
          READ(21,*) (TEMP(J),J=1,NATT)
          R=RANDNUM()
          IF (R*100.0.LT.FLOAT(SUB(1))) THEN  
            COUNT=COUNT+1
            IF (SUB(2).NE.NATT) THEN 
             DO 80 J=1,SUB(2)
              X(COUNT,J)=TEMP(SUB(2+J))
80           CONTINUE
            ELSE
             DO 90 J=1,NATT
              X(COUNT,J)=TEMP(J)
90           CONTINUE
            ENDIF
          ENDIF
100	 CONTINUE
        ENDIF

        IF (SUB(1).NE.100) THEN
         NIND=COUNT
        ENDIF
         
        WRITE (*,*)  ' - Read in of data sucessful'
        WRITE (22,*) ' - Read in of data sucessful'
	RETURN
      END

      SUBROUTINE READPARA(NATT,NG,XMU,XVAR,T,SUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      INTEGER NATT,      ! Number of attributes/variables/dimensions
     &        NG         ! Number of components or groups
      INTEGER SUB(2+MNATT*10)
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          XMU(MAXNG,MNATT),
     &          T(MAXNG)
      DO 400 K=1,NG
        READ (21,*) (XMU(K,J),J=1,NATT)
        write (*,*) '  - Read in mean for group ',K
        write (22,*) '  - Read in mean for group ',K
        DO 404 I=1,NATT
         READ (21,*) (XVAR(K,IC(I,J)),J=1,I)
404     CONTINUE
        write (*,*) '  - Read in covariance matrix for group ',K
        write (22,*) '  - Read in covariance matrix for group ',K
400   CONTINUE
      READ (21,*) (T(K),K=1,NG)
        write (*,*) '  - Read in mixing proportions'
        write (22,*) '  - Read in mixing proportions'
        WRITE (*,*)  ' - Read in of parameters sucessful'
        WRITE (22,*) ' - Read in of parameters sucessful'
	RETURN
	END

      SUBROUTINE READALLOC(NIND,NATT,NG0,XMU,XVAR,T,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      REAL X(MNIND,MNATT)
      INTEGER NIND,
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG,        ! Number of components or groups
     &        N(MAXNG)
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          XMU(MAXNG,MNATT),
     &          T(MAXNG)
	     DO 224 K=1,NG0
	      N(K)=0
	      DO 224 I=1,NATT
	       XMU(K,I)=0.0
	       DO 224 L=1,I
               XVAR(K,IC(I,L))=0.0
224        CONTINUE

 
	     DO 227 I=1,NIND
	      READ (21,*) IDT
	      N(IDT)=N(IDT)+1
	      DO 226 J=1,NATT 
	       XMU(IDT,J)=XMU(IDT,J)+DBLE(X(I,J))
	       DO 226 L=1,J
	         XVAR(IDT,IC(J,L))=XVAR(IDT,IC(J,L))+DBLE(X(I,J)*X(I,L))
226          CONTINUE
227        CONTINUE
         DO 213 KK=1,NG0
           T(KK)=FLOAT(N(KK))/FLOAT(NIND)
             DO 213 JJ=1,NATT
              XMU(KK,JJ)=XMU(KK,JJ)/(FLOAT(N(KK)))
              DO 213 II=1,JJ
                XVAR(KK,IC(JJ,II))=
     &          XVAR(KK,IC(JJ,II))-(XMU(KK,II)*XMU(KK,JJ)*FLOAT(N(KK)))
                XVAR(KK,IC(JJ,II))=XVAR(KK,IC(JJ,II))/FLOAT(N(KK)-1)
212       CONTINUE	
213       CONTINUE	
        WRITE (*,*)  ' - Read in of allocation sucessful'
        WRITE (22,*) ' - Read in of allocation sucessful'
	RETURN
	END

      SUBROUTINE MAIN(NIND,NATT,NCOV,X,NG0,NG1,XMU,XVAR,T,MAXITS,
     & TOLS,AIT,OPTION,START,NRANDS,NKMEANS,FACT,RDSB,NATTQ)

C   PURPOSE
C     This is the main subroutine which controls the program
C     and branches according to the options specified by FLAGS.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'

C    Global Parameters
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO

C   INPUT PARAMETERS
      INTEGER   FLAGS(40),     ! Flags summarising options chosen (see below)
     &          FYLENO,        ! File id for output file
     &          NIND,          ! Number of points (samples, or observations)
     &          NATT,          ! Number of dimensions (variates or attributes)
     &          NG0,           ! Minimum number of groups to be tested
     &          NG1,           ! Maximum number of groups to be tested
     &          NCOV,          ! Structure of covariance matrices
     &          IX,IY,IZ,      ! Random seeds
     &          FACT(3)        ! Control mixtures of Factor analyzers

      REAL X(MNIND,MNATT) ! Data or sample.
      DIMENSION XMU(MAXNG,MNATT),
     &          XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          T(MAXNG),
     &          TXML(MITER)
      DIMENSION D(MAXNG,MNATT),
     &          B(MAXNG,MNATT,MNATTQ)


      DOUBLE PRECISION TOLS(2)
      INTEGER MAXITS(2)
       IF (START.NE.3) THEN
         MODE=2
	 SETUPON=1
	 FORMOT=1
         CALL EM(NIND,NATT,NG1,NCOV,X,XMU,XVAR,T,TXML,TOLS(2),
     &          MAXITS(2),MODE,FACT,B,D,NATTQ,FORMOT,SETUPON,IER)
         DO 7010 K=1,NG1
            write (24,*) (XMU(K,I),I=1,NATT)
            DO 7011 I=1,NATT
              write (24,*) (B(K,I,J),J=1,NATTQ)
7011        CONTINUE
            write (24,*) (D(K,J),J=1,NATT)
7010     CONTINUE
         write (24,*) (T(K),K=1,NG1)
         RETURN
       ENDIF        

       IF (OPTION.EQ.2) THEN
         IF (START.EQ.3) THEN
          CALL AUTO(NIND,NATT,NCOV,X,NG1,XMU,XVAR,T,MAXITS(1),
     &              NRANDS,NKMEANS,TOLS(1),AIT,RDSB,FACT,B,D,NATTQ)
         ENDIF   
         MODE=2
	 SETUPON=0
	 FORMOT=1
         CALL EM(NIND,NATT,NG1,NCOV,X,XMU,XVAR,T,TXML,TOLS(2),
     &          MAXITS(2),MODE,FACT,B,D,NATTQ,FORMOT,SETUPON,IER)
         CALL ESTEP(NIND,NATT,NCOV,NG1,X,XMU,V,DV,T,W,
     &            XLOGL,IER)
      
         DO 8010 K=1,NG1
            write (24,*) (XMU(K,I),I=1,NATT)
            DO 8011 I=1,NATT
              write (24,*) (B(K,I,J),J=1,NATTQ)
8011        CONTINUE
            write (24,*) (D(K,J),J=1,NATT)
8010     CONTINUE
         write (24,*) (T(K),K=1,NG1)

	     
       ELSEIF (OPTION.EQ.3) THEN
        DO 200 K=NG0,NG1
           CALL AUTO(NIND,NATT,NCOV,X,K,XMU,XVAR,T,MAXITS(1),
     &      NRANDS,NKMEANS,TOLS(1),AIT,RDSB,FACT,B,D,NATTQ)
           MODE=2
	 SETUPON=0
	  FORMOT=1
           CALL EM(NIND,NATT,K,NCOV,X,XMU,XVAR,T,TXML,TOLS(2),
     &      MAXITS(2),MODE,FACT,B,D,NATTQ,FORMOT,SETUPON,IER)
200     CONTINUE
       ENDIF
 
         DO 9010 K=1,NG1
            write (24,*) (XMU(K,I),I=1,NATT)
            DO 9011 I=1,NATT
              write (24,*) (B(K,I,J),J=1,NATTQ)
9011        CONTINUE
            write (24,*) (D(K,J),J=1,NATT)
9010     CONTINUE
         write (24,*) (T(K),K=1,NG1)

       RETURN
       END

      SUBROUTINE AUTO(NIND,NATT,NCOV,X,NG,XMU,XVAR,T,MAXIT,
     &         NRANDS,NKMEANS,TOL,AIT,RDSB,FACT,B,D,NATTQ)
C    Constants that define array sizes at compilation time.
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INCLUDE 'EMMIX-f1.max'
	COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	INTEGER FACT(3)
	REAL X(MNIND,MNATT)
	DIMENSION XMU(MAXNG,MNATT),
     &            XVAR(MAXNG,MNATT*(MNATT+1)/2)
	DIMENSION XVAR2(MAXNG,MNATT*(MNATT+1)/2),
     &            XMU2(MAXNG,MNATT),T2(MAXNG),SLIKS(MITER)
	DIMENSION D(MAXNG,MNATT),
     &            B(MAXNG,MNATT,MNATTQ)


	  write (22,*) '  --------------------------' 
	TMP=1.0 
	NS=0      
	DO I=1,NRANDS      
          write (*,*) '  Calling Random Start ',I
	  write (22,*) '  Random Start ',I
          write (22,*) '     Seeds ', IX,IY,IZ
	  CALL RANDCON(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2
     &             ,RDSB)
	  CALL FIT(NIND,NATT,NCOV,X,NG,XMU,XVAR,T,XMU2,XVAR2,T2,
     &   	   TOL,MAXIT,AIT,XLIKE,TLIKE,TMP,FACT,B,D,NATTQ,IER)
          write(22,*)'     Log Likelihhod = ',TLIKE
	  write (22,*) '  --------------------------' 

	  IF (IER.LE.0) THEN
	   NS=NS+1
	   SLIKS(NS)=TLIKE
          ELSEIF (IER.GT.0) THEN
	 WRITE(22,*)'Due to error code ',IER,' this solution ignored'
	  ENDIF
c	  write (22,*) 'End of Random start ',I
	  TMP=0.0
	END DO
	DO I=1,NKMEANS
          write (*,*) '  Calling K-means Start ',I
          write (22,*) '  K-means Start ',I
          write (22,*) '     Seeds ', IX,IY,IZ
	  CALL KMEANS(NIND,NATT,NCOV,NG,X,XMU2,XVAR2,T2,IER)
	  CALL FIT(NIND,NATT,NCOV,X,NG,XMU,XVAR,T,XMU2,XVAR2,T2,
     &	           TOL,MAXIT,AIT,XLIKE,TLIKE,TMP,FACT,B,D,NATTQ,IER)
          write(22,*)'     Log Likelihhod = ',TLIKE
	  write (22,*) '  --------------------------' 
	  IF (IER.LE.0) THEN
	   NS=NS+1
	   SLIKS(NS)=TLIKE
          ELSEIF (IER.GT.0) THEN
	 WRITE(22,*)'Due to error code ',IER,' this solution ignored'
	  ENDIF
c	  write (22,*) 'End of K-means start ',I
	END DO
        write (22,*)
        write (22,*) 'Possible starting solutions found'
        write (22,*) (SLIKS(I),I=1,NS) 
	write (22,*) '  --------------------------' 
        write (22,*)
	RETURN
	END


      SUBROUTINE UPDSOL(NG,NATT,NCOV,XMU,XMU2,XVAR,XVAR2,T,
     &                  T2,B,D,BT,DT,FACT,NATTQ)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'
      INTEGER FACT(3)
	DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),XMU(MAXNG,MNATT),
     &            T(MAXNG)
        DIMENSION XVAR2(MAXNG,MNATT*(MNATT+1)/2),
     &            XMU2(MAXNG,MNATT),
     &            T2(MAXNG)
	DIMENSION D(MAXNG,MNATT),
     &            B(MAXNG,MNATT,MNATTQ)
	DIMENSION DT(MAXNG,MNATT),
     &            BT(MAXNG,MNATT,MNATTQ)
	IF (NCOV.EQ.3) THEN
	  DO I=1,NATT
	    XVAR(1,I)=XVAR2(1,I)
	  END DO
	ELSEIF (NCOV.EQ.4) THEN
	  DO K=1,NG
         DO I=1,NATT
	    XVAR(1,I)=XVAR2(1,I)
	   END DO
	  ENDDO
	ELSEIF (NCOV.EQ.1) THEN
         DO I=1,NATT*(NATT+1)/2
	     XVAR(1,I)=XVAR2(1,I)
	   END DO
	ELSEIF (NCOV.EQ.2) THEN
	  DO K=1,NG
         DO I=1,NATT*(NATT+1)/2
	     XVAR(K,I)=XVAR2(K,I)
	   END DO
	  ENDDO	 
      ENDIF	

	DO K=1,NG
	 T(K)=T2(K)
	 DO I=1,NATT
	  XMU(K,I)=XMU2(K,I)
       END DO
	END DO

	IF (FACT(1).EQ.1) THEN
	 DO 200 K=1,NG
	  DO 200 I=1,NATT
	   D(K,I)=DT(K,I)
	   DO 200 J=1,NATTQ
	    B(K,I,J)=BT(K,I,J)
200        CONTINUE 
        ENDIF
	RETURN
	END



	SUBROUTINE RANDCON(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2,
     &                   RDSB)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'

	REAL X(MNIND,MNATT)
	DIMENSION XMU2(MAXNG,MNATT),
     &            XVAR2(MAXNG,MNATT*(MNATT+1)/2),T2(MAXNG)
	IF ((NCOV.EQ.1).OR.(NCOV.EQ.2)) THEN
	  CALL RANDSTD(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2,RDSB)
	ELSEIF ((NCOV.EQ.3).OR.(NCOV.EQ.4)) THEN
	  CALL RANDDIAG(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2,RDSB)
        ENDIF
        RETURN
	END

      SUBROUTINE RANDSTD(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2,
     &                   RDSB)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'

        COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	REAL X(MNIND,MNATT)
	DIMENSION XMU2(MAXNG,MNATT),
     &            XVAR2(MAXNG,MNATT*(MNATT+1)/2),T2(MAXNG)
	DIMENSION N(MAXNG)
	DO 177 KK=1,NG
        N(KK)=0
        DO 177 II=1,NATT
          XMU2(KK,II)=0.0
            DO 176 JJ=1,II
	      IND=II*(II-1)/2+JJ
              XVAR2(KK,IND)=0.0
176         CONTINUE
177     CONTINUE
        IF (RDSB.EQ.100) THEN
         DO 111 I=1,NIND
          R=RANDNUM()
          JG=INT(R*FLOAT(NG))+1
          N(JG)=N(JG)+1
          DO 1091 JJ=1,NATT
            XMU2(JG,JJ)=XMU2(JG,JJ)+DBLE(X(I,JJ))	
            DO 1091 II=1,JJ
  	      IND=JJ*(JJ-1)/2+II
              XVAR2(JG,IND)=XVAR2(JG,IND)+DBLE(X(I,II)*X(I,JJ))
1091        CONTINUE
111      CONTINUE
         IPTS=NIND
        ELSE
         IPTS=INT(RDSB/100.0*FLOAT(NIND))
         DO 211 I=1,IPTS
           P=INT(RANDNUM()*FLOAT(NIND))+1
          R=RANDNUM()
          JG=INT(R*FLOAT(NG))+1
          N(JG)=N(JG)+1
          DO 2091 JJ=1,NATT
            XMU2(JG,JJ)=XMU2(JG,JJ)+DBLE(X(P,JJ))	
            DO 2091 II=1,JJ
		    IND=JJ*(JJ-1)/2+II
              XVAR2(JG,IND)=XVAR2(JG,IND)+DBLE(X(P,II)*X(P,JJ))
2091        CONTINUE
211      CONTINUE
        ENDIF
          DO 213 KK=1,NG
           T2(KK)=FLOAT(N(KK))/FLOAT(IPTS)
             DO 213 JJ=1,NATT
              XMU2(KK,JJ)=XMU2(KK,JJ)/(FLOAT(N(KK)))
              DO 213 II=1,JJ
	        IND=JJ*(JJ-1)/2+II
                XVAR2(KK,IND)=
     &          XVAR2(KK,IND)-(XMU2(KK,II)*XMU2(KK,JJ)*FLOAT(N(KK)))
                XVAR2(KK,IND)=XVAR2(KK,IND)/FLOAT(N(KK)-1)
212       CONTINUE	
213       CONTINUE	
	RETURN
	END


      SUBROUTINE RANDDIAG(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2,RDSB)
c         Have to modify 
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'

	COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	REAL X(MNIND,MNATT)
	DIMENSION XMU2(MAXNG,MNATT),
     &            XVAR2(MAXNG,MNATT*(MNATT+1)/2),
     &            T2(MAXNG)
	DIMENSION N(MAXNG)

	DO 177 KK=1,NG
        N(KK)=0
        DO 177 II=1,NATT
          XMU2(KK,II)=0.0
          XVAR2(KK,II)=0.0
177     CONTINUE
        IPTS=INT(RDSB/100.0*FLOAT(NIND))
        DO 111 I=1,IPTS
          P=INT(RANDNUM()*FLOAT(NIND))+1
          R=RANDNUM()
          JG=INT(R*FLOAT(NG))+1
          N(JG)=N(JG)+1
          DO 1091 JJ=1,NATT
            XMU2(JG,JJ)=XMU2(JG,JJ)+DBLE(X(P,JJ))	
            XVAR2(JG,JJ)=XVAR2(JG,JJ)+DBLE(X(P,JJ)*X(P,JJ))
1091      CONTINUE
111     CONTINUE

        DO 112 KK=1,NG
           T2(KK)=FLOAT(N(KK))/FLOAT(IPTS)
             DO 112 JJ=1,NATT
              XMU2(KK,JJ)=XMU2(KK,JJ)/(FLOAT(N(KK)))
                 XVAR2(KK,IND)=
     &              XVAR2(KK,JJ)-(XMU2(KK,JJ)*XMU2(KK,JJ))*FLOAT(N(KK))
                 XVAR2(KK,JJ)=XVAR2(KK,JJ)/FLOAT(N(KK)-1)
112     CONTINUE	
	RETURN
	END



	SUBROUTINE KMEANCON(NIND,NATT,NCOV,X,NG,XMU2,XVAR2,T2)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'

        COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	REAL X(MNIND,MNATT)
	DIMENSION XVAR2(MAXNG,MNATT*(MNATT+1)/2),
     &            XMU2(MAXNG,MNATT),T2(MAXNG)

	RETURN
	END

      SUBROUTINE FIT(NIND,NATT,NCOV,X,NG,XMU,XVAR,T,XMU2,XVAR2,T2,
     &               TOL,MAXIT,AIT,XLIKE,TLIKE,TMP,FACT,B,D,NATTQ,IER)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'
      INTEGER FACT(3)
	REAL X(MNIND,MNATT)
	DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &            XMU(MAXNG,MNATT),T(MAXNG)
        DIMENSION XVAR2(MAXNG,MNATT*(MNATT+1)/2),
     &            XMU2(MAXNG,MNATT),T2(MAXNG)
	DIMENSION D(MAXNG,MNATT),
     &            B(MAXNG,MNATT,MNATTQ)
	DIMENSION DT(MAXNG,MNATT),
     &            BT(MAXNG,MNATT,MNATTQ)
        MODE=1 
	SETUPON=1
	FORMOT=0
        write (*,*) 'Calling EM from Fit'
        CALL EM(NIND,NATT,NG,NCOV,X,XMU2,XVAR2,T2,TLIKE,
     &    TOL,MAXIT,MODE,FACT,BT,DT,NATTQ,FORMOT,SETUPON,IER)
        write (*,*) 'Return to fit from EM '
 

	IF (IER.GT.0) RETURN

	IF ((TLIKE.GT.XLIKE).OR.(TMP.EQ.1.0)) THEN
         write (*,*) 'Update'
	 XLIKE=TLIKE
	 CALL UPDSOL(NG,NATT,NCOV,XMU,XMU2,XVAR,XVAR2,T,T2,
     &	 B,D,BT,DT,FACT,NATTQ)
         TMP=0.0
	ENDIF
 

	RETURN
	END

      
      SUBROUTINE SETUP(NIND,NATT,NCOV,NG0,NG1,NKMEANS,NRANDS,RDSB,
     &  INFYLE,OFYLE,MFYLE,MODFYLE,OPTION,START,TOLS,SUB,SCAL,FACT,
     &  NATTQ,REV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      INTEGER SUB(2+MNATT),FACT(3)
      DOUBLE PRECISION TOLS(2)
      CHARACTER*255 INFYLE
      CHARACTER*255 OFYLE
      CHARACTER*255 MFYLE
      CHARACTER*255 MODFYLE
      
C	WRITE (*,*) 'Factor Analysis (0=No,1-Yes)'
C        READ (*,*) FACT(1) 
        FACT(1)=1
	IER=0
        WRITE(*,*)'      ---------------------------------' 
        WRITE(*,*)'      VERSION OF EMMIX TO FIT MIXTURES'
        WRITE(*,*)'            OF FACTORS ANALYZERS'
        WRITE(*,*)'            Jun 7th 2001 - Model file'
        WRITE(*,*)'            added for cross validation'
	WRITE(*,*)'            Jun 18th 2001 -Add common D'
        WRITE(*,*)'      ---------------------------------' 
C	CALL GETOPT(OPTION)
       OPTION=2.0
       CALL GETIO(INFYLE,OFYLE,IER)
	write (*,*) 'What is the name of the matlab file'
	write (*,*) 'you wish to output to (must end  in .m)'
        READ (*,'(A)') MFYLE
	write (*,*) 'What is the name of the model file'
	write (*,*) 'you wish to output to store the model '
        READ (*,'(A)') MODFYLE
	IF (IER.NE.0) RETURN
       CALL GETNIND(NIND,SUB)
       CALL GETNATT(NATT,SUB)
       WRITE (*,*) 'Do you wish to transpose the data'
       WRITE (*,*) ' so that you cluster the columns' 
       WRITE (*,*) ' (0=No ,1=Yes)'
       READ (*,*) REV
       WRITE (*,*) 'How many factor dimensions'
       READ (*,*) NATTQ
       WRITE (*,*) 'Scale Data (0-No,1-Yes)'
       READ (*,*) SCAL 
	CALL GETNG(NG0,NG1,OPTION)
	CALL GETNCOV(NCOV)
        CALL GETST(START)
	IF (START.EQ.3) THEN
 	 CALL GETSRCH(NRANDS,NKMEANS,RDSB)
        ENDIF
	write (*,*) 'Use 0-Unrestricted Ds, 1-Common Ds'
	READ (*,*) FACT(2)
	write (*,*) 'Intialization method for B'
	write (*,*) ' 0-Bishop (based on eigenvalues/vectors)'
        write (*,*) ' 1-Hinton (random)'
	READ (*,*) FACT(3)

C	write (*,*) 'What form do you want the output'
C	write (*,*) ' (0-Parameters, 1-Partition, 3-Both)'
C	read (*,*) FORMOT
	
	CALL GETSEEDS(IX,IY,IZ)
	
	
	      
c	OPTION=2
c	INFYLE='green'
c	OFYLE='green.out'
c	NIND=50
c	SUB(1)=100
c	NATT=4
c	SUB(2)=4
c	NG0=2
c	NG1=2
c	NCOV=2
c	START=3
c	NRANDS=1
c	RDSB=100
c	NKMEANS=0
c	FACT=1
c	NATTQ=2
c	SCAL=0
c     IX=1
c	IY=2
c	IZ=3
	 
      RETURN
      END

      SUBROUTINE GETNIND(NIND,SUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'
      INTEGER SUB(2+MNATT)
      SUB(1)=100
10    WRITE(*,*) 'Number of data points: '
      READ (*,*) NIND
C      WRITE (*,*) 'Percentage of data to be used'
C      READ (*,*) SUB(1) 
       SUB(1)=100
      IF (NIND*SUB(1).GT.MNIND*100) THEN
       WRITE(*,*) 'You may need to Change EMMIX-f1.max' 
      ENDIF
      RETURN
      END
 
      SUBROUTINE GETNATT(NATT,SUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX-f1.max'
      INTEGER SUB(2+MNATT)
10    SUB(2)=0
      WRITE(*,*) 'Total number of variables/dimensions'
      READ (*,*) NATT
      WRITE(*,*)'Use a subset of the data variables '
      WRITE(*,*)'(New number of varables):'

      READ (*,*) SUB(2) 
      IF (SUB(2).GT.MNATT) THEN
       WRITE (*,*) 'You may need to CHANGE EMMIX-f1.max'
      ENDIF
      IF (SUB(2).NE.NATT) THEN
        WRITE (*,*) 'Enter variable column numbers on a single line'
        READ (*,*) (SUB(2+I),I=1,SUB(2))
      ENDIF
      RETURN
      END

    
      SUBROUTINE GETIO(INFYLE,OFYLE,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*255 INFYLE
      CHARACTER*255 OFYLE

C    Read in input file and open
      WRITE (*,*)'Enter name of input file: '
      READ (*,'(A)') INFYLE
      OPEN (UNIT=21,FILE=INFYLE,STATUS = 'OLD',ERR=505)
      GOTO 506

505   WRITE (*,*) 'Cannot locate ',INFYLE
      WRITE (*,*) 'please re-enter file name: '
      READ (*,'(A)') INFYLE
      OPEN (UNIT=21,FILE=INFYLE,STATUS = 'OLD',ERR=505)

506   CONTINUE
      CLOSE(21)

C    Read in output file and open^M
      WRITE (*,*)'Enter name of output file: '
      READ (*,'(A)') OFYLE
c     OPEN (UNIT=22,FILE=OFYLE,STATUS = 'UNKNOWN')
      
      RETURN
      END

      SUBROUTINE GETOPT(OPTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      WRITE(*,519)
519   FORMAT(
     &/,12X,' Do you wish to:',
     &/,12X,'  1. Fit a g-component normal mixture model for a',
     &/,12X,'     specified g',
     &/,12X,'  2. Fit a g-component normal mixture model for a',
     &/,12X,'     range of values of g',
     &/,12X,'  3. Perform discriminant analysis')

      READ (*,*) TEMP
      IF (TEMP.EQ.1) THEN
       OPTION=2.0
      ELSE 
       OPTION=3.0
      ENDIF
      RETURN
      END

      SUBROUTINE GETNG(NG0,NG1,OPTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (OPTION.EQ.3) THEN
502    WRITE(*,*)'What is the minimum number of groups you want to fit'
       READ (*,*) NG0
       WRITE(*,*) 'What is the maximum number of groups you want to fit'
       READ (*,*) NG1
       IF (NG0.GT.NG1) THEN
         WRITE (*,*) '   Min must be less than max try again'
         GOTO 502
       ENDIF
      ELSE
       WRITE(*,*) 'How many groups do you want to fit: '
       READ (*,*) NG0
       NG1=NG0
      ENDIF
      RETURN
      END
      

      SUBROUTINE GETST(START)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
520   WRITE(*,*) 'Switch for initialisation'
      WRITE(*,*)'    1 = initial partition of data,'
      WRITE(*,*)'    2 = initial parameter estimates,'
      WRITE(*,*)'    3 = automatic initial grouping '
      READ (*,*) START
      IF ((START.NE.1).AND.(START.NE.2).AND.(START.NE.4)
     &                              .AND.(START.NE.3)) THEN
        WRITE (*,*) '  ERROR expecting a 1,2 or a 3 repeat answer: '
        GOTO 520
      ENDIF
      RETURN
      END

      SUBROUTINE GETNCOV(NCOV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
510   WRITE(*,*)'Covariance matrix option (1 = equal,2 = unrestricted,'
      WRITE(*,*)'   3 = diagonal equal,4 = diagonal unrestricted): '
      READ (*,*) NCOV
      IF ((NCOV.LT.1).AND.(NCOV.GT.4)) THEN
        WRITE (*,*) '  ERROR expecting a 1,2,3 or 4 repeat answer: '
        GOTO 510
      ENDIF
      RETURN
      END

	SUBROUTINE GETSRCH(NRANDS,NKMEANS,RDSB)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	WRITE (*,*) 'How many Random starts'
	READ (*,*) NRANDS
	WRITE (*,*) 'What percentage of the data to be used (eg 10=10%)'
	READ (*,*) RDSB
	WRITE (*,*) 'How many K-means starts'
	READ (*,*) NKMEANS
	RETURN
	END

	SUBROUTINE GETSEEDS(IX,IY,IZ)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	WRITE (*,*) 'Enter 0 to auto seed'
	WRITE (*,*) 'Random Seed 1'
	READ (*,*) IX
	IF (IX.NE.0) THEN
	  WRITE (*,*) 'Random Seed 2'
	  READ (*,*) IY
	  WRITE (*,*) 'Random Seed 3'
	  READ (*,*) IZ
      ELSE
	  WRITE (*,*) 'Not done yet'
	ENDIF
	RETURN
	END
C
      SUBROUTINE EM(NIND,NATT,NG,NCOV,X,XMU,XVAR,T,
     &       TLIKE,TOL,MAXIT,MODE,FACT,B,D,NATTQ,FORMOT,SETUPON,IER)
C      The purpose of the subroutine is to implement the EM algorithm
C      of Dempster et al. (1977).

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'         
      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG,        ! Number of components or groups
     &        NCOV,      ! Covariance structure
     &        MODE,      ! 1-(silent) no output 2-output
     &        FACT(3),   !control for mix of Fact
     &        IDT(MNIND),N(MAXNG)
      REAL X(MNIND,MNATT)! Data set/sample
C     
      DIMENSION XMU(MAXNG,MNATT),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          T(MAXNG),
     &          TXML(MITER) ,XLA(MITER)
     
      DIMENSION W(MNIND,MAXNG),DV(MAXNG)

      DIMENSION D(MAXNG,MNATT),
     &          B(MAXNG,MNATT,MNATTQ)

      DIMENSION GAMM(MAXNG,MNATT,MNATTQ)
 
      IER=0
      IOUNT=0
      IF (FACT(1).EQ.1) THEN
      write (*,*) 'Calling Initial Factor Setup'
      CALL FSASET2(NIND,NATT,NG,NCOV,X,DV,XVAR,V,
     &             D,B,XMU,NATTQ,MODE,FACT,SETUPON,T,IER)
      ENDIF
 
      write (*,*) '   Applying EM'
      IF (MODE.EQ.2) THEN
       WRITE(22,*) '   Progress at Each Iteration'
       WRITE(22,*) '   --------------------------' 
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    MAIN ITERATIVE LOOP ENTRY POINT
900   CONTINUE     
      
      IOUNT=IOUNT+1
      IF (FACT(1).EQ.1) THEN
       CALL ESTEP(NIND,NATT,NCOV,NG,X,XMU,V,DV,T,W,
     &            XLOGL,IER)
       CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,XMU,T,W,FACT,IER)
      ELSE
       CALL ESTEP(NIND,NATT,NCOV,NG,X,XMU,V,DV,T,W,
     &            XLOGL,IER)
      ENDIF
      IF (IER.GT.0) RETURN 
      TXML(IOUNT)=XLOGL 
      write (*,*) 'Log-Likelihood at iteration',IOUNT,' = ',TXML(IOUNT)
      write (*,*) ' Determinants=',(DV(KKK),KKK=1,NG)
      IF (MODE.EQ.2) THEN
      write (22,*) 'Log-Likelihood at iteration',IOUNT,' = ',TXML(IOUNT)
      write (22,*) ' Determinants=',(DV(KKK),KKK=1,NG)
      ENDIF
      IF (XIT(IOUNT,MAXIT,TOL,TXML,XLA).EQ.1) GOTO 1099
      IF (FACT(1).EQ.1) THEN
       CALL FSAMSP(NIND,NATT,NG,NCOV,XVAR,V,D,B,
     &                XMU,GAMM,NATTQ,DV,FACT,T,IER)
      ELSE
       CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,XMU,T,W,FACT,IER)
      ENDIF
      IF (IER.GT.0) GOTO 1099
      
      GOTO 900
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

1099  CONTINUE
      TLIKE=TXML(IOUNT)
       CALL ESTEP(NIND,NATT,NCOV,NG,X,XMU,V,DV,T,W,
     &            XLOGL,IER)
C        OUTPUT



         IF (FORMOT.GT.0) THEN
3012      FORMAT('T=[',1000G15.5)
          write (23,3012) (T(K),K=1,NG)
          write (23,*) '];'
	  IF (FACT(1).EQ.1) THEN
	  DO 3011 K=1,NG
3013        FORMAT('MU(',I4,',:)=[',1000G15.5)
	    write (23,3013) K,(XMU(K,I),I=1,NATT)
            write (23,*) '];'
3011      CONTINUE
	   DO 1012 K=1,NG
	    DO 1011 I=1,NATT
3014          FORMAT('B(',I4,',:,',I4,')=[',1000G15.5)
              write(23,3014)I,K,(B(K,I,J),J=1,NATTQ)
              write (23,*) '];'
1011        CONTINUE
1012       CONTINUE
	   DO 2012 K=1,NG
	    DO 2011 I=1,NATT
3015         FORMAT('GAMM(',I4,',:,',I4,')=[',1000G15.5)
             write(23,3015)I,K,(GAMM(K,I,J),J=1,NATTQ)
             write (23,*) '];'
2011        CONTINUE
2012       CONTINUE
	    DO 1013 K=1,NG
            DO 1013 J=1,NATT
3016         FORMAT('D(',I4,',',I4,')=[',1000G15.5)
	    write(23,3016) K,J,D(K,J)
              write (23,*) '];'
1013        CONTINUE
	  ENDIF

C      Calculate multivariate density of each point for every group
	DO 949 K=1,NG
	 N(K)=0
949     CONTINUE
        DO 950 JJ=1,NIND
	  MAXW=1
          DO 920 K=2,NG
	   IF (W(JJ,K).GT.W(JJ,MAXW)) THEN
	    MAXW=K
	   ENDIF
920       CONTINUE
        write (23,*) 'idt(',JJ,',1)=',MAXW,';'
        IDT(JJ)=MAXW
	N(MAXW)=N(MAXW)+1
	DO 940 I=1,NG
	 write (23,*) '   idt(',JJ,',',I+1,')=',W(JJ,I),';'
940     CONTINUE
950     CONTINUE	     
       write (22,*) 'Implied Grouping' 
       write (22,201) (IDT(I),I=1,NIND)
201    FORMAT (2X,10I4)
       write (22,*) 'Number in each group'
       write (22,202) (N(K),K=1,NG)
202    FORMAT (2X,10I4) 
	   ENDIF



      RETURN
      END


      SUBROUTINE ESTEP(NIND,NATT,NCOV,NG,X,XMU,V,DV,T,W,
     &            XLOGL,IER)
C     This Subroutine implements the E-step of the EM algorithm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'         
      
      COMMON /STORE2/ FLAGS,FYLENO

C   INPUT PARAMETERS
      INTEGER   FLAGS(40),     ! Flags summarising options chosen (see below)
     &          FYLENO         ! File id for output file
	INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG,        ! Number of components or groups
     &        NCOV       ! Covariance structure
      REAL X(MNIND,MNATT)  ! Data set/sample

      DIMENSION XMU(MAXNG,MNATT),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          T(MAXNG),
     &          TXML(MITER)

      DIMENSION W(MNIND,MAXNG),DV(MAXNG)
      DIMENSION AL(MAXNG)
      IER=0
      PI=3.141592653
      XLOGL=0.0

C      Calculate multivariate density of each point for every group
       DO 950 JJ=1,NIND
        GUM=0.0
        DO 920 K=1,NG
	    IF (DV(K).EQ.0.0) THEN
	      IER=2
	      RETURN
	    ENDIF
          AL(K)=0.0
          DO 910 I=1,NATT
            XJJI=DBLE(X(JJ,I))
            ALTEMP=(XJJI-XMU(K,I))
CCC            DO 910 J=1,NATT
            DO 905 J=1,I-1
              XJJJ=DBLE(X(JJ,J))
CCC              AL(K)=AL(K)+ALTEMP*V(K,IC(I,J))*
CCC     &              (XJJJ-XMU(K,J))
              AL(K)=AL(K)+2.0*ALTEMP*V(K,IC(I,J))*
     &              (XJJJ-XMU(K,J))
905         CONTINUE
	    AL(K)=AL(K)+ALTEMP*ALTEMP*V(K,IC(I,I))
910       CONTINUE
C         Check as if AL(K) too large under flow may occur
          IF (AL(K).GT.DENMAX) THEN
            AL(K)=0.0
          ELSE
C           Calculate component density
            AL(K)=(-0.5)*AL(K)
            AL(K)=EXP(AL(K))/(SQRT(DV(K))*(2.*PI)**(NATT/2.0))
          ENDIF
C         Check to stop underflow
          IF ((T(K).GT.XLOWEM).AND.(AL(K).GT.XLOWEM)) THEN
C           Calculate mixture density
            GUM=GUM+T(K)*AL(K)
          ELSE
          ENDIF
920     CONTINUE


        IF (GUM.EQ.0.0) THEN
          DO 930 K=1,NG
930         W(JJ,K)=0.0
            IER=-111
        ELSE
          DO 940 K=1,NG
C           Check to catch numerical underflow
            IF (T(K).LT.XLOWEM.OR.AL(K).LT.XLOWEM) THEN
              W(JJ,K)=0.0
            ELSE
C             Calculate posterior probabilities
              W(JJ,K)=T(K)*AL(K)/GUM
            ENDIF
940       CONTINUE
C         Calculate Log-Likelihood contribution from point JJ
          XLOGL=XLOGL+LOG(GUM)
       ENDIF
950   CONTINUE
CC      write (*,*) 'Log-like=',XLOGL
CC      write (*,*) 'Det=',DV(1),DV(2)
      RETURN
      END

      FUNCTION IC(I,J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (I.LE.J) THEN
        IC=J*(J-1)/2+I
      ELSE
        IC=I*(I-1)/2+J
      ENDIF
      RETURN
      END


      SUBROUTINE MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,
     &           XMU,T,W,FACT,IER)
C     This Subroutine implements the M-step of the EM algorithm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'         
      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG,        ! Number of components or groups
     &        NCOV,      ! Covariance structure
     &        FACT(3)    ! Control Factor Analysers
      REAL X(MNIND,MNATT)  ! Data set/sample

      DIMENSION XMU(MAXNG,MNATT),
     &          XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          T(MAXNG),
     &          TXML(MITER)

      DIMENSION W(MNIND,MAXNG),DV(MAXNG)
         IER=0
C     Compute new estimates of mixing proportions 
         CALL CALC_T(NIND,NATT,NG,W,T,IER)
         IF (IER.NE.0) RETURN
C     Compute new estimates of group means
         CALL CALC_XMU(NIND,NATT,NG,X,W,XMU,T)
         IF (IER.NE.0) RETURN
C     Compute new estimates of covariance matrices
         CALL CALC_XVAR(NIND,NATT,NG,NCOV,X,W,XMU,T,XVAR)

      IER=0
	IF (FACT(1).EQ.0) THEN
      CALL GDET(NCOV,NATT,NG,XVAR,V,DV,IER)
       IF (IER.NE.0) THEN
        write (*,*) '  ERROR: Problem Inverting Covariance Matrix'
	RETURN
       ENDIF
	ENDIF
         RETURN
         END       
      
      SUBROUTINE CALC_XMU(NIND,NATT,NG,X,W,XMU,T)
C     This Subroutine calculates estimates of the group means 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG         ! Number of components or groups
      DIMENSION T(MAXNG),XMU(MAXNG,MNATT)
      DIMENSION W(MNIND,MAXNG)
      REAL X(MNIND,MNATT)
C     Compute new estimates of group means (XMU)
      DO 1310 K=1,NG
        DO 1310 J=1,NATT
          XMU(K,J)=0.0
          DO 1300 JJ=1,NIND
           XMU(K,J)=XMU(K,J)+DBLE(X(JJ,J))*W(JJ,K)
1300      CONTINUE
          XMU(K,J)=XMU(K,J)/(T(K)*FLOAT(NIND))
1310  CONTINUE
      RETURN
      END

      SUBROUTINE CALC_XVAR(NIND,NATT,NG,NCOV,X,W,XMU,T,XVAR)
C     This Subroutine calculates estimates of the group covariances 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG         ! Number of components or groups
      DIMENSION T(MAXNG),XMU(MAXNG,MNATT)
      DIMENSION W(MNIND,MAXNG),
     &          XVAR(MAXNG,MNATT*(MNATT+1)/2)
      REAL X(MNIND,MNATT)
C     Compute new estimate of covariance matrix for each group

       IF (NCOV.EQ.4) THEN
        DO 60 K=1,NG
          DO 20 I=1,NATT
              XVAR(K,IC(I,I))=0.0
20        CONTINUE
          DO 40 JJ=1,NIND
           DO 30 I=1,NATT
       XVAR(K,IC(I,I))=XVAR(K,IC(I,I))+(DBLE(X(JJ,I))-XMU(K,I))**2
     &                  *W(JJ,K)/(T(K)*FLOAT(NIND))
30          CONTINUE
40        CONTINUE
60      CONTINUE

       ELSEIF (NCOV.EQ.3) THEN
        DO 160 K=1,NG
          DO 120 I=1,NATT
              XVAR(K,IC(I,I))=0.0
120       CONTINUE
          DO 140 JJ=1,NIND
           DO 130 I=1,NATT
      XVAR(K,IC(I,I))=XVAR(K,IC(I,I))+(DBLE(X(JJ,I))-XMU(K,I))**2
     &                  *W(JJ,K)/(T(K)*FLOAT(NIND))
130         CONTINUE
140       CONTINUE
160     CONTINUE
        DO 180 I=1,NATT
        XVAR(1,IC(I,I))=XVAR(1,IC(I,I))*FLOAT(NIND)
          DO 170 K=2,NG
        XVAR(1,IC(I,I))=XVAR(1,IC(I,I))+XVAR(K,IC(I,I))*FLOAT(NIND)
170       CONTINUE
        XVAR(1,IC(I,I))=XVAR(1,IC(I,I))/FLOAT(NIND)
180     CONTINUE
         

       ELSEIF (NCOV.EQ.2) THEN

        DO 260 K=1,NG
          DO 220 J=1,NATT
            DO 220 I=1,J
              XVAR(K,IC(I,J))=0.0
220       CONTINUE
          DO 240 JJ=1,NIND
            DO 230 J=1,NATT
             DO 230 I=1,J
        XVAR(K,IC(I,J))=XVAR(K,IC(I,J))+(DBLE(X(JJ,I))-XMU(K,I))
     &        *(DBLE(X(JJ,J))-XMU(K,J))*W(JJ,K)/(T(K)*FLOAT(NIND))
230         CONTINUE
240       CONTINUE
260     CONTINUE

       ELSEIF (NCOV.EQ.1) THEN

        DO 360 K=1,NG
          DO 320 J=1,NATT
            DO 320 I=1,J
              XVAR(K,IC(I,J))=0.0
320       CONTINUE
          DO 340 JJ=1,NIND
            DO 330 J=1,NATT
              DO 330 I=1,J
         XVAR(K,IC(I,J))=XVAR(K,IC(I,J))+(DBLE(X(JJ,I))-XMU(K,I))
     &             *(DBLE(X(JJ,J))-XMU(K,J))*W(JJ,K)
330         CONTINUE
340       CONTINUE
360     CONTINUE

        DO 390 J=1,NATT
          DO 390 I=1,J
            XVAR(1,IC(I,J))=XVAR(1,IC(I,J))*FLOAT(NIND)
            DO 380 K=2,NG
       XVAR(1,IC(I,J))=XVAR(1,IC(I,J))+XVAR(K,IC(I,J))*FLOAT(NIND)
380         CONTINUE
        XVAR(1,IC(I,J))=XVAR(1,IC(I,J))/NIND
390     CONTINUE

      ENDIF    

      RETURN
      END 



      SUBROUTINE CALC_T(NIND,NATT,NG,W,T,IER)
C     This Subroutine calculates estimates of the mixing proportions
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG         ! Number of components or groups
      DIMENSION T(MAXNG)
      DIMENSION W(MNIND,MAXNG),WSUM(MAXNG)

      DO 1370 K=1,NG
        WSUM(K)=0
        DO 1369 J=1,NIND
          WSUM(K)=WSUM(K)+W(J,K)
1369    CONTINUE
1370  CONTINUE

      DO 970 K=1,NG
         IF (WSUM(K).EQ.0) THEN
          IER=22
        WRITE(22,*)'  Problem: No points allocated to component ',K
        WRITE(22,*)'           during EM iteration '
        RETURN
         ENDIF
C       Calculate mixing proportion for group K
        T(K)=WSUM(K)/NIND
970   CONTINUE
      RETURN
      END



C
      FUNCTION XIT(IOUNT,MAXIT,TOL,TXML,XLA)
C     Test for exit from EM Algorithm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)       
      INCLUDE 'EMMIX-f1.max'  
      INTEGER FLAGS(40),FYLENO,USA(MNIND)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION TXML(MITER),XLA(MITER)
      XIT=0.0
      LAST=IOUNT-10
      IF (IOUNT.LE.10) RETURN


      IF (IOUNT.GE.MAXIT) THEN      
         WRITE (22,115) MAXIT 
115      FORMAT (/2X,'Note: This sample did not converge in ',I3, 
     &     ' iterations.',/8X,'However the program will continue ', 
     &     'to print results ',/8X,'obtained from the last cycle ', 
     &     'estimates.') 
         XIT=1.0
      ENDIF 
      ALIM=TOL*TXML(LAST)
      DIFF=TXML(IOUNT)-TXML(LAST)
      IF (ABS(DIFF).LE.ABS(ALIM)) THEN
        XLA(IOUNT)=0
        XIT=1.0 
      ENDIF

C      IF (FLAGS(21).LT.0) THEN
C        XLAOLD=XLA(LAST)
C        XNUM=TXML(IOUNT)-TXML(IOUNT-1)
C        DEM=TXML(IOUNT-1)-TXML(IOUNT-2)
C        IF (DEM.LT.1E-35) THEN
C          XLANEW=0
CC        ELSE
C          C=XNUM/DEM
C        ENDIF
C        IF ((C.LT.(.99)).OR.(C.GT.(1.01))) THEN
C         XLANEW=TXML(IOUNT-1)+XNUM*1/(1-C)
C        ELSE
C         XLANEW=0
C        ENDIF
C
C        XLA(IOUNT)=XLANEW

C        IF (XLA(IOUNT).NE.0) THEN
C          ALIM=TOL*XLAOLD
C          DIFF=XLA(IOUNT)-XLAOLD
C        ENDIF
C      IF ((ABS(DIFF).LE.ABS(ALIM)).AND.(XLA(IOUNT).GE.TXML(IOUNT)))
C     &      THEN
C          TXML(IOUNT)=XLA(IOUNT)
C          XIT=1.0 
C         ENDIF
C      ENDIF
      RETURN
      END

 
      SUBROUTINE GDET(NCOV,NATT,NG,XVAR,V,DV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          DV(MAXNG)

      IER=0
	IF (NCOV.EQ.1) CALL GDET1(NATT,NG,XVAR,V,DV,IER)
	IF (NCOV.EQ.2) CALL GDET2(NATT,NG,XVAR,V,DV,IER)
	IF (NCOV.EQ.3) CALL GDET3(NATT,NG,XVAR,V,DV,IER)
	IF (NCOV.EQ.4) CALL GDET4(NATT,NG,XVAR,V,DV,IER)
      
	IF (IER.NE.0) THEN
895   FORMAT (/2X,'Terminal error in matrix inversion for group ',I3,
     &        ': error code ',I3)
899   CONTINUE
c      WRITE (FYLENO,895) K, IER
c      WRITE (*,895)  K, IER
	ENDIF

      RETURN
      END


      SUBROUTINE GDET1(NATT,NG,XVAR,V,DV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.

      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          DV(MAXNG) 
        NNULL=0
        IT=0
        TOL=0.0
        DO 810 I=1,NATT
	    IND=I*(I+1)/2
          TOL=TOL+SQRT(XVAR(1,IND))
810     CONTINUE
        TOL=(TOL/NATT)*0.000001
        CALL SYMINV (1,XVAR,NATT,V,NULL,IER,RMAX,DV,TOL,NG)
        IF (IER.GT.0) RETURN
        IF (NULL.NE.0) THEN
c          WRITE (FYLENO,815) NULL
815       FORMAT (/2X,'Rank deficiency of common covariance matrix ',
     &                ' is ',I3)         
          NNULL=NULL+1
        ENDIF
        IF (DV(1).EQ.0) THEN
c         WRITE(FYLENO,*)'  Determinant is equal to zero for common cov.'
          IER=5
          RETURN
        ENDIF
      NULL=NNULL
      RETURN
	END

      SUBROUTINE GDET2(NATT,NG,XVAR,V,DV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          DV(MAXNG)	 
      NNULL=0
      DO 830 K=1,NG
        IT=0
        TOL=0.0
        DO 810 I=1,NATT
	    IND=I*(I+1)/2
          TOL=TOL+SQRT(XVAR(K,IND))
810     CONTINUE
        TOL=(TOL/NATT)*0.000001
        CALL SYMINV (K,XVAR,NATT,V,NULL,IER,RMAX,DV,TOL,NG)
        IF (IER.GT.0) RETURN
        IF (NULL.NE.0) THEN
c          WRITE (FYLENO,815)  K, NULL
815       FORMAT (/2X,'Rank deficiency of covariance matrix ',I3,
     &                ' is ',I3)         
          NNULL=NULL+1
        ENDIF
        IF (DV(K).EQ.0) THEN
c         WRITE (FYLENO,*)'  Determinant is equal to zero for grp ',K
          IER=5
          RETURN
        ENDIF
830   CONTINUE
      NULL=NNULL
      RETURN
	END

      SUBROUTINE GDET3(NATT,NG,XVAR,V,DV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          DV(MAXNG)
	 DV(1)=1.0
	 DO 830 I=1,NATT
	  IF (XVAR(1,I).LT.XMACHMAX) THEN
	    V(1,I)=1.0/XVAR(1,I)
	    DV(1)=DV(1)*XVAR(1,I)
        ELSE
	    IER=5
	    RETURN
	  ENDIF
830   CONTINUE
      RETURN
	END


      SUBROUTINE GDET4(NATT,NG,XVAR,V,DV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          DV(MAXNG)	 
      DO 830 K=1,NG
	 DV(K)=1.0
	 DO 830 I=1,NATT
	  IF (XVAR(K,I).LT.XMACHMAX) THEN
	    V(K,I*(I+1)/2)=1.0/XVAR(K,I*(I+1)/2)
	    DV(K)=DV(K)*XVAR(K,I*(I+1)/2)
        ELSE
	    IER=5
	    RETURN
	  ENDIF
830   CONTINUE
      RETURN
	END

      SUBROUTINE SYMINV(KNG,XVAR,N,V,NULLTY,IFAULT,RMAX,DV,
     &                  TOL,NG)
C       Modified from
C       Algorithm AS7, Applied Statistics, Vol.17, 1968.
C
C       ARGUMENTS:-
C       XVAR()     = input, the symmetric matrix to be inverted, stored in
C                 lower triangular form
C       N       = input, order of the matrix
C       V()     = output, the inverse of a (a generalized inverse if c is
C                 singular), also stored in lower triangular.
C                 c and a may occupy the same locations.
C       W()     = workspace, dimension at least n.
C       NULLTY  = output, the rank deficiency of a.
C       IFAULT  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 0 otherwise
C       RMAX    = output, approximate bound on the accuracy of the diagonal
C                 elements of c.  e.g. if rmax = 1.e-04 then the diagonal
C                 elements of c will be accurate to about 4 dec. digits.
C
C       Latest revision - 18 April 1981

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          DV(MAXNG)
      DIMENSION W(MNATT)

        NROW=N
        IFAULT=1
        IF(NROW.LE.0) GO TO 100
        IFAULT=0
C       Cholesky factorization of A, result in C
        CALL CHOLA(KNG,XVAR,NROW,V,NULLTY,IFAULT,RMAX,W,TOL)
        IF(IFAULT.NE.0) GO TO 100
C       Invert C & form the product (CINV)'*CINV, where CINV is the inverse
C       of C, row by row starting with the last row.
C       IROW = the row number, NDIAG = location of last element in the row.
        NN=NROW*(NROW+1)/2
cc        DO 200 IMY=1,NN
cc200       CMY(IMY)=V(KNG,IMY)
          DV(KNG)=1.0
		DO 200 J=1,NROW
	      JJ=J*(J+1)/2
	      DV(KNG)=DV(KNG)*V(KNG,JJ)*V(KNG,JJ)
200       CONTINUE
        IROW=NROW
        NDIAG=NN
10      IF(V(KNG,NDIAG).EQ.0.0) GO TO 60
        L=NDIAG
        DO 20 I=IROW,NROW
          W(I)=V(KNG,L)
          L=L+I
20      CONTINUE
        ICOL=NROW
        JCOL=NN
        MDIAG=NN
30      L=JCOL
        X=0.0
        IF(ICOL.EQ.IROW) X=1.0/W(IROW)
        K=NROW
40      IF(K.EQ.IROW) GO TO 50
        X=X-W(K)*V(KNG,L)
        K=K-1
        L=L-1
        IF(L.GT.MDIAG) L=L-K+1
        GO TO 40
50      V(KNG,L)=X/W(IROW)
        IF(ICOL.EQ.IROW) GO TO 80
        MDIAG=MDIAG-ICOL
        ICOL=ICOL-1
        JCOL=JCOL-1
        GO TO 30
60      L=NDIAG
        DO 70 J=IROW,NROW
          V(KNG,L)=0.0
          L=L+J
70      CONTINUE
80      NDIAG=NDIAG-IROW
        IROW=IROW-1
        IF(IROW.NE.0) GO TO 10
100     CONTINUE
        RETURN
        END


        SUBROUTINE CHOLA(KNG,XVAR,N,V,NULLTY,IFAULT,
     &                   RMAX,R,TOL)
C       Modified from
C       Algorithm AS6, Applied Statistics, Vol.17, 1968, with
C       modifications by A.J.Miller
C
C       ARGUMENTS:-
C       XVAR()     = input, a +ve definite matrix stored in lower-triangular
C                 form.
C       N       = input, the order of A
C       V()     = output, a lower triangular matrix such that U*U' = A.
C                 A & U may occupy the same locations.
C       NULLTY  = output, the rank deficiency of A.
C       IFAULT  = output, error indicator
C                       = 1 if N < 1
C                       = 2 if A is not +ve semi-definite
C                       = 0 otherwise
C       RMAX    = output, an estimate of the relative accuracy of the
C                 diagonal elements of U.
C       R()     = output, array containing bounds on the relative accuracy
C                 of each diagonal element of U.
C
C       Latest revision - 18 April 1981
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          V(MAXNG,MNATT*(MNATT+1)/2),
     &          R(MNATT)

C       ETA should be set equal to the smallest +ve value such that
C       1.0 + eta is calculated as being greater than 1.0 in the accuracy
C       being used.
C        DATA ETA/1.e-07/
        ETA=TOL
        IFAULT=1
        IF(N.LE.0) GO TO 100
        IFAULT=2
        NULLTY=0
        RMAX=ETA
        R(1)=ETA
        J=1
        K=0
C       Factorize column by column, ICOL = Column No.
        DO 80 ICOL=1,N
          L=0
C       IROW = row number within column ICOL
          DO 40 IROW=1,ICOL
            K=K+1
            W=XVAR(KNG,K)
            IF(IROW.EQ.ICOL) RSQ=(W*ETA)**2
            M=J
            DO 10 I=1,IROW
              L=L+1
              IF(I.EQ.IROW) GO TO 20
              W=W-V(KNG,L)*V(KNG,M)
              IF(IROW.EQ.ICOL) RSQ=RSQ+(V(KNG,L)**2*R(I))**2
              M=M+1
10          CONTINUE
20          IF(IROW.EQ.ICOL) GO TO 50
            IF(V(KNG,L).EQ.0.0) GO TO 30
            V(KNG,K)=W/V(KNG,L)
            GO TO 40
30          V(KNG,K)=0.0
            IF(ABS(W).GT.ABS(RMAX*XVAR(KNG,K))) GO TO 100
40        CONTINUE
C         End of row, estimate relative accuracy of diagonal element.
50        RSQ=SQRT(RSQ)
          IF(ABS(W).LE.5.*RSQ) GO TO 60
          IF(W.LT.0.0) GO TO 100
          V(KNG,K)=SQRT(W)
          R(I)=RSQ/W
          IF(R(I).GT.RMAX) RMAX=R(I)
          GO TO 70
60        V(KNG,K)=0.0
          NULLTY=NULLTY+1
70        J=J+ICOL
80      CONTINUE
        IFAULT=0.0
100     RETURN
        END


      SUBROUTINE GDETQ(NATT,XVAR,V,DV,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
      INCLUDE 'EMMIX-f1.max'
      DIMENSION XVAR(MNATTQ*(MNATTQ+1)/2),
     &          V(MNATTQ*(MNATTQ+1)/2)
      NNULL=0
	IT=0
	TOL=0.0
	DO 810 I=1,NATT
	  TOL=TOL+SQRT(XVAR(IC(I,I)))
810     CONTINUE
	TOL=(TOL/NATT)*0.000001
	CALL SYMINVQ(XVAR,NATT,V,NULL,IER,RMAX,DV,TOL,NG)
	IF (IER.GT.0) RETURN
	IF (NULL.NE.0) THEN
	  write (*,*) 'NULL=',NULL
815       FORMAT (/2X,'Rank deficiency of covariance matrix ',I3,
     &                ' is ',I3)
	  NNULL=NULL+1
        ENDIF
        IF (DV.EQ.0) THEN
	 write (*,*) 'DV=',DV
	 IER=5
	 RETURN
	ENDIF
830   CONTINUE
      NULL=NNULL
      RETURN
      END

      SUBROUTINE SYMINVQ(XVAR,N,V,NULLTY,IFAULT,RMAX,DV,
     &                  TOL,NG)
C       Modified from
C       Algorithm AS7, Applied Statistics, Vol.17, 1968.
C
C       ARGUMENTS:-
C       XVAR()     = input, the symmetric matrix to be inverted, stored in
C                 lower triangular form
C       N       = input, order of the matrix
C       V()     = output, the inverse of a (a generalized inverse if c is
C                 singular), also stored in lower triangular.
C                 c and a may occupy the same locations.
C       W()     = workspace, dimension at least n.
C       NULLTY  = output, the rank deficiency of a.
C       IFAULT  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 0 otherwise
C       RMAX    = output, approximate bound on the accuracy of the diagonal
C                 elements of c.  e.g. if rmax = 1.e-04 then the diagonal
C                 elements of c will be accurate to about 4 dec. digits.
C
C       Latest revision - 18 April 1981

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MNATT*(MNATT+1)/2),
     &          V(MNATT*(MNATT+1)/2)
      DIMENSION W(MNATT)

        NROW=N
        IFAULT=1
        IF(NROW.LE.0) GO TO 100
        IFAULT=0
C       Cholesky factorization of A, result in C
        CALL CHOLAQ(XVAR,NROW,V,NULLTY,IFAULT,RMAX,W,TOL)
        IF(IFAULT.NE.0) GO TO 100
C       Invert C & form the product (CINV)'*CINV, where CINV is the inverse
C       of C, row by row starting with the last row.
C       IROW = the row number, NDIAG = location of last element in the row.
        NN=NROW*(NROW+1)/2
cc        DO 200 IMY=1,NN
cc200       CMY(IMY)=V(IMY)
          DV=1.0
		DO 200 J=1,NROW
	      JJ=J*(J+1)/2
	      DV=DV*V(JJ)*V(JJ)
200       CONTINUE
        IROW=NROW
        NDIAG=NN
10      IF(V(NDIAG).EQ.0.0) GO TO 60
        L=NDIAG
        DO 20 I=IROW,NROW
          W(I)=V(L)
          L=L+I
20      CONTINUE
        ICOL=NROW
        JCOL=NN
        MDIAG=NN
30      L=JCOL
        X=0.0
        IF(ICOL.EQ.IROW) X=1.0/W(IROW)
        K=NROW
40      IF(K.EQ.IROW) GO TO 50
        X=X-W(K)*V(L)
        K=K-1
        L=L-1
        IF(L.GT.MDIAG) L=L-K+1
        GO TO 40
50      V(L)=X/W(IROW)
        IF(ICOL.EQ.IROW) GO TO 80
        MDIAG=MDIAG-ICOL
        ICOL=ICOL-1
        JCOL=JCOL-1
        GO TO 30
60      L=NDIAG
        DO 70 J=IROW,NROW
          V(L)=0.0
          L=L+J
70      CONTINUE
80      NDIAG=NDIAG-IROW
        IROW=IROW-1
        IF(IROW.NE.0) GO TO 10
100     CONTINUE
        RETURN
        END


        SUBROUTINE CHOLAQ(XVAR,N,V,NULLTY,IFAULT,
     &                   RMAX,R,TOL)
C       Modified from
C       Algorithm AS6, Applied Statistics, Vol.17, 1968, with
C       modifications by A.J.Miller
C
C       ARGUMENTS:-
C       XVAR()     = input, a +ve definite matrix stored in lower-triangular
C                 form.
C       N       = input, the order of A
C       V()     = output, a lower triangular matrix such that U*U' = A.
C                 A & U may occupy the same locations.
C       NULLTY  = output, the rank deficiency of A.
C       IFAULT  = output, error indicator
C                       = 1 if N < 1
C                       = 2 if A is not +ve semi-definite
C                       = 0 otherwise
C       RMAX    = output, an estimate of the relative accuracy of the
C                 diagonal elements of U.
C       R()     = output, array containing bounds on the relative accuracy
C                 of each diagonal element of U.
C
C       Latest revision - 18 April 1981
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max' 
      DIMENSION XVAR(MNATT*(MNATT+1)/2),
     &          V(MNATT*(MNATT+1)/2),
     &          R(MNATT)

C       ETA should be set equal to the smallest +ve value such that
C       1.0 + eta is calculated as being greater than 1.0 in the accuracy
C       being used.
C        DATA ETA/1.e-07/
        ETA=TOL
        IFAULT=1
        IF(N.LE.0) GO TO 100
        IFAULT=2
        NULLTY=0
        RMAX=ETA
        R(1)=ETA
        J=1
        K=0
C       Factorize column by column, ICOL = Column No.
        DO 80 ICOL=1,N
          L=0
C       IROW = row number within column ICOL
          DO 40 IROW=1,ICOL
            K=K+1
            W=XVAR(K)
            IF(IROW.EQ.ICOL) RSQ=(W*ETA)**2
            M=J
            DO 10 I=1,IROW
              L=L+1
              IF(I.EQ.IROW) GO TO 20
              W=W-V(L)*V(M)
              IF(IROW.EQ.ICOL) RSQ=RSQ+(V(L)**2*R(I))**2
              M=M+1
10          CONTINUE
20          IF(IROW.EQ.ICOL) GO TO 50
            IF(V(L).EQ.0.0) GO TO 30
            V(K)=W/V(L)
            GO TO 40
30          V(K)=0.0
            IF(ABS(W).GT.ABS(RMAX*XVAR(K))) GO TO 100
40        CONTINUE
C         End of row, estimate relative accuracy of diagonal element.
50        RSQ=SQRT(RSQ)
          IF(ABS(W).LE.5.*RSQ) GO TO 60
          IF(W.LT.0.0) GO TO 100
          V(K)=SQRT(W)
          R(I)=RSQ/W
          IF(R(I).GT.RMAX) RMAX=R(I)
          GO TO 70
60        V(K)=0.0
          NULLTY=NULLTY+1
70        J=J+ICOL
80      CONTINUE
        IFAULT=0.0
100     RETURN
        END
C
C
C   This group of subroutines deal with generating random numbers
C   
C   The subroutines in this file  basically act as a interface between
C   how this program calls the random number generator ( R=RANDNUM(SEED)
C  

C   D.Peel Nov 1995

       SUBROUTINE DETERRANDOM(IER)
C      This subroutine is called at the very beginning of the program
C      to determine if the random number generator is working 
C
C      The result of the test is stored in the common random 
C      variable RANDTYPE
C      This is now defunct as only one random number generator is tried
C      but is left in the code in case a number of generators are available
C      in future versions

	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	   COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	   DIMENSION XISEEDS(1000)
           IER=0
	   DO 100 K=1,1000
             XISEEDS(K)=0
100    CONTINUE
           IX=10023
           IY=324
           IZ=54367
           XISEEDS(1)=.435543
           DO 106 I=2,1000
             XISEEDS(I)=RANDOM(IX,IY,IZ)
             IF (XISEEDS(I).EQ.0) GO TO 107
             DO 105 J=1,I-1
               IF (XISEEDS(I).EQ.XISEEDS(J)) GO TO 107
105      CONTINUE
106    CONTINUE
       RANDTYPE=1
           RETURN

107    WRITE(*,*)'PROBLEM: The random number generator'
       WRITE(*,*)'  does not seem to produce random numbers.'
       WRITE(*,*)'  Applied statistics random number generator failed'
       WRITE(*,*)'  to produce suitable random numbers as well. Modify'
       WRITE(*,*)'  the file nmmrand.f to make this program'
       WRITE(*,*)'  compatible with your inbuilt random number'
       WRITE(*,*)'  generator or alternatively use another random'
       WRITE(*,*)'  number generator.'
       WRITE(*,*)      
       WRITE(*,*)' At present MIXCLUS will still function but you will'
       WRITE(*,*)' be unable to use features that incorporate random'
       WRITE(*,*)' numbers'
       RANDTYPE=0
       IER=40
	   RETURN
	   END
	   
	   FUNCTION RANDNUM()
C          This is the function called by the program NMM. If you
C          wish to use your own portable random number generator 
C          then it should be used in place of this function.

	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	   COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	   IF (RANDTYPE.EQ.1) THEN
             RANDNUM=RANDOM(IX,IY,IZ)
	   ELSE
      WRITE(*,*)'ERROR: As previously described due to random number'
      WRITE(*,*)'       generator problems features utilising'
      WRITE(*,*)'       random numbers are unavailable'
	   STOP
	   ENDIF
           RETURN
	   END

      DOUBLE PRECISION FUNCTION RANDOM(IX,IY,IZ)
C
C     Algorithm AS 183 Appl. Statist. (1982) vol.31, no.2
C
C     Returns a pseudo-random number rectangularly distributed
C     between 0 and 1.   The cycle length is 6.95E+12 (See page 123
C     of Applied Statistics (1984) vol.33), not as claimed in the
C     original article.
C
C     IX, IY and IZ should be set to integer values between 1 and
C     30000 before the first entry.
C
C     Integer arithmetic up to 30323 is required.
C
      INTEGER IX, IY, IZ
c      COMMON /RANDC/ IX, IY, IZ


C
      IX = 171 * MOD(IX, 177) - 2 * (IX / 177)
      IY = 172 * MOD(IY, 176) - 35 * (IY / 176)
      IZ = 170 * MOD(IZ, 178) - 63 * (IZ / 178)
C
      IF (IX .LT. 0) IX = IX + 30269
      IF (IY .LT. 0) IY = IY + 30307
      IF (IZ .LT. 0) IZ = IZ + 30323
C
C     If integer arithmetic up to 5212632 is available, the preceding
C     6 statements may be replaced by:
C
C     IX = MOD(171 * IX, 30269)
C     IY = MOD(172 * IY, 30307)
C     IZ = MOD(170 * IZ, 30323)
C
      RANDOM = MOD(FLOAT(IX) / 30269. + FLOAT(IY) / 30307. +
     +                        FLOAT(IZ) / 30323., 1.0)
      RETURN
      END
C
C
C  This group of subroutines implements the K-means Clustering algorithm.
C  Implemented by David Peel May 1994                      

      SUBROUTINE KMEANS(NIND,NATT,NCOV,NG,X,XMU2,XVAR2,T2,IER)
C      Main subroutine
      
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER T
       EXTERNAL RANDNUM
       INTEGER FLAGS(40),FYLENO
       COMMON /STORE2/ FLAGS,FYLENO
       DOUBLE PRECISION RANDNUM
       INCLUDE 'EMMIX-f1.max'
       REAL X(MNIND,MNATT)
       DIMENSION XKOLD(MAXNG,MNATT),
     &           XSTAN(MNATT),
     &           XMU2(MAXNG,MNATT),
     &           XM(MNATT),XVR(MNATT),XUS(MNATT),
     &           XVAR2(MAXNG,MNATT*(MNATT+1)/2),
     &           T2(MAXNG)
      INTEGER N(MAXNG)
      IER=0
      DO 1 K=1,NG
       N(K)=0.0
       DO 1 I=1,NATT
	XMU2(K,I)=0.0
1     CONTINUE
      CALL KSTAND(NIND,NATT,X,XM,XVR) 
      CALL KSEED(NIND,NATT,NG,X,XKOLD,XM,XVR,IER)

      DO 19 KK=1,NIND
        DO 17 J=1,NATT
          XSTAN(J)=(DBLE(X(KK,J))-XM(J))/SQRT(XVR(J))
17      CONTINUE
       CALL WINNER(NATT,NG,XSTAN,XKOLD,GRP,IER)
       CALL INIT(NIND,NATT,NG,XSTAN,XMU2,GRP,N,IER)
       IF (IER.NE.0) RETURN 
19    CONTINUE
      DO 30 T=1,MAXKM

      DO 430 KK=1,NG
	N(KK)=0
        DO 420 LL=1,NATT 
            XKOLD(KK,LL)=XMU2(KK,LL)
	    XMU2(KK,LL)=0.0
420   	  CONTINUE
430     CONTINUE
	DO 20 KK=1,NIND
          DO 220 J=1,NATT
            XSTAN(J)=(DBLE(X(KK,J))-XM(J))/SQRT(XVR(J))
220       CONTINUE
          CALL WINNER(NATT,NG,XSTAN,XKOLD,GRP,IER)
          CALL INIT(NIND,NATT,NG,XSTAN,XMU2,GRP,N,IER)
	  IF (IER.NE.0) RETURN
20      CONTINUE
          ET=RULE(NG,NATT,XKOLD,XMU2)
          IF (ET.LE.TOLKM) GO TO 99
30    CONTINUE
      
      WRITE (FYLENO,*) 'REACHED MAXIMUM NUMBER OF ',MAXKM,' ITERATIONS'
      IER=-41
99    CONTINUE 
       DO 300 K=1,NG
       DO 300 I=1,NATT
	 XKOLD(K,I)=XMU2(K,I)
         XMU2(K,I)=XMU2(K,I)*SQRT(XVR(I))+XM(I)
300    CONTINUE
      CALL CALC_PARA(NIND,NATT,NCOV,NG,X,XMU2,XKOLD,XVAR2,T2,XM,XVR)
      RETURN
      END
     
      SUBROUTINE CALC_PARA(NIND,NATT,NCOV,NG,X,
     &               XMU2,XKOLD,XVAR2,T2,XM,XVR)
C     This Subroutine calculates estimates of the group covariance
C      and mixing proportions.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      INTEGER NIND,      ! Number of sample points
     &        NATT,      ! Number of attributes/variables/dimensions
     &        NG         ! Number of components or groups
      DIMENSION T2(MAXNG),XMU2(MAXNG,MNATT)
      DIMENSION XKOLD(MAXNG,MNATT)
      DIMENSION XUS(MNATT),XM(MNATT),XVR(MNATT) ,
     &          XVAR2(MAXNG,MNATT*(MNATT+1)/2)
      REAL X(MNIND,MNATT)

      DO 300 K=1,NG
       T2(K)=0.0
       DO 300 I=1,NATT
        DO 300 J=1,I
         XVAR2(K,IC(I,J))=0.0
300   CONTINUE

      DO 410 KK=1,NIND
       DO 400 II=1,NATT
        XUS(II)=(DBLE(X(KK,II))-XM(II))/SQRT(XVR(II))
400    CONTINUE
       CALL WINNER(NATT,NG,XUS,XKOLD,GRP,IER)
       T2(GRP)=T2(GRP)+1.0
       DO 401 II=1,NATT
        XUS(II)=DBLE(X(KK,II))
401    CONTINUE
       DO 330 J=1,NATT
	 DO 330 I=1,J
	  XVAR2(GRP,IC(I,J))=XVAR2(GRP,IC(I,J))+(XUS(I)
     &            	 -XMU2(GRP,I))*(XUS(J)-XMU2(GRP,J))
330    CONTINUE

410   CONTINUE

      DO 500 K=1,NG
       DO 499 I=1,NATT
        DO 499 J=1,I
	 XVAR2(K,IC(I,J))=XVAR2(K,IC(I,J))/T2(K)
499    CONTINUE	 
       T2(K)=T2(K)/FLOAT(NIND)
500   CONTINUE
      RETURN
      END



      
      


      SUBROUTINE KSTAND(NIND,NATT,X,XM,XVR) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX-f1.max'
      REAL X(MNIND,MNATT)
      DIMENSION  XVR(MNATT),XM(MNATT)
      DO 200 J=1,NATT
       XM(J)=0
       DO 200 I=1,NIND
        XM(J)=XM(J)+DBLE(X(I,J))/FLOAT(NIND) 
200   CONTINUE       
      DO 210 J=1,NATT
       XVR(J)=0
       DO 210 I=1,NIND
        XVR(J)=XVR(J)+(DBLE(X(I,J))-XM(J))*
     &            (DBLE(X(I,J))-XM(J))/FLOAT(NIND-1)
210   CONTINUE
      RETURN
      END

      SUBROUTINE KSEED(NIND,NATT,NG,X,XK,XM,XVR,IER)               
c     This Subroutine chooses the initial K seeds (Means of clusters)
c     for the algorithm. At present they are chosen from data set at 
c     random.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER CHOICE
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM
       INCLUDE 'EMMIX-f1.max'
      REAL  X(MNIND,MNATT)
      DIMENSION XK(MAXNG,MNATT)
      DIMENSION  XVR(MNATT),XM(MNATT)
      DO 210 I=1,NG
        R=RANDNUM()
        R=R*NIND
c       Convert CHOICE to integer
	CHOICE=INT(R)+1
	DO 200 J=1,NATT
          XK(I,J)=(DBLE(X(CHOICE,J))-XM(J))/SQRT(XVR(J))
200	CONTINUE
210   CONTINUE
      RETURN
      END

      SUBROUTINE WINNER(NATT,NG,XSTAN,XK,GRP,IER)
c     This subroutine determines the allocation of the KKth point 
c     ie which mean is closest to the given data point (Euclidean).

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX-f1.max'
      DIMENSION XSTAN(MNATT),XK(MAXNG,MNATT)
      DO 310 I=1,NG
        DIST=0
  	DO 300 J=1,NATT
          DIST=DIST+(XSTAN(J)-XK(I,J))**2
300 	CONTINUE
        IF (I.EQ.1) DISTB=DIST 
        IF (DIST.LE.DISTB) THEN
	  GRP=I
          DISTB=DIST
        ENDIF
310   CONTINUE
      RETURN
      END 
          
      SUBROUTINE INIT(NIND,NATT,NG,XSTAN,XK,GRP,N,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX-f1.max'
      DIMENSION XK(MAXNG,MNATT),
     &          XSTAN(MNATT)
      INTEGER N(MAXNG)
c         Update rules
          DO 440 LL=1,NATT
          XK(GRP,LL)=(XK(GRP,LL)*FLOAT(N(GRP))
     &     	      +XSTAN(LL))/(FLOAT(N(GRP)+1))
440       CONTINUE
          N(GRP)=N(GRP)+1
      RETURN
      END

      SUBROUTINE UPDATE(NIND,NATT,NG,XSTAN,XK,GRPOLD,GRP,N,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX-f1.max'
       DIMENSION XSTAN(MNATT),N(MAXNG)
       DIMENSION XK(MAXNG,MNATT)
          IF (N(GRPOLD).GT.1) THEN
c         Update rules
          DO 440 LL=1,NATT
          XK(GRP,LL)=(XK(GRP,LL)
     &  	     *FLOAT(N(GRP))+XSTAN(LL))/FLOAT(N(GRP)+1)
          XK(GRPOLD,LL)=(XK(GRPOLD,LL)
     &	             *FLOAT(N(GRPOLD))-XSTAN(LL))/FLOAT(N(GRPOLD)-1)
440       CONTINUE
          N(GRPOLD)=N(GRPOLD)-1
          N(GRP)=N(GRP)+1
          ELSE
           WRITE (*,*) 'Not enough points in K-means group',GRP
	   IER=1
           RETURN
          ENDIF
      RETURN
      END


      FUNCTION RULE(NG,NATT,XKOLD,XK)
c     This function returns the value used to determine if the algorithm
c     has converged it is a measure of the change in the nodes from iteration 
c     to iteration. 
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX-f1.max'
      INTEGER R
      DIMENSION XK(MAXNG,MNATT),XKOLD(MAXNG,MNATT)
      RULE=0.0
      DO 510 KK=1,NATT
        DO 500 R=1,NG
           RULE=RULE+abs(XK(R,KK)-XKOLD(R,KK))
500     CONTINUE
510   CONTINUE
      RETURN
      END

      SUBROUTINE FACTE(NIND,NATT,NG,X,W,GAMM,XMU,NATTQ)
C     This Subroutine determine the factor weights `e'
C     that can be plotted 
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX-f1.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      REAL X(MNIND,MNATT)
      DIMENSION XMU(MAXNG,MNATT),
     &          W(MNIND,MAXNG),
     &          GAMM(MAXNG,MNATT,MNATTQ),
     &          E(MNIND,NATTQ)
      DO 300 JJ=1,NIND
	DO 100 J=1,NATTQ
	 E(JJ,J)=0
100     CONTINUE
	DO 110 K=1,NG
	  DO 110 I=1,NATTQ
	    DO 110 J=1,NATT
	      E(JJ,I)=E(JJ,I)+W(JJ,K)*GAMM(K,J,I)*
     &              (X(JJ,J)-XMU(K,J))
110       CONTINUE
	WRITE(26,*) (E(JJ,J),J=1,NATTQ)
300   CONTINUE
      RETURN
      END
 

      SUBROUTINE FSAMSP(NIND,NATT,NG,NCOV,XVAR,V2,D,B,
     &                XMU,GAMM,NATTQ,DV2,FACT,T,IER)
C     This Subroutine implements the M-step of the EM algorithm 
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX-f1.max'
      INTEGER FLAGS(40),FYLENO,FACT(3)
      COMMON /STORE2/ FLAGS,FYLENO
      REAL X(MNIND,MNATT)
      DIMENSION XMU(MAXNG,MNATT),
     &        XVAR(MAXNG,MNATT*(MNATT+1)/2),T(MAXNG),
     &        B(MAXNG,MNATT,MNATTQ),D(MAXNG,MNATT),
     &        V2(MAXNG,MNATT*(MNATT+1)/2),DV2(MAXNG),
     &        GAMM(MAXNG,MNATT,MNATTQ),ALPA(MNATT*(MNATT+1)/2)
      DIMENSION TEMP1(MNATT,MNATTQ),
     &          TEMP2(MNATTQ*(MNATTQ+1)/2),
     &          TEMP3(MNATTQ*(MNATTQ+1)/2),
     &          TEMP4(MNATT,MNATTQ),
     &          TD(MNATT)


        DO 2000 K=1,NG

C      Calculate Gamma 
C--------------------------------------------------------------
C	CALL VMULT(V2,0,B,0,NATT,NATT,NATTQ,NG,GAMM)       !inv(SigmaN)*B
	DO 10 L=1,NATTQ
	  DO 10 I=1,NATT
            GAMM(K,I,L)=0.0D0
 	    DO 10 J=1,NATT
	      GAMM(K,I,L)=GAMM(K,I,L)+V2(K,IC(I,J))*B(K,J,L)
10     CONTINUE
C--------------------------------------------------------------
C      Calculate Alpha
C--------------------------------------------------------------
C	CALL MULT(B,1,V2,0,NATTQ,NATT,NATT,NG,TEMP1)       ! (B'*SigmaN
	DO 20 L=1,NATT
	  DO 20 I=1,NATTQ
            TEMP1(L,I)=0.0D0
 	    DO 20 J=1,NATT
	      TEMP1(L,I)=TEMP1(L,I)+B(K,J,I)*V2(K,IC(J,L))
20     CONTINUE
C--------------------------------------------------------------
C	CALL MULT(ALPA,0,B,0,NATTQ,NATT,NATTQ,NG,TEMP)    !  *B)
	DO 30 L=1,NATTQ
	  DO 30 I=1,NATTQ
            ALPA(IC(I,L))=0.0D0
 	    DO 30 J=1,NATT
	      ALPA(IC(I,L))=ALPA(IC(I,L))+TEMP1(J,I)*B(K,J,L)
30     CONTINUE

C--------------------------------------------------------------
C	CALL EYEmA(TEMP,NATTQ,NG,ALPA)                    ! I- ^
	DO 40 I=1,NATTQ
	 DO 40 J=1,I
	   IF (I.EQ.J) THEN
	    ALPA(IC(I,J))=1.0D0-ALPA(IC(I,J))
	   ELSE
	    ALPA(IC(I,J))=-ALPA(IC(I,J))
	   ENDIF
40      CONTINUE
C--------------------------------------------------------------
C      Calculate BT (and hence B)
C--------------------------------------------------------------
C	CALL MULT(GAMM,1,XVAR,0,NATTQ,NATT,NATT,NG,B)  !(Gamm'*Sigma
	DO 50 L=1,NATT
	  DO 50 I=1,NATTQ
            B(K,L,I)=0.0D0
 	    DO 50 J=1,NATT
	      B(K,L,I)=B(K,L,I)+GAMM(K,J,I)*XVAR(K,IC(J,L))
50     CONTINUE
C--------------------------------------------------------------
C	CALL MULT(B,0,GAMM,0,NATTQ,NATT,NATTQ,NG,TEMP)    ! *Gamm)
	DO 55 L=1,NATTQ
	  DO 55 I=1,NATTQ
            TEMP2(IC(I,L))=0.0D0
 	    DO 55 J=1,NATT
	      TEMP2(IC(I,L))=TEMP2(IC(I,L))+B(K,J,I)*GAMM(K,J,L)
55     CONTINUE
C--------------------------------------------------------------
C	CALL ADD(TEMP,0,ALPA,0,NATTQ,NATTQ,NG,B)             ! +alpha
	DO 60 I=1,NATTQ
	  DO 60 J=1,I
            TEMP2(IC(I,J))=TEMP2(IC(I,J))+ALPA(IC(I,J))
60      CONTINUE
C--------------------------------------------------------------
C       Invert result
	  IER=0
	  CALL GDETQ(NATTQ,TEMP2,TEMP3,DV,IER)
	  IF (IER.GT.0) THEN
	    WRITE (*,*) 'ERROR inverting matrix to calculate B'
	    RETURN
	  ENDIF
C--------------------------------------------------------------
C	CALL MULT(TEMP,0,GAMM,1,NATTQ,NATTQ,NATT,NG,B)   ! ^ *Gamm'
        DO 70 L=1,NATT
         DO 70 I=1,NATTQ
          TEMP1(L,I)=0.0D0
          DO 70 J=1,NATTQ
           TEMP1(L,I)=TEMP1(L,I)+TEMP3(IC(I,J))*GAMM(K,L,J)
70      CONTINUE
C--------------------------------------------------------------
C	CALL MULT(B,0,XVAR,0,NATTQ,NATT,NATT,NG,TEMP)     ! *Sigma
        DO 80 L=1,NATTQ
         DO 80 I=1,NATT
          B(K,I,L)=0.0D0
 	  DO 80 J=1,NATT
           B(K,I,L)=B(K,I,L)+XVAR(K,IC(J,I))*TEMP1(J,L)
80      CONTINUE

C--------------------------------------------------------------
C	CALL TRANS(TEMP,NATTQ,NATT,NG,B)              ! Transpose result
C    Done in mult above
C--------------------------------------------------------------

C      Calculate D 
C--------------------------------------------------------------
C	CALL MULT(XVAR,0,GAMM,0,NATT,NATT,NATTQ,NG,TEMP)  ! (Sigma*Gamm
        DO 90 L=1,NATTQ
         DO 90 I=1,NATT
          TEMP1(I,L)=0.0D0
	  DO 90 J=1,NATT
           TEMP1(I,L)=TEMP1(I,L)+XVAR(K,IC(I,J))*GAMM(K,J,L)
90      CONTINUE
C--------------------------------------------------------------
C	CALL MULT(TEMP,0,B,1,NATT,NATTQ,NATT,NG,D)        !  *B)
        DO 95 L=1,NATT
         DO 95 I=1,NATT
          V2(K,IC(I,L))=0.0D0
	  DO 95 J=1,NATTQ
           V2(K,IC(I,L))=V2(K,IC(I,L))+TEMP1(I,J)*B(K,L,J)
95      CONTINUE
C--------------------------------------------------------------
C	CALL MINUS(XVAR,0,D,0,NATT,NATT,NG,TEMP)          ! Sigma - ^
C	CALL DIAG(TEMP,0,NATT,NG,D)                       ! diag{ ^ }
	DO 100 I=1,NATT
         D(K,I)=XVAR(K,IC(I,I))-V2(K,IC(I,I))
100     CONTINUE
C--------------------------------------------------------------

2000    CONTINUE

C      Common D
	IF (FACT(2).EQ.1) THEN
	 DO 3000 JJ=1,NATT
	 TD(JJ)=0.0
	 DO 3000 II=1,NG
	  TD(JJ)=TD(JJ)+T(II)*D(II,JJ)
3000     CONTINUE
	DO 3001 II=1,NG
	DO 3001 JJ=1,NATT
	  D(II,JJ)=TD(JJ)
3001    CONTINUE
	ENDIF

C      PPDA Part
        IF (FACT(2).EQ.2) THEN
          DO 1999 II=1,NG
          DSUM=0.0
          DO 1998 JJ=1,NATT
           DSUM=DSUM+D(II,JJ)/FLOAT(NATT)
1998      CONTINUE
          DO 1999 JJ=1,NATT
           D(II,JJ)=DSUM
1999      CONTINUE
         ENDIF

        DO 2001 K=1,NG
C--------------------------------------------------------------

C        Invert part of Sigma New        
C      Calculate inversion part of Sigma New
C        calculate BT * inv(D) 
C--------------------------------------------------------------
C        calculate BT * inv(D)
	  DO 260 I=1,NATTQ
           DO 260 J=1,NATT
 	    TEMP1(J,I)=B(K,J,I)*1.0/D(K,J)
260       CONTINUE
C--------------------------------------------------------------
C	CALL MULT(BD,0,B,0,NATTQ,NATT,NATTQ,NG,V2)    !   ^ * B
C       CALL MULT2V(TEMP2,1,B,0,NATTQ,NATT,NATTQ,NG,V2)    !   ^ * B
	DO 265 L=1,NATTQ
         DO 265 I=1,L
	   TEMP2(IC(L,I))=0.0D0
           DO 265 J=1,NATT
            TEMP2(IC(L,I))=TEMP2(IC(L,I))+TEMP1(J,I)*B(K,J,L)
265     CONTINUE
C--------------------------------------------------------------
C	CALL EYEaA(V2,NATTQ,NG,TEMP)                  ! I + ^
        DO 266 I=1,NATTQ
	   TEMP2(IC(I,I))=1.0D0+TEMP2(IC(I,I))
266     CONTINUE
C--------------------------------------------------------------
	  IER=0
	  CALL GDETQ(NATTQ,TEMP2,TEMP3,DV,IER)
           
	  IF (IER.GT.0) THEN
	   WRITE (*,*) 'ERROR inverting Sigma New matrix in Mstep'
	   RETURN
	  ENDIF
C--------------------------------------------------------------
C	  CALL MULT(BD,1,V2,0,NATT,NATTQ,NATT,NG,TEMP)
	  DO 120 I=1,NATT
           DO 120 L=1,NATTQ
	    TEMP4(I,L)=0.0D0
	    DO 120 J=1,NATTQ
CC	     TEMP4(I,L)=TEMP4(I,L)+TEMP1(J,I)*V2(K,IC(J,L))
	     TEMP4(I,L)=TEMP4(I,L)+TEMP1(I,J)*TEMP3(IC(J,L))
120     CONTINUE
C--------------------------------------------------------------
C	  CALL MULT(TEMP,0,BD,0,NATT,NATTQ,NATT,NG,V2)
	  DO 130 I=1,NATT
           DO 130 L=1,I
	    V2(K,IC(I,L))=0.0D0
	    DO 130 J=1,NATTQ
	     V2(K,IC(I,L))=V2(K,IC(I,L))+TEMP4(I,J)*TEMP1(L,J)
130     CONTINUE
C--------------------------------------------------------------
	  DO 270 I=1,NATT
	   DO 270 J=1,I
	    IF (I.EQ.J) THEN
	     V2(K,IC(I,J))=1.0/D(K,J)-V2(K,IC(I,J))
	    ELSE
	     V2(K,IC(I,J))=-1*V2(K,IC(I,J))
	    ENDIF
270       CONTINUE
C--------------------------------------------------------------

C       calculate det of Sigma New
C--------------------------------------------------------------
C	  CALL MULT(B,1,V2,0,NATTQ,NATT,NATT,NG,BD)    !B'*inv(B'B+D)
	  DO 140 I=1,NATTQ
           DO 140 L=1,NATT
	    TEMP1(L,I)=0.0D0
	    DO 140 J=1,NATT
	     TEMP1(L,I)=TEMP1(L,I)+B(K,J,I)*V2(K,IC(J,L))
140     CONTINUE
C--------------------------------------------------------------
C	  CALL MULT(BD,0,B,0,NATTQ,NATT,NATTQ,NG,TEMP)  !B'*inv(B'B+D)*B
	  DO 150 I=1,NATTQ
           DO 150 L=1,I
	    TEMP2(IC(I,L))=0.0D0
	    DO 150 J=1,NATT
	     TEMP2(IC(I,L))=TEMP2(IC(I,L))+TEMP1(J,I)*B(K,J,L)
150     CONTINUE
C--------------------------------------------------------------
C	  CALL EYEmA(TEMP,NATTQ,NG,BD)               !I-B'*inv(B'B+D)*B
	  DO 170 I=1,NATTQ
	   DO 160 J=1,I
	    TEMP2(IC(I,J))=-TEMP2(IC(I,J))
160        CONTINUE
	   TEMP2(IC(I,I))=1.0+TEMP2(IC(I,I))
170       CONTINUE
C--------------------------------------------------------------
	  IER=0
           
	  CALL GDETQ(NATTQ,TEMP2,TEMP3,DV,IER)
	  DV2(K)=DV
	  IF (IER.GT.0) THEN
	   WRITE (*,*) 'ERROR Calculating det of Sigma New'
	   RETURN
	  ENDIF
C--------------------------------------------------------------
C         calculate det of D
	   DV2(K)=1/DV2(K)
	   DO 1000 I=1,NATT        
	    DV2(K)=DV2(K)*D(K,I)
1000       CONTINUE

2001    CONTINUE
      RETURN
      END

      SUBROUTINE FSASET2(NIND,NATT,NG,NCOV,X,DV2,XVAR,V2,
     &                  D,B,XMU,NATTQ,MODE,FACT,SETUPON,T,IER)
C     This Subroutine implements the M-step of the EM algorithm 
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX-f1.max'
      INTEGER FLAGS(40),FYLENO,FACT(3)
      COMMON /STORE2/ FLAGS,FYLENO
      REAL X(MNIND,MNATT)
      DIMENSION XMU(MAXNG,MNATT),
     &          T(MAXNG),
     &          B(MAXNG,MNATT,MNATTQ),D(MAXNG,MNATT),
     &          TD(MNATT)

      DIMENSION TEMP1(MNATTQ),
     &          TEMP2(MNATT,MNATTQ),
     &          TEMP3(MNATTQ*(MNATTQ+1)/2),
     &          TEMP4(MNATT,MNATTQ),
     &          TEMP5(MNATTQ*(MNATTQ+1)/2)

      DIMENSION V2(MAXNG,MNATT*(MNATT+1)/2),EIG(MNATT),
     &          XVAR(MAXNG,MNATT*(MNATT+1)/2),
     &          DV2(MAXNG),
     &          EIGVAL(MNATT),EIGVEC(MNATT,MNATT)
      IF (SETUPON.EQ.1) THEN

C     Initialse method by B=0 and D= diag(component covariances)
       IF (FACT(3).EQ.1) THEN 
	write (22,*) 'Using Random'
	CALL GDET(NCOV,NATT,NG,XVAR,V2,DV2,IER)
        DPROD=1
        DO 1200 K=1,NG
 	 DO 1110 I=1,NATT         
	  D(K,I)=V2(K,IC(I,I))
         DPROD=DPROD*D(K,I)/FLOAT(NATT)
1110     CONTINUE
         DO 1200 I=1,NATT
          DO 1200 J=1,NATTQ
           XR=RANDNUM()
           CALL NORM(B(K,I,J),XR)
           IER=0
           B(K,I,J)=B(K,I,J)*SQRT(DPROD)**(1/FLOAT(NATT))/FLOAT(NATTQ) 
1200     CONTINUE
C      Common D
	IF (FACT(2).EQ.1) THEN
	 DO 5000 JJ=1,NATT
	 TD(JJ)=0.0
	 DO 5000 II=1,NG
	  TD(JJ)=TD(JJ)+T(II)*D(II,JJ)
5000     CONTINUE
	DO 5001 II=1,NG
	DO 5001 JJ=1,NATT
	  D(II,JJ)=TD(JJ)
5001    CONTINUE
	ENDIF

       ELSE
	write (22,*) 'Using Eigenvalues'
        
C--------------------------------------------------------------
	DO 110 K=1,NG
 	 DO 110 I=1,NATT         
	  D(K,I)=XVAR(K,IC(I,I))
110      CONTINUE
C      Common D
	IF (FACT(2).EQ.1) THEN
	 DO 4000 JJ=1,NATT
	 TD(JJ)=0.0
	 DO 4000 II=1,NG
	  TD(JJ)=TD(JJ)+T(II)*D(II,JJ)
4000     CONTINUE
	DO 4001 II=1,NG
	DO 4001 JJ=1,NATT
	  D(II,JJ)=TD(JJ)
4001    CONTINUE
	ENDIF

	DO 3000 K=1,NG
C--------------------------------------------------------------
	 DO 100 I=1,NATT       
	  DO 100 J=1,I     
	   V2(K,IC(I,J))=XVAR(K,IC(I,J))*SQRT((1.0/D(K,I))*(1.0/D(K,J)))
100      CONTINUE
C--------------------------------------------------------------
	 CALL TDIAG(NATT,K,V2,EIGVAL,EIG,EIGVEC,IER)
	 CALL LRVT(NATT,EIGVAL,EIG,EIGVEC,IER)
	 IF (NATT.LE.30) THEN
          WRITE(22,*)'Eigenvalues of sample covariance grp ',K 
          WRITE (22,*) (EIGVAL(ii),ii=1,NATT)
	 ELSE
          WRITE(22,*)'First 30 Eigenvalues of sample covariance grp ',K 
          WRITE (22,*) (EIGVAL(ii),ii=1,30)
	 ENDIF
         Esigma=0.0
         DO 102 I=NATTQ+1,NATT
         Esigma=Esigma+EIGVAL(I)/FLOAT(NATT-NATTQ)
102      CONTINUE 	
C--------------------------------------------------------------
	DO 111 I=1,NATTQ
	   TEMP1(I)=(EIGVAL(I)-Esigma)**(0.5)      ! (EIGVAL-I )^1/2
111     CONTINUE
C--------------------------------------------------------------
C	MULT(EIGVEC,0,TEMP1,0,NATT,NATTQ,NATTQ,TEMP)   
	DO 208 I=1,NATT
	DO 208 J=1,NATTQ
	 B(K,I,J)=EIGVEC(I,J)*TEMP1(J)
208     CONTINUE
C--------------------------------------------------------------
	DO 210 I=1,NATT  
	DO 210 J=1,NATTQ  
	   B(K,I,J)=B(K,I,J)*(D(K,I))**(0.5)
210     CONTINUE
C--------------------------------------------------------------
3000   CONTINUE
       ENDIF

       ENDIF

        DO 2000 K=1,NG
C        Invert part of Sigma New        
C      Calculate inversion part of Sigma New
C--------------------------------------------------------------
C        calculate BT * inv(D) 
	  DO 260 I=1,NATTQ
	   DO 260 J=1,NATT
 	    TEMP2(J,I)=B(K,J,I)*1.0/D(K,J) 
260       CONTINUE
C--------------------------------------------------------------
C	CALL MULT2V(TEMP2,1,B,0,NATTQ,NATT,NATTQ,NG,V2)    !   ^ * B
 	DO 265 L=1,NATTQ
         DO 265 I=1,L
          TEMP5(IC(L,I))=0.0D0
          DO 265 J=1,NATT
	   TEMP5(IC(L,I))=TEMP5(IC(L,I))+TEMP2(J,I)*B(K,J,L)
265     CONTINUE
C--------------------------------------------------------------
C	CALL EYEaV(V2,NATTQ,NG,TEMP5)                  ! I + ^
	DO 266 I=1,NATTQ
           TEMP5(IC(I,I))=1.0D0+TEMP5(IC(I,I))
266     CONTINUE
C--------------------------------------------------------------
	  IER=0
	  CALL GDETQ(NATTQ,TEMP5,TEMP3,DV,IER)
	  IF (IER.GT.0) THEN
	   WRITE (*,*) 'ERROR inverting Sigma New matrix in Mstep'
	   RETURN
	  ENDIF
C--------------------------------------------------------------
C	  CALL MULTV(TEMP2,0,TEMP3,0,NATT,NATTQ,NATTQ,TEMP4)
          DO 267 L=1,NATTQ
           DO 267 I=1,NATT
            TEMP4(I,L)=0.0D0
            DO 267 J=1,NATTQ
             TEMP4(I,L)=TEMP4(I,L)+TEMP2(I,J)*TEMP3(IC(J,L))
267     CONTINUE
C--------------------------------------------------------------
C	  CALL VMULT(TEMP,0,BD,0,NATT,NATTQ,NATT,NG,V2)
C	  CALL VMULT(TEMP4,0,TEMP2,1,NATT,NATTQ,NATT,NG,V2)****
	DO 268 L=1,NATT
          DO 268 I=1,NATT
            V2(K,IC(I,L))=0.0D0
            DO 268 J=1,NATTQ
	     V2(K,IC(I,L))=V2(K,IC(I,L))+TEMP4(I,J)*TEMP2(L,J)
268     CONTINUE
C--------------------------------------------------------------
	  DO 270 I=1,NATT
	   DO 270 J=1,I
	    IF (I.EQ.J) THEN
	     V2(K,IC(I,J))=1.0/D(K,J)-V2(K,IC(I,J))
	    ELSE
	     V2(K,IC(I,J))=-1*V2(K,IC(I,J))
	    ENDIF
270       CONTINUE
C--------------------------------------------------------------
C       calculate det of Sigma New
C	  CALL MULTV(B,1,V2,0,NATTQ,NATT,NATT,NG,BD)    !B'*(B'B+D)
C	  CALL MULTV(V2,1,B,0,NATT,NATT,NATTQ,NG,TEMP2)    !B'*(B'B+D)
	   DO 271 L=1,NATTQ
            DO 271 I=1,NATT
             TEMP2(I,L)=0.0D0
             DO 271 J=1,NATT
              TEMP2(I,L)=TEMP2(I,L)+V2(K,IC(I,J))*B(K,J,L)
271       CONTINUE
C--------------------------------------------------------------
C	  CALL MULT(BD,0,B,0,NATTQ,NATT,NATT,NG,TEMP3)  !B'*(B'B+D)*B
C	  CALL VMULT(TEMP2,1,B,0,NATTQ,NATT,NATT,NG,TEMP3)  !B'*(B'B+D)*B
	   DO 272 L=1,NATTQ
            DO 272 I=1,NATTQ
             TEMP3(IC(I,L))=0.0D0
             DO 272 J=1,NATT
              TEMP3(IC(I,L))=TEMP3(IC(I,L))+TEMP2(J,I)*B(K,J,L)
272       CONTINUE
C--------------------------------------------------------------
C	  CALL EYEmA(TEMP,NATTQ,NG,BD)               !I-B'*(B'B+D)*B
C	  CALL EYEmV(TEMP3,NATTQ,NG,TEMP5)               !I-B'*(B'B+D)*B
	  DO 273 I=1,NATTQ
	   DO 273 J=1,I
	    IF (I.EQ.J) THEN
	     TEMP3(IC(I,J))=1.0-TEMP3(IC(I,J))
	    ELSE
	     TEMP3(IC(I,J))=-1*TEMP3(IC(I,J))
	    ENDIF
273       CONTINUE
C--------------------------------------------------------------
	  IER=0

	  CALL GDETQ(NATTQ,TEMP3,TEMP5,DV,IER)
	  DV2(K)=DV
	  IF (IER.GT.0) THEN
	   WRITE (*,*) 'ERROR Calculating det of Sigma New'
	   RETURN
	  ENDIF
C--------------------------------------------------------------
C         calculate det of D
	   DV2(K)=1/DV2(K)
	   DO 1000 I=1,NATT        
	    DV2(K)=DV2(K)*D(K,I)
1000       CONTINUE
2000    CONTINUE
      RETURN
      END    


c      SUBROUTINE LOOP(NIND,NATT,NG,X,XMU,V,XVAR,DV,T,NCOV,IER,TXML,IDT,
c     &   WL,W,XUU,USA,TOLS,NATTQ,B,D)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cC     This subroutine uses the EM algorithm from a specified starting
cC     value to find a solution of the likelihood equation.
c         CALL FSASET2(NIND,NATT,NG,NCOV,X,XVAR,DV,DV2,V2,D,B,
c     &               XMU,NATTQ,T,IER)
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cC       MAIN ITERATIVE LOOP
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       IF ((FLAGS(22).EQ.100).OR.(FLAGS(23).EQ.2)) THEN
c       CALL ESTEP(NIND,NATT,NG,X,XMU,V2,T,DEN,WL,W,XUU,USA,DV2,
c     &            XLOGL,IOUNT,XMAH,IER)
cC       RE-Calculate Mu and Pi and XVAR
c      CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,
c     &           XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER)
c       ELSE
c       CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
c     &            XLOGL,IOUNT,XMAH,IER)
c       ENDIF
c
c      IF ((FLAGS(22).EQ.100).OR.(FLAGS(23).EQ.2)) THEN
c        IER=0
c        CALL FSAMSP(NIND,NATT,NG,NCOV,XVAR,V2,D,B,
c     &                XMU,GAMM,NATTQ,DV2,IER)
c        IF (IER.NE.0) RETURN
c       KK=1
c13010   CONTINUE
c13020   CONTINUE
c      ELSE
c      CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,
c     &           XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER)
cccc      ENDIF
      SUBROUTINE TDIAG (N,KK,A,D,E,Z,IFAULT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C 
       INCLUDE 'EMMIX-f1.max'

      DIMENSION A(MAXNG,MNATT*(MNATT+1)/2), D(MNATT), E(MNATT) 
      DIMENSION Z(MNATT,MNATT)
      DATA ZERO/0.0/,ONE/1.0/
      DATA ETA/1.0D-37/,EPS/1.0D-14/
C 
C       Algorithm as 60.1 appl.statist. (1973) vol.22 no.2
C 
C        reduces real symmetric matrix to tridiagonal form
C 
C       tol is a machine dependent constant , tol = eta/eps , where
C       eta is the smallest positive number representable in the
C       computer and eps is the smallest positive number for which
C       1+eps.ne.1.
C 
C         eta=eps*tol
C         eps=0.7105427358e-14
C         tol=0.3131513063e-293
C         precis=1.0e-14
C 
C         nb
C           real constants must be le 15 decimal digits
C           the range of a real constant is from 1.0e-293 to 1.0e+322
C 
      TOL=ETA/EPS
      IFAULT=1
      IF (N.LE.1) RETURN
      IFAULT=0
      DO 10 I=1,N
        DO 10 J=1,I
 10     Z(I,J)=A(KK,IC(I,J))
      I=N
      DO 110 I1=2,N
        L=I-2
        F=Z(I,I-1)
        G=ZERO
        IF (L.LT.1) GO TO 30
        DO 20 K=1,L
 20       G=G+Z(I,K)**2
 30     H=G+F*F
C 
C       if g is too small for orthogonality to be guaranteed, the
C       transformation is skipped
C 
        IF (G.GT.TOL) GO TO 40
        E(I)=F
        D(I)=ZERO
        GO TO 100
 40     L=L+1
	G=SQRT(H)
        IF (F.GE.ZERO) G=-G
        E(I)=G
        H=H-F*G
        Z(I,I-1)=F-G
        F=ZERO
        DO 80 J=1,L
          Z(J,I)=Z(I,J)/H
          G=ZERO
C 
C       form element of a * u
C 
          DO 50 K=1,J
 50         G=G+Z(J,K)*Z(I,K)
          IF (J.GE.L) GO TO 70
          J1=J+1
          DO 60 K=J1,L
 60         G=G+Z(K,J)*Z(I,K)
C 
C       form element of p
C 
 70       E(J)=G/H
          F=F+G*Z(J,I)
 80     CONTINUE
C 
C       form k
C 
        HH=F/(H+H)
C 
C       form reduced a
C 
        DO 90 J=1,L
          F=Z(I,J)
          G=E(J)-HH*F
          E(J)=G
          DO 90 K=1,J
          Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
 90     CONTINUE
        D(I)=H
 100    I=I-1
 110  CONTINUE
      D(1)=ZERO
      E(1)=ZERO
C 
C       accumulation of transformation matrices
C 
      DO 160 I=1,N
        L=I-1
        IF (D(I).EQ.ZERO.OR.L.EQ.0) GO TO 140
        DO 130 J=1,L
          G=ZERO
          DO 120 K=1,L
 120        G=G+Z(I,K)*Z(K,J)
          DO 130 K=1,L
          Z(K,J)=Z(K,J)-G*Z(K,I)
 130    CONTINUE
 140    D(I)=Z(I,I)
        Z(I,I)=ONE
        IF (L.EQ.0) GO TO 160
        DO 150 J=1,L
          Z(I,J)=ZERO
          Z(J,I)=ZERO
 150    CONTINUE
 160  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE LRVT (N,D,E,Z,IFAULT)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       INCLUDE 'EMMIX-f1.max'
C 
C       algorithm as 60.2 appl.statist. (1973) vol.22, no.2
C 
C        finds latent roots and vectors of tridiagonal matrix
C 
      DIMENSION D(MNATT), E(MNATT), Z(MNATT,MNATT)
      DATA MITS/30/,ZERO/0.0/,ONE/1.0/,TWO/2.0/
C 
      PRECIS=1.0D-14
      IFAULT=2
      IF (N.LE.1) RETURN
      IFAULT=1
      N1=N-1
      DO 10 I=2,N
 10     E(I-1)=E(I)
      E(N)=ZERO
      B=ZERO
      F=ZERO
      DO 100 L=1,N
        JJ=0
	H=PRECIS*(ABS(D(L))+ABS(E(L)))
        IF (B.LT.H) B=H
C 
C       look for small sub-diagonal element
C 
        DO 20 M1=L,N
          M=M1
	  IF (ABS(E(M)).LE.B) GO TO 30
 20     CONTINUE
 30     IF (M.EQ.L) GO TO 90
 40     IF (JJ.EQ.MITS) RETURN
        JJ=JJ+1
C 
C       form shift
C 
        P=(D(L+1)-D(L))/(TWO*E(L))
	R=SQRT(P*P+ONE)
        PR=P+R
        IF (P.LT.ZERO) PR=P-R
        H=D(L)-E(L)/PR
        DO 50 I=L,N
 50       D(I)=D(I)-H
        F=F+H
C 
C       ql transformation
C 
        P=D(M)
        C=ONE
        S=ZERO
        M1=M-1
        I=M
        DO 80 I1=L,M1
          J=I
          I=I-1
          G=C*E(I)
          H=C*P
	  IF (ABS(P).GE.ABS(E(I))) GO TO 60
          C=P/E(I)
	  R=SQRT(C*C+ONE)
          E(J)=S*E(I)*R
          S=ONE/R
          C=C/R
          GO TO 70
 60       C=E(I)/P
          R=SQRT(C*C+ONE)
          E(J)=S*P*R
          S=C/R
          C=ONE/R
 70       P=C*D(I)-S*G
          D(J)=H+S*(C*G+S*D(I))
C 
C       form vector
C 
          DO 80 K=1,N
          H=Z(K,J)
          Z(K,J)=S*Z(K,I)+C*H
          Z(K,I)=C*Z(K,I)-S*H
 80     CONTINUE
        E(L)=S*P
        D(L)=C*P
	IF (ABS(E(L)).GT.B) GO TO 40
 90     D(L)=D(L)+F
 100  CONTINUE
C 
C       order latent roots and vectors
C 
      DO 130 I=1,N1
        K=I
        P=D(I)
        I1=I+1
        DO 110 J=I1,N
          IF (D(J).LE.P) GO TO 110
          K=J
          P=D(J)
 110    CONTINUE
        IF (K.EQ.I) GO TO 130
        D(K)=D(I)
        D(I)=P
        DO 120 J=1,N
          P=Z(J,I)
          Z(J,I)=Z(J,K)
          Z(J,K)=P
 120    CONTINUE
 130  CONTINUE
      IFAULT=0
      RETURN
      END
      SUBROUTINE NORM(Q,R)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X,Y,G3X,V1,V2,SS,CL,Z
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM,R
      R=RANDNUM()
      IF(R.LE..8638D0) THEN
10       R=RANDNUM()
         X=R
         R=RANDNUM()
         X=R+X
         R=RANDNUM()
         X=X+R
         Z=2.0D0*(X-1.50D0)
      ELSEIF(R.LE.0.9745D0) THEN
         R=RANDNUM()
         X=R
         R=RANDNUM()
         X=X+R
         Z=1.5D0*(X-1.D0)
      ELSEIF(R.LE.0.9973002039D0) THEN
20       R=RANDNUM()
         X=6.0D0*R-3.0D0
         R=X
         Y=0.358D0*RANDNUM()
         R=X
         IF(ABS(X).LT.1.0D0) THEN
          G3X=17.49731196D0*DEXP(-0.5D0*X*X)-4.73570326D0*(3D0-X*X)
     &        -2.15787544D0*(1.5D0-ABS(X))
         ELSEIF(ABS(X).LT.1.5D0) THEN
          G3X=17.49731196D0*DEXP(-0.5D0*X*X)-2.36785163D0*(3.0D0-ABS(X))
     &        **2.0D0-2.15787544D0*(1.5D0-ABS(X))
         ELSEIF(ABS(X).LT.3.0D0) THEN
          G3X=17.49731196D0*DEXP(-0.5D0*X*X)-2.36785163D0*(3.0D0-ABS(X))
     &       **2.0D0
         ELSE
          G3X=0.0D0
         ENDIF
         IF(Y.LT.G3X)GO TO 30
         GO TO 20
30       CONTINUE
         Z=X
        ELSE
40      V1=RANDNUM()
        V2=RANDNUM()
        V1=2.0D0*V1-1.0D0
        V2=2.0D0*V2-1.0D0
        SS=V1*V1+V2*V2
        IF(SS.GT.1.0D0)GO TO 40
        CL=SQRT(9.0-2*LOG(SS)/SS)
        X=CL*V1
        Y=CL*V2
        IF(ABS(X).GT.3.0D0) THEN
          Z=X
          GOTO 45
        ELSEIF(ABS(Y).GT.3.0D0) THEN
          Z=Y
          Q=Z
          RETURN
        ENDIF
        GO TO 40
       ENDIF
45     Q=Z
      RETURN
      END

