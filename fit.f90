!************ FORTRAN DEMO PROGRAM ***********************************************
PROGRAM DLSFIT_Demo
   REAL,ALLOCATABLE::x(:),y(:),sig(:),a(:),da(:),covar(:,:)
   INTEGER,ALLOCATABLE::spos(:),ia(:)
   EXTERNAL::FPOLY_lin,FPOLY_non_lin,FGAUSJ  !User supplied routines
   EXTERNAL::FCL
   LOGICAL lf,sigs_available,rescale_sigs
   INTEGER::sbno,sno

   REAL,PARAMETER::k=2.      ! DLS exponent k
   REAL,PARAMETER::res=1.D-6  ! Parameter used in indefinite case; irrelevant if expk=2;
      ! otherwise should be of the same order of magnitude as resolution of measurement;
      ! in WLS has to be scaled by minimal sig(i).

   n0=11
   ALLOCATE(x(n0));ALLOCATE(y(n0));ALLOCATE(sig(n0));ALLOCATE(spos(n0))
   !=========================================================================
   ! Data used in Example 1 (see Fig.2)
   !=========================================================================
   x=(/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11./)
   y=(/1.,2.1,3.1,3.9,4.8,7.5,6.85,6.5,9.15,11.,11.1/)
   ! sig contains uncertainties of y-data points; if unknown, set all sig(i)=1
   sig=(/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)  ! 'sig' must be supplied
   ! For demo purpose, these 'sig' are overestimated.
   !=========================================================================
   ! End of experimental data
   !=========================================================================

   !=========================================================================
   ! Example of input parameters
   !=========================================================================
   lf=.false.     ! .TRUE. for linear fit, .FALSE. for non-linear fit
   ma=2           ! Fit to straight line;
   ALLOCATE(a(ma));ALLOCATE(da(ma));ALLOCATE(ia(ma));ALLOCATE(covar(ma,ma))
   a=(/1.D0,0.D0/)   ! Trial parameters; irrelevant for linear fit
   ia=1             ! No frozen parameters
   rp=0.9            ! remove parameter
   sigs_available=.TRUE.
   !=========================================================================
   ! End of input parameters
   !=========================================================================

   CALL Echo_Of_Program_Name  ! Writes the program name on the standard output unit (screen)


   ! The next program line calls DLSFIT;
   ! The user must supply the name of model function subroutine which must agree with the value of 'lf'
   ! ('FPOLY_lin' if 'lf=.true.', and 'FPOLY_non_lin' if 'lf=.false' in this example)
   CALL DLSFIT(x,y,sig,n0,a,ma,rp,chisq,dls,db,ns,spos,sbno,covar,ia,FPOLY_non_lin,lf,k,res)

   ! NOTE: errors in the best-fit parameters are not handled by DLSFIT !!!
   ! However, it returns matrix covar which can be used for calculation of errors.

   ! Error estimation in the best-fit parameters:
   !=============================================
   !
   ! Case 1: experimental errors are NOT available; in that case enable the following
   !         program line to calculate errors and rescale 'chisq', see Sec. 4(i)
   !
   ! sig0=db/zbrent(FCL,.1,3.,1D-11);FORALL(i=1:ma) da(i)=SQRT(covar(i,i))*sig0;chisq=chisq/(sig0*sig0)
   !
   !
   ! Case 2: experimental errors are AVAILABLE;
   !
   !        in that case calculate firstly sig0 by enabling the next line
   ! sig0=db/zbrent(FCL,.1,3.,1D-11)
   !
   !        In ideal case sig0=1
   !
   !        If individual experimental errors are correctly estimated, then sig0 should be
   !        ---------------------------------------------------------
   !        close to 1 (say in interval 0.5 to 1); in that case calculate the errors by the
   !        following line:
   !
   ! FORALL(i=1:ma)  da(i)=SQRT(covar(i,i))
   !
   !
   !        If sig0 is notably different from 1 (say less than 0.5 or greater than 1.5) then
   !        -----------------------------------
   !        either reconsider experimental errors sig(i) - which is preferable
   !
   !        or
   !
   !        apply the next line:
   !
   ! FORALL(i=1:ma)  da(i)=SQRT(covar(i,i))*sig0
   !
   !        which effectively calculates errors 'da' with rescaled sig(i)*sig0 - see Sec. 4(ii)
   !        and rescale 'chisq' accordingly
   !
   ! chisq=chisq/(sig0*sig0)


   ! The following code efficiently summarize above options.
   ! =======================================================
   !
   sig0=db/zbrent(FCL,.1,3.,1D-11)
   IF(sigs_available)THEN
      IF ((sig0<0.5).OR.(sig0>1.5))THEN
         PRINT"(A,F8.4,A)"," sig0 = ",REAL(sig0,4)," is notably different from 1"
         PRINT"("" Rescale experimental errors (T/F) = ""\)"
         READ(*,*) rescale_sigs
         PRINT*
         IF(.NOT.rescale_sigs) sig0=1.
         chisq=chisq/(sig0*sig0)
      END IF
   ELSE
      chisq=chisq/(sig0*sig0)
   END IF
   FORALL(i=1:ma) da(i)=SQRT(covar(i,i))*sig0


   ! The following code illustrates how to obtain fitting parameters on arbitrary
   ! ============================================================================
   ! subset from the ordered collection C; 'sno' is ordinal number of subset
   ! =======================================================================
   !
   ! CALL fit(x,y,sig,spos(sno),a,ia,ma,covar,chisq,FPOLY_non_lin,lf)


   CALL DLS_Fit_Report  ! Writes the DLS fit report on the standard output unit (screen)
   STOP
CONTAINS
   SUBROUTINE Echo_Of_Program_Name
      ! Writes the program name on the standard output unit (screen)
      OPEN(UNIT=6,CARRIAGECONTROL='FORTRAN') ! Standard output unit (screen)
      WRITE(6,*)""
      WRITE(6,*)"Program DLSFIT_AA_Demo for DLS robust data fitting"
      WRITE(6,*)"--------------------------------------------------"
      WRITE(6,*)""
   END SUBROUTINE Echo_Of_Program_Name
   SUBROUTINE DLS_Fit_Report
      ! Writes the DLS fit report on the standard output unit (screen)
      CHARACTER(10)::int2str
      WRITE(6,"(""+Best fit parameters: A=                      "")")
       DO i=1,ma
         WRITE(6,*)i,":",a(i),achar(241),da(i)
       END DO
      WRITE(6,*)""
      WRITE(6,*)"DLS is max on subset ",TRIM(int2str(sbno))," (from ",TRIM(int2str(ns)),")"
      WRITE(6,*)"Number of close points = ",TRIM(int2str(spos(sbno)))
      WRITE(6,*)"DLSmax =",dls
      WRITE(6,*)"Width of best subset Db =",db
      WRITE(6,*)"CHISQ on best subset =",chisq
      WRITE(6,*)""
      WRITE(6,"("" Program DLSFIT_AA_Demo is finished. Press any key...""\)")
      READ(*,*)
      CLOSE(UNIT=6)
   END SUBROUTINE DLS_Fit_Report
END PROGRAM DLSFIT_Demo
!*********************************************************************************
SUBROUTINE DLSFIT(x,y,sig,n0,a,m,rp,&
                & chisq,dls,db,ns,spos,sbno,&
                & covar,ia,funcs,&
                & lf,k,res)

   ! Input parameters:
   ! =================
   ! x(1:n0), y(1:n0) - set of n0 experimental data points;
   ! sig(1:n0)        - individual standard deviations sig(1:n0) of their y_i's
   !                    (if not available, set all sig(i)=1 prior to DLSFIT call);
   ! a(1:m)           - array of m initial values for model function parameters
   !                    (irrelevant for linear fit, but must be supplied);
   ! rp               - removal parameter
   ! ia               - for frozen a(j) set ia(j)=0, otherwise ia(j) is non-zero;
   ! funcs            - the name of a user supplied subroutine. Following Numerical Recipes
   !                    methods, funcs returns m basis functions evaluated at x in the array 
   !                    afunc(1:m) for general-linear WLS. Otherwise, it returns the calculated 
   !                    value f(x;a) of model function in variable y, together with the values 
   !                    of all model function partial derivatives, returned in array dyda(1:m).
   ! lf               - .TRUE. for general-linear WLS data fitting, otherwise .FALSE.
   ! k                - DLS exponent;
   ! res              - the original resolution of measurement;
   !
   ! Output parameters:
   ! ==================    
   ! a(1:m)           - the best-fit parameters associated with the best subset;
   ! chisq            - chi-square obtained on the best subset;
   ! dls              - the density of least-squares obtained on the best subset;
   ! db               - the width of the best subset;
   ! ns               - total number of subsets in ordered collection C
   ! spos             - integer array containing the number of data in each subset from C;
   ! sbno             - the ordinal number of the best subset in C
   ! covar(1:m,1:m)   - covariance matrix
   !
   ! Note: x(1:n0), y(1:n0) are reordered on output so that the first spos(i) data points belong 
   ! to the subset s_i in C.            

!***********************************************************************************
   ! Remove the next line if Compaq compiler is not used
!   USE numerical_libraries  !Compaq compiler directive
!***********************************************************************************
   REAL,PARAMETER:: Pi=3.1415926535897932384626433832795D0, q=1./3.
   INTEGER ia(m),spos(n0),sbno
   REAL x(n0),y(n0),sig(n0),k
   REAL a(m),covar(m,m),dyda(m),afunc(m),da(m)
   EXTERNAL funcs
   LOGICAL lf
   IF(k>2.9.OR.k<2.)RETURN
   dtiny=(MAXVAL(y)-MINVAL(y))/10.**PRECISION(db);nt=n0;ns=1;dls=0.
   CALL fit(x,y,sig,nt,a,ia,m,covar,chisq,funcs,lf)
   DO  WHILE (nt>m+3)
      CALL NewMeritVals;IF(dbt==0) EXIT
      CALL ArrangeLayer
   END DO
   CALL fit(x,y,sig,spos(sbno),a,ia,m,covar,chisq,funcs,lf)
CONTAINS
   SUBROUTINE ArrangeLayer
      ! Rearranges arrays x,y,sig so that their first spos(i)
      ! elements belong to the i-th subset in ordered collection
      dr=rp*dbt;ns=ns+1
      DO
         i=1;nto=nt
         DO  WHILE (i<=nt) !Remove sublayer
            yc=fun(x(i))
            IF(ABS(y(i)-yc)/sig(i)>=dr)THEN
               yc=fun(x(nt))
               IF (ABS(y(nt)-yc)/sig(i)<dr)THEN
                  swap=x(i);x(i)=x(nt);x(nt)=swap
                  swap=y(i);y(i)=y(nt);y(nt)=swap
                  swap=sig(i);sig(i)=sig(nt);sig(nt)=swap
               END IF
               nt=nt-1
            ELSE
               i=i+1
            END IF
         END DO
         IF((nt==nto).OR.(nt<=m+3))EXIT
         CALL fit(x,y,sig,nt,a,ia,m,covar,chisq,funcs,lf)
      END DO
   END SUBROUTINE ArrangeLayer
   SUBROUTINE NewMeritVals
      ! Calculates new DLS value
      spos(ns)=nt;dbt=0
      DO i=1,nt
         d=ABS(y(i)-fun(x(i)))/sig(i); IF(dbt<d) dbt=d
      END DO
      IF(dbt>dtiny)THEN
         dlst=chisq/(dbt**k)
      ELSE
         os2=SUM((/(1./(sig(i)*sig(i)),i=1,nt)/))/nt
         dlst=os2*(1.+q*(nt-1))*res**(2-k);dbt=0
      END IF
      IF(dlst>dls)THEN
         dls=dlst;sbno=ns;db=dbt
      END IF
   END SUBROUTINE NewMeritVals
   FUNCTION fun(x)   ! For given x, calculates value of model function
      IF(lf)THEN
         CALL funcs(x,afunc,m);fun=DOT_PRODUCT(a,afunc)
      ELSE
         CALL funcs(x,a,fun,dyda,m)
      END IF
   END FUNCTION fun
END SUBROUTINE DLSFIT
!*********************************************************************************
 SUBROUTINE fit(x,y,sig,n0,a,ia,m,covar,chisq,funcs,lf)
   ! If lf=.TRUE. performs linear fit by calling lfit
   ! Otherwise, drives mrqmin which performs non-linear fit
   DIMENSION x(n0),y(n0),sig(n0)
   DIMENSION a(m),ia(m),covar(m,m),al(m,m)
   REAL:: lamda,lamdaO
   EXTERNAL funcs
   LOGICAL lf
   PARAMETER (MaxIt=100,RelErr=1.D-5)
   IF(lf)THEN
      CALL lfit(x,y,sig,n0,a,ia,m,covar,m,chisq,funcs);RETURN
   END IF
   lamda=-1;chisqO=0;lamdaO=0;it=0
   DO WHILE (it.LE.MaxIt)
      CALL mrqmin(x,y,sig,n0,a,ia,m,covar,al,m,chisq,funcs,lamda)
      IF(chisq<=0)EXIT
      IF((ABS(chisqO-chisq)/chisq<RelErr).AND.(lamda<lamdaO))EXIT
      chisqO=chisq;lamdaO=lamda;it=it+1
   END DO
   lamda=0.
   CALL mrqmin(x,y,sig,n0,a,ia,m,covar,al,m,chisq,funcs,lamda)
 END SUBROUTINE fit
!*********************************************************************************
FUNCTION FCL(x)
   ! FCL(x) is the function which appears on the left-hand side of Eq. 6.
   REAL,PARAMETER::k=2.      ! DLS exponent k
   PARAMETER (Pi=3.1415926535897932384626433832795D0)
   FCL=x*(x*x+k)*EXP(-x*x/2)-k*SQRT(Pi/2)*DERF(x/SQRT(2.))
END FUNCTION FCL
!*********************************************************************************
SUBROUTINE FPOLY_lin(x,p,np)
   ! A user supplied routine for linear fit to polynomial model function;
   ! for input value x, it returns the vector of basis functions - here,
   ! powers of x: p(1)=1, p(2)=x, p(3)=x*x,..,p(np)=x**(np-1).
   REAL p(np)
   p(1)=1.
   DO j=2,np
      p(j)=p(j-1)*x
   END DO
END
!*********************************************************************************
SUBROUTINE FPOLY_non_lin(X,A,Y,DYDA,MA)
   ! A user supplied routine for non-linear fit to polynomial model function;
   ! for input value X, it returns Y=(polynomial value) and
   ! DYDA(1)=1, DYDA(2)=X,..,DYDA(MA)=X**(MA-1).
   DIMENSION A(MA),DYDA(MA)
   CALL FPOLY_lin(X,DYDA,MA)
   DYDA(1)=1.;Y=A(1)
   DO I=2,MA
      DYDA(I)=DYDA(I-1)*X
      Y=Y+A(I)*DYDA(I)
   END DO
END
!*********************************************************************************
SUBROUTINE FGAUSJ(X,A,Y,DYDA,MA)
   ! A user supplied routine for non-linear fit to Gaussian model function;
   ! for input value X, it returns the value Y and the vector of partial
   ! derivatives DYDA of Gaussian model function, specified by MA=4 parameters:
   ! height A(1), center A(2), width A(3), and offset A(4).
   DIMENSION A(MA),DYDA(MA)
   REAL,PARAMETER::TWOSQRLN2=1.6651092223153955127063292897904
   B=TWOSQRLN2*(X-A(2))/A(3)
   DYDA(1)=EXP(-B*B)
   Y=A(1)*DYDA(1)
   DYDA(2)=2.*Y*B*TWOSQRLN2/A(3)
   DYDA(3)=2.*Y*B*B/A(3)
   Y=Y+A(4)
   DYDA(4)=1.
END SUBROUTINE FGAUSJ
!*********************************************************************************
CHARACTER*(*) FUNCTION int2str(i)
   ! Converts integer to string
   CHARACTER(10)::dummy
   WRITE(dummy,'(I10)')i
   int2str=adjustl(dummy)
END FUNCTION int2str
!*********************************************************************************
SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,funcs,alamda)
   ! From Numerical Recipes Sec. 15-5
   ! Levenberg-Marquardt method, attempting to reduce the value chisq of a fit between a set
   ! of data points x(1:ndata), y(1:ndata) with individual standard deviations sig(1:ndata),
   ! and a nonlinear function dependent on ma coefficients a(1:ma). The input array ia(1:ma)
   ! indicates by nonzero entries those components of a(1:ma) that should be fitted for, and
   ! by zero entries those components that should be held fixed at their input values. The
   ! program returns current best-fit values for the parameters a(1:ma), and chisq. The arrays
   ! covar(1:nca,1:nca), alpha(1:nca,1:nca) with physical dimension nca (>= the number
   ! of fitted parameters) are used as working space during most iterations. Supply a
   ! subroutine funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit, and its
   ! derivatives dyda with respect to the fitting parameters a(1:ma) at x. On the first call
   ! provide an initial guess for the parameters a(1:ma), and set alamda<0 for initialization
   ! (which then sets alamda=.001). If a step succeeds chisq becomes smaller and alamda
   ! decreases by a factor of 10. If a step fails alamda grows by a factor of 10. You must
   ! call this routine repeatedly until convergence is achieved. Then, make one final call
   ! with alamda=0, so that covar(1:ma,1:ma) returns the covariance matrix, and alpha the
   ! curvature matrix. (Parameters held fixed will return zero covariances.)
   INTEGER ma,nca,ndata,ia(ma),MMAX
   REAL alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),sig(ndata),x(ndata),y(ndata)
   EXTERNAL funcs
   PARAMETER (MMAX=20)
   INTEGER j,k,l,mfit
   REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
   SAVE ochisq,atry,beta,da,mfit
   if(alamda.lt.0.)then
      mfit=0
      do j=1,ma
         if (ia(j).ne.0) mfit=mfit+1
      end do
      alamda=0.001
      call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
      ochisq=chisq
      do j=1,ma
         atry(j)=a(j)
      end do
   end if
   do j=1,mfit
      do k=1,mfit
      covar(j,k)=alpha(j,k)
      end do
      covar(j,j)=alpha(j,j)*(1.+alamda)
      da(j)=beta(j)
   end do
   call gaussj(covar,mfit,nca,da,1,1)
   if(alamda.eq.0.)then
      call covsrt(covar,nca,ma,ia,mfit)
      call covsrt(alpha,nca,ma,ia,mfit)
      return
   end if
   j=0
   do l=1,ma
      if(ia(l).ne.0) then
         j=j+1
         atry(l)=a(l)+da(j)
      endif
   end do
   call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
   if(chisq.lt.ochisq)then
      alamda=0.1*alamda
      ochisq=chisq
      do j=1,mfit
         do k=1,mfit
            alpha(j,k)=covar(j,k)
         end do
         beta(j)=da(j)
      end do
      do l=1,ma
         a(l)=atry(l)
      end do
   else
      alamda=10.*alamda
      chisq=ochisq
   endif
   return
END
!**********************************************************************************
SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,funcs)
   ! From Numerical Recipes Sec. 15-5. Computes the matrix alpha and vector beta.
   ! mrqcof calls the user-supplied routine funcs(x,a,y,dyda), which for input values x
   ! and parameters a(1:ma) returns the value of model function y=y(x;a) and the vector
   ! of partial derivatives dyda.
   INTEGER ma,nalp,ndata,ia(ma),MMAX
   REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),y(ndata)
   EXTERNAL funcs
   PARAMETER (MMAX=20)
   INTEGER mfit,i,j,k,l,m
   REAL dy,sig2i,wt,ymod,dyda(MMAX)
   mfit=0
   do j=1,ma
      if (ia(j).ne.0) mfit=mfit+1
   end do
   do j=1,mfit
      do k=1,j
         alpha(j,k)=0.
      end do
      beta(j)=0.
   end do
   chisq=0.
   do i=1,ndata
      call funcs(x(i),a,ymod,dyda,ma)
      sig2i=1./(sig(i)*sig(i))
      dy=y(i)-ymod
      j=0
      do l=1,ma
         if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do m=1,l
               if(ia(m).ne.0) then
                  k=k+1
                  alpha(j,k)=alpha(j,k)+wt*dyda(m)
               end if
            end do
            beta(j)=beta(j)+dy*wt
         end if
      end do
      chisq=chisq+dy*dy*sig2i
   end do
   do j=2,mfit
      do k=1,j-1
         alpha(k,j)=alpha(j,k)
      end do
   end do
   return
END
!**********************************************************************************
SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq,funcs)
   ! From Numerical Recipes Sec. 15-4.
   ! Given a set of data points x(1:ndat), y(1:ndat) with individual standard deviations
   ! sig(1:ndat), use chi-square minimization to fit for some or all of the coefficients
   ! a(1:ma)" of function that depends linearly on coefficients a(1:ma). The input array
   ! ia(1:ma) indicates those components that should be held fixed at their input values.
   ! The program returns values for a(1:ma), chisq, and the covariance matrix covar(1:ma,1:ma).
   ! (Parameters held fixed will return zero covariances). npc is the physical dimension of
   ! covar(npc,npc) in the calling routine. The user supplies a subroutine funcs(x,afunc,ma)
   ! that returns the ma basis functions evaluated at x in the array afunc.
   INTEGER ma,ia(ma),npc,ndat,MMAX
   REAL chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
   EXTERNAL funcs
   PARAMETER (MMAX=50)
   INTEGER i,j,k,l,m,mfit
   REAL sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)
   mfit=0
   do j=1,ma
      if(ia(j).ne.0) mfit=mfit+1
   end do
   if(mfit.eq.0) pause "lfit: no parameters to be fitted"
   do j=1,mfit
      do k=1,mfit
         covar(j,k)=0.
      end do
      beta(j)=0.
   end do
   do i=1,ndat
      call funcs(x(i),afunc,ma)
      ym=y(i)
      if(mfit.lt.ma)then
         do j=1,ma
            if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
         end do
      endif
      sig2i=1./sig(i)**2
      j=0
      do l=1,ma
         if(ia(l).ne.0) then
            j=j+1
            wt=afunc(l)*sig2i
            k=0
            do m=1,l
               if (ia(m).ne.0) then
                  k=k+1
                  covar(j,k)=covar(j,k)+wt*afunc(m)
               endif
            end do
            beta(j)=beta(j)+ym*wt
         endif
      end do
   end do
   do j=2,mfit
      do k=1,j-1
         covar(k,j)=covar(j,k)
      end do
   end do
   call gaussj(covar,mfit,npc,beta,1,1)
   j=0
   do l=1,ma
      if(ia(l).ne.0)then
         j=j+1
         a(l)=beta(j)
      endif
   end do
   chisq=0.
   do i=1,ndat
      call funcs(x(i),afunc,ma)
      sum=0.
         do j=1,ma
            sum=sum+a(j)*afunc(j)
         end do
      chisq=chisq+((y(i)-sum)/sig(i))**2
   end do
   call covsrt(covar,npc,ma,ia,mfit)
   return
END
!*********************************************************************************
SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
   ! From Numerical Recipes Sec. 15-4.
   ! Subroutine covsrt spreads the covariances back into the full ma by ma covariance
   ! matrix, in the proper rows and columns and with zero variances and covariances set
   ! for variables which were held frozen.
   INTEGER ma,mfit,npc,ia(ma)
   REAL covar(npc,npc)
   INTEGER i,j,k
   REAL swap
   do i=mfit+1,ma
      do j=1,i
         covar(i,j)=0.
         covar(j,i)=0.
      end do
   end do
   k=mfit
   do j=ma,1,-1
      if(ia(j).ne.0)then
         do i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
         end do
         do i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
         end do
         k=k-1
      endif
   end do
   return
END
!*********************************************************************************
SUBROUTINE gaussj(a,n,np,b,m,mp)
   ! From Numerical Recipes Sec. 2-1.
   ! Linear equation solution by Gauss-Jordan elimination. a(1:n,1:n) is an input matrix stored
   ! in an array of physical dimensions np by np. b(1:n,1:m) is an input matrix containing the m
   ! right-hand side vectors, stored in an array of physical dimensions np by mp. On output, a(1:n,1:n)
   ! is replaced by its matrix inverse, and b(1:n,1:m) is replaced by the corresponding set of solution
   ! vectors. Parameter: NMAX is the largest anticipated value of n.
   INTEGER m,mp,n,np,NMAX
   REAL a(np,np),b(np,mp)
   PARAMETER (NMAX=50)
   INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
   REAL big,dum,pivinv
   do j=1,n
      ipiv(j)=0
   end do
   do i=1,n
      big=0.
      do j=1,n
         if(ipiv(j).ne.1)then
            do k=1,n
               if(ipiv(k).eq.0) then
                  if(abs(a(j,k)).ge.big)then
                     big=abs(a(j,k))
                     irow=j
                     icol=k
                  endif
               endif
            end do
         endif
      end do
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
         do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
         end do
         do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
         end do
      endif
      indxr(i)=irow
      indxc(i)=icol
      if(a(icol,icol).eq.0.) pause "singular matrix in gaussj"
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n
         a(icol,l)=a(icol,l)*pivinv
      end do
      do l=1,m
         b(icol,l)=b(icol,l)*pivinv
      end do
      do ll=1,n
         if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do l=1,n
               a(ll,l)=a(ll,l)-a(icol,l)*dum
            end do
            do l=1,m
               b(ll,l)=b(ll,l)-b(icol,l)*dum
            end do
         endif
      end do
   end do
   do l=n,1,-1
      if(indxr(l).ne.indxc(l))then
         do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
         end do
      endif
   end do
   return
END
!*********************************************************************************
FUNCTION zbrent(func,x1,x2,tol)
   ! From Numerical Recipes Sec. 9-3.
   ! Using Brent’s method, find the root of a function func known to lie between x1 and x2.
   ! The root, returned as zbrent, will be refined until its accuracy is tol.
   ! Parameters: Maximum allowed number of iterations, and machine floating-point precision.
   INTEGER ITMAX
   REAL zbrent,tol,x1,x2,func,EPS
   EXTERNAL func
   PARAMETER (ITMAX=1000,EPS=3.e-17)
   INTEGER iter
   REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=func(a)
   fb=func(b)
   if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause "root must be bracketed for zbrent"
   c=b
   fc=fb
   do iter=1,ITMAX
      if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
         c=a
         fc=fa
         d=b-a
         e=d
      endif
      if(abs(fc).lt.abs(fb)) then
         a=b
         b=c
         c=a
         fa=fb
         fb=fc
         fc=fa
      endif
      tol1=2.*EPS*abs(b)+0.5*tol
      xm=.5*(c-b)
      if(abs(xm).le.tol1 .or. fb.eq.0.)then
         zbrent=b
         return
      endif
      if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
         s=fb/fa
         if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
         else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
         endif
         if(p.gt.0.) q=-q
         p=abs(p)
         if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
         else
            d=xm
            e=d
         endif
      else
         d=xm
         e=d
      endif
      a=b
      fa=fb
      if(abs(d) .gt. tol1) then
         b=b+d
      else
         b=b+sign(tol1,xm)
      endif
      fb=func(b)
   end do
   pause "zbrent exceeding maximum iterations"
   zbrent=b
   return
END
!*********************************************************************************
