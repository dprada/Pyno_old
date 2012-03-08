MODULE WAT

INTEGER::switch
INTEGER::nw,natw
INTEGER::nparts,nparts2
INTEGER::hb_def

DOUBLE PRECISION::sk_param,roh_param,roo_param,cos_angooh_param


DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE :: Lbox,Lbox2

DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: XARR    ! posiciones (molecula,atomo,coordenada)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DARR    !! (index_water,Hi,num_neights)
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IARR    !! (index_water,Hi,num_neights)
DOUBLE PRECISION,    ALLOCATABLE, DIMENSION(:,:,:,:) :: vect_norm_htoo    !! (index_water,Hi,num_neights)
DOUBLE PRECISION,    ALLOCATABLE, DIMENSION(:,:) :: wat_perp

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: num_h2o, o2h, o2which
INTEGER, ALLOCATABLE, DIMENSION(:) :: num_o2h
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: h2o
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: strength_h2o
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: strength_o2h

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: hbsmol,shell_w
INTEGER, ALLOCATABLE, DIMENSION(:) :: ms_short,ms_short2,mss_ind_wat

DOUBLE PRECISION::pi

PARAMETER (pi=acos(-1.0d0))
PARAMETER (natw=3)
PARAMETER (nparts=8)
PARAMETER (nparts2=nparts*nparts*2*2)


CONTAINS


!##################################################################################
!##################################################################################


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! GENERAL FUNCTIONS FOR EVERY HBOND DEFINITION

!!!
!!! MEMORY
!!!

SUBROUTINE INITIALIZE_COORS_MEMORY()
  ALLOCATE(xarr(nw,3,3),iarr(nw,2,nparts),darr(nw,2,nparts),Lbox(3,3),Lbox2(3,3))
  xarr=0.0d0
  iarr=0
  darr=0.0d0
  Lbox=0.0d0
  Lbox2=0.0d0

END SUBROUTINE INITIALIZE_COORS_MEMORY

SUBROUTINE FREE_COORS_MEMORY()
  DEALLOCATE(xarr,iarr,darr,Lbox,Lbox2)
END SUBROUTINE FREE_COORS_MEMORY

SUBROUTINE INITIALIZE_HBONDS_MEMORY()
  ALLOCATE(num_h2o(nw,2),num_o2h(nw),h2o(nw,2,6),o2h(nw,6),o2which(nw,6))
  ALLOCATE(strength_o2h(nw,6),strength_h2o(nw,2,6))
  num_h2o=0
  num_o2h=0
  o2h=0
  h2o=0
  o2which=0
  strength_h2o=0.0d0
  strength_o2h=0.0d0
  ALLOCATE(vect_norm_htoo(nw,2,nparts,3))
  ALLOCATE(wat_perp(nw,3))
  vect_norm_htoo=0.0d0
  wat_perp=0.0d0
END SUBROUTINE INITIALIZE_HBONDS_MEMORY

SUBROUTINE FREE_HBONDS_MEMORY()
  DEALLOCATE(num_h2o,num_o2h,h2o,o2h,o2which)
  DEALLOCATE(strength_o2h,strength_h2o)
  DEALLOCATE(vect_norm_htoo,wat_perp)
END SUBROUTINE FREE_HBONDS_MEMORY

SUBROUTINE INITIALIZE_MSS_MEMORY()
  ALLOCATE(hbsmol(nw,6),shell_w(nw,17),ms_short(17),ms_short2(17),mss_ind_wat(17))     ! numero maximo de moleculas vecinas en las dos capas
  hbsmol=0
  shell_w=0
  ms_short=0
  ms_short2=0
  mss_ind_wat=0
END SUBROUTINE INITIALIZE_MSS_MEMORY

SUBROUTINE FREE_MSS_MEMORY()
  DEALLOCATE(hbsmol,shell_w,ms_short2,ms_short,mss_ind_wat)
END SUBROUTINE FREE_MSS_MEMORY


!!!
!!! AUXILIARY GENERAL FUNCTIONS
!!!

SUBROUTINE PBC(vector)

  implicit none
  integer::i
  DOUBLE PRECISION,dimension(3),intent(INOUT)::vector
  
  DO i=1,3
     IF (abs(vector(i))>Lbox2(i,i)) THEN
        IF (vector(i)>Lbox2(i,i)) THEN
           vector(i)=vector(i)-Lbox(i,i)
        ELSE
           vector(i)=vector(i)+Lbox(i,i)
        END IF
     END IF
  END DO
  
END SUBROUTINE PBC

SUBROUTINE ALL_NORM_WATER()

  IMPLICIT NONE
  INTEGER::i
  DOUBLE PRECISION,DIMENSION(3)::norm

  DO i=1,NW
     CALL PERPENDICULAR_WATER(i,norm)
     wat_perp(i,:)=norm
  END DO

END SUBROUTINE ALL_NORM_WATER

SUBROUTINE PERPENDICULAR_WATER (a,vect)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::a
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::vect
  DOUBLE PRECISION,DIMENSION(3)::posh1,posh2

  vect=0.0d0

  posh1=xarr(a,1,:)-xarr(a,2,:)
  posh2=xarr(a,1,:)-xarr(a,3,:)

  CALL PRODUCT_VECT(posh1,posh2,vect)
  CALL NORMALIZE_VECT (vect)

END SUBROUTINE PERPENDICULAR_WATER

SUBROUTINE PRODUCT_VECT(a,b,normal)

  DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::a,b
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::normal
  DOUBLE PRECISION::norm
  
  normal(1)=a(2)*b(3)-a(3)*b(2)
  normal(2)=-a(1)*b(3)+a(3)*b(1)
  normal(3)=a(1)*b(2)-a(2)*b(1)

END SUBROUTINE PRODUCT_VECT

SUBROUTINE NORMALIZE_VECT (a)

  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::a
  DOUBLE PRECISION::norm

  norm=sqrt(dot_product(a,a))
  a=a/norm

END SUBROUTINE NORMALIZE_VECT


SUBROUTINE DMAT()

  IMPLICIT NONE
  INTEGER::i,j,jj,g
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter_Hs
  INTEGER,DIMENSION(:),ALLOCATABLE::lista
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dHs
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_aux
  DOUBLE PRECISION,DIMENSION(3)::aux,aux2
  DOUBLE PRECISION::norm,val

  ALLOCATE(filter_Hs(NW),dHs(NW),vect_aux(NW,3),lista(nparts))

  filter_Hs=.true.
  lista=0
  dHs=0.0d0
  vect_aux=0.0d0
  aux=0.0d0
  aux2=0.0d0

  DO i=1,NW

     filter_Hs(i)=.false.

     DO jj=1,2
        aux2=XARR(i,jj+1,:)
        DO j=1,NW
           aux=XARR(j,1,:)-aux2
           CALL PBC (aux)
           dHs(j)=sqrt(dot_product(aux,aux))
           vect_aux(j,:)=aux
        END DO

        DO j=1,nparts
           g=MINLOC(dHs(:),DIM=1,MASK=filter_Hs(:))
           lista(j)=g
           norm=dHs(g)
           DARR(i,jj,j)=norm
           vect_norm_htoo(i,jj,j,:)=vect_aux(g,:)/norm
           filter_Hs(g)=.false.
        END DO
        IARR(i,jj,:)=lista(:)
        
        DO j=1,nparts
           g=lista(j)
           filter_Hs(g)=.true.
        END DO
     END DO
     filter_Hs(i)=.true.

  END DO

  DEALLOCATE(lista)
  DEALLOCATE(filter_Hs,dHs)
  
END SUBROUTINE DMAT

SUBROUTINE DMAT_EF()

  IMPLICIT NONE
  INTEGER::i,j,jj,g,oo,h,hh,gg
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter_Hs
  LOGICAL,DIMENSION(:),ALLOCATABLE::aux_filter
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dHs
  DOUBLE PRECISION,DIMENSION(3)::vect_aux,aux2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_aux2
  INTEGER,DIMENSION(:),ALLOCATABLE::list,list2
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE::iarr2
  DOUBLE PRECISION::norm,val

  ALLOCATE(list(nparts2),list2(nparts),filter_Hs(nparts2),aux_filter(NW),dHs(nparts2),iarr2(NW,2,nparts))
  iarr2=0
  dHs=0.0d0
  list=0
  list2=0
  aux_filter=.false.
  DARR=0.0d0
  filter_Hs=.false.

  DO i=1,NW

     !! List of second neighbors of molecule i:
     list=0
     oo=0
     DO j=1,2
        DO h=1,nparts
           g=IARR(i,j,h)
           aux_filter(g)=.true.

           DO jj=1,2
              DO hh=1,nparts
                 gg=IARR(g,jj,hh)
                 aux_filter(gg)=.true.
              END DO
           END DO

        END DO
     END DO

     aux_filter(i)=.false.

     DO j=1,NW
        IF (aux_filter(j)==.true.) THEN
           oo=oo+1
           list(oo)=j
           aux_filter(j)=.false.
        END IF
     END DO

     !! Checking the second neighbors of molecule i:

     ALLOCATE(vect_aux2(oo,3))
     vect_aux2=0.0d0
     filter_Hs(1:oo)=.true.
     DO j=1,2
        aux2=XARR(i,j+1,:)

        DO hh=1,oo
           h=list(hh)
           vect_aux=XARR(h,1,:)-aux2
           CALL PBC (vect_aux)
           vect_aux2(hh,:)=vect_aux(:)
           dHs(hh)=sqrt(dot_product(vect_aux,vect_aux))
        END DO

        DO jj=1,nparts
           g=MINLOC(dHs(:),DIM=1,MASK=filter_Hs(:))
           list2(jj)=g
           IARR2(i,j,jj)=list(g)
           norm=dHs(g)
           DARR(i,j,jj)=norm
           vect_norm_htoo(i,j,jj,:)=vect_aux2(g,:)/norm
           filter_Hs(g)=.false.
        END DO
        DO jj=1,nparts
           g=list2(jj)
           filter_Hs(g)=.true.
        END DO
     END DO
     filter_Hs(1:oo)=.false.

     DEALLOCATE(vect_aux2)

END DO


DEALLOCATE(list,list2,filter_Hs,aux_filter,dHs)

IARR=IARR2 
DEALLOCATE(iarr2)

END SUBROUTINE DMAT_EF

!##################################################################################
!##################################################################################

!!!
!!! HBONDS OF A BOX OF WATER
!!!

SUBROUTINE hbonds_box ()

  Lbox2=Lbox/2.0d0

  !! Optimization for hbonds

  IF (switch==0) THEN
     CALL DMAT()
     switch=1
  ELSE
     CALL DMAT_EF()
  END IF

  SELECT CASE (hb_def)
  CASE (1)
     CALL ALL_NORM_WATER ()
     CALL HBONDS_SKINNER()
     CALL HBONDS_COMPLETE_OXYGENS()
     CALL SORTING_HBONDS_MAXLOC()
  CASE (2)
     CALL HBONDS_DIST_OH ()
     CALL HBONDS_COMPLETE_OXYGENS()
     CALL SORTING_HBONDS_MINLOC()
  CASE (3)
     CALL HBONDS_DIST_OO_ANG_OOH ()
     CALL HBONDS_COMPLETE_OXYGENS()
     CALL SORTING_HBONDS_MINLOC()
  CASE (4)
     CALL HBONDS_DONOR_ACCEPTOR_NUMBER()
  CASE (5)
     CALL HBONDS_TOPOLOGICAL()
  CASE (6)
     CALL HBONDS_DONOR_ANG_OOH()
  CASE (7)
     CALL HBONDS_NEAREST_NEIGHBOUR()
  CASE DEFAULT
     PRINT*, 'Error: Hbond definition unknown'
  END SELECT

END SUBROUTINE hbonds_box



SUBROUTINE HBONDS_COMPLETE_OXYGENS()

  IMPLICIT NONE
  INTEGER::i,j,jj,g,h

  DO i=1,NW
     DO j=1,2
        DO jj=1,num_h2o(i,j)
           g=h2o(i,j,jj)
           h=num_o2h(g)+1
           num_o2h(g)=h
           o2h(g,h)=i
           o2which(g,h)=j
           strength_o2h(g,h)=strength_h2o(i,j,jj)
        END DO
     END DO
  END DO

END SUBROUTINE HBONDS_COMPLETE_OXYGENS

SUBROUTINE SORTING_HBONDS_MAXLOC()

  IMPLICIT NONE
  INTEGER::i,j,jj,g,h
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::back1,back2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::back3

  ALLOCATE(filtro(6),back1(6),back2(6),back3(6))

  filtro=.false.

  DO i=1,NW
     ! Hydrogens
     DO j=1,2
        g=num_h2o(i,j)
        IF (g>1) THEN
           back1=h2o(i,j,:)
           back3=strength_h2o(i,j,:)
           filtro(1:g)=.true.
           DO jj=1,g
              h=MAXLOC(back3(:),DIM=1,MASK=filtro)
              h2o(i,j,jj)=back1(h)
              strength_h2o(i,j,jj)=back3(h)
              filtro(h)=.false.
           END DO
        END IF
     END DO
     !Oxygens
     g=num_o2h(i)
     IF (g>1) THEN
        back1=o2h(i,:)
        back2=o2which(i,:)
        back3=strength_o2h(i,:)
        filtro(1:g)=.true.
        DO jj=1,g
           h=MAXLOC(back3(:),DIM=1,MASK=filtro)
           o2h(i,jj)=back1(h)
           o2which(i,jj)=back2(h)
           strength_o2h(i,jj)=back3(h)
           filtro(h)=.false.
        END DO
     END IF

  END DO

  DEALLOCATE(filtro,back1,back2,back3)

END SUBROUTINE SORTING_HBONDS_MAXLOC

SUBROUTINE SORTING_HBONDS_MINLOC()

  IMPLICIT NONE
  INTEGER::i,j,jj,g,h
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::back1,back2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::back3

  ALLOCATE(filtro(6),back1(6),back2(6),back3(6))

  filtro=.false.

  DO i=1,NW
     ! Hydrogens
     DO j=1,2
        g=num_h2o(i,j)
        IF (g>1) THEN
           back1=h2o(i,j,:)
           back3=strength_h2o(i,j,:)
           filtro(1:g)=.true.
           DO jj=1,g
              h=MINLOC(back3(:),DIM=1,MASK=filtro)
              h2o(i,j,jj)=back1(h)
              strength_h2o(i,j,jj)=back3(h)
              filtro(h)=.false.
           END DO
        END IF
     END DO
     !Oxygens
     g=num_o2h(i)
     IF (g>1) THEN
        back1=o2h(i,:)
        back2=o2which(i,:)
        back3=strength_o2h(i,:)
        filtro(1:g)=.true.
        DO jj=1,g
           h=MINLOC(back3(:),DIM=1,MASK=filtro)
           o2h(i,jj)=back1(h)
           o2which(i,jj)=back2(h)
           strength_o2h(i,jj)=back3(h)
           filtro(h)=.false.
        END DO
     END IF

  END DO

  DEALLOCATE(filtro,back1,back2,back3)

END SUBROUTINE SORTING_HBONDS_MINLOC



!!!
!!! MICROSTATES OF A BOX OF WATER
!!!

SUBROUTINE MICROSTATES (num_wat,mss)

!  IMPLICIT NONE 

  INTEGER,INTENT(IN)::num_wat
  INTEGER,DIMENSION(num_wat,17),INTENT(OUT)::mss
  INTEGER::j
  
  CALL INITIALIZE_HBONDS_MEMORY()
  CALL INITIALIZE_MSS_MEMORY()

  CALL HBONDS_BOX()

  CALL BUILD_HBONDS_LIMIT_NOSIMETRIC ()
  DO j=1,NW      
     CALL STATE_SHORT_NOSIMETRIC(j)   
     CALL REMOVE_PERMUT_SHORT_NOSIMETRIC(j)
     mss(j,:)=ms_short2(:)
  END DO

  CALL FREE_HBONDS_MEMORY()
  CALL FREE_MSS_MEMORY()

END SUBROUTINE MICROSTATES


SUBROUTINE MICROSTATES_IND_WAT (num_wat,mss)

!  IMPLICIT NONE 

  INTEGER,INTENT(IN)::num_wat
  INTEGER,DIMENSION(num_wat,17),INTENT(OUT)::mss
  INTEGER::j
  
  CALL INITIALIZE_HBONDS_MEMORY()
  CALL INITIALIZE_MSS_MEMORY()

  CALL HBONDS_BOX()

  CALL BUILD_HBONDS_LIMIT_NOSIMETRIC ()
  DO j=1,NW      
     CALL STATE_SHORT_NOSIMETRIC(j)   
     CALL REMOVE_PERMUT_SHORT_NOSIMETRIC(j)
     mss(j,:)=mss_ind_wat(:)-1
  END DO

  CALL FREE_HBONDS_MEMORY()
  CALL FREE_MSS_MEMORY()

END SUBROUTINE MICROSTATES_IND_WAT


!##################################################################################
!##################################################################################

!!!
!!! AUXILIARY FUNCTIONS FOR MICROSTATES
!!!


SUBROUTINE BUILD_HBONDS_LIMIT_NOSIMETRIC ()

  IMPLICIT NONE

  INTEGER::i,ii,j,jj,g,gg,h,hh,b,iii

  !pongo los hidrogenos

  DO i=1,NW
     DO j=1,2
        g=h2o(i,j,1)
        hbsmol(i,j)=g
     END DO
     DO j=1,2
        hbsmol(i,j+2)=o2h(i,j)
        hbsmol(i,j+4)=o2which(i,j)
     END DO
  END DO

END SUBROUTINE BUILD_HBONDS_LIMIT_NOSIMETRIC

SUBROUTINE STATE_SHORT_NOSIMETRIC(a)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::a
  INTEGER::i,j,g,b
  INTEGER,DIMENSION(17)::microstate,aux
  LOGICAL,DIMENSION(17)::filtro

  microstate=0
  microstate(1)=a
  
  microstate(2:5)=(/ hbsmol(a,1:4) /)

  g=6
  DO i=2,3
     b=microstate(i)
     IF (b>0) THEN 
        microstate(g:g+1)=(/ hbsmol(b,1:2) /)
        IF ((hbsmol(b,3)==a).or.(hbsmol(b,4)==a)) THEN
           IF (hbsmol(b,3)==a) THEN
              microstate(g+2)=hbsmol(b,4)
           ELSE
              microstate(g+2)=hbsmol(b,3)
           END IF
           IF ((hbsmol(b,3)==a).and.(hbsmol(b,4)==a)) THEN
              !print*,'problema ssn1'     !!! Comprobar esto
              microstate(g+2)=hbsmol(b,3)
           END IF
        ELSE
           microstate(g+2)=hbsmol(b,3)
        END IF
        g=g+3
     ELSE
        microstate(g:g+2)=0
        g=g+3
     END IF
  END DO
  DO i=4,5
     b=microstate(i)
     IF (b>0) THEN
        microstate(g+1:g+2)=(/ hbsmol(b,3:4) /)
        IF ((hbsmol(b,1)==a).or.(hbsmol(b,2)==a)) THEN
           IF (hbsmol(b,1)==a) THEN
              microstate(g)=hbsmol(b,2)
           ELSE
              microstate(g)=hbsmol(b,1)
           END IF
           IF ((hbsmol(b,1)==a).and.(hbsmol(b,2)==a)) THEN
!              print*,'pproblema ssn2' !!! Comprobar esto
!              print*,a,'|',hbsmol(b,:)
              microstate(g)=hbsmol(b,1)
!              print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
!                 & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)
           END IF
        ELSE
           microstate(g)=hbsmol(b,1) !!!! esto esta mal, hay que elegir entre el 1 y el 2
           !print*,'siiiiiiiiiiiiiiiiii'
        END IF
        g=g+3
     ELSE
        microstate(g:g+2)=0
        g=g+3
     END IF
  END DO

!print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
!                 & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)

  !Quito indices de mols

  filtro=.true.
  aux=microstate
  mss_ind_wat=microstate
  DO i=1,17
     IF (aux(i)<=0) filtro(i)=.false.
  END DO

  shell_w(a,:)=0
  g=0

  DO j=1,17
     IF (filtro(j)==.true.) THEN
        b=aux(j)
        g=g+1
        shell_w(a,g)=b
        DO i=1,17
           IF (filtro(i)==.true.) THEN
              IF (aux(i)==b) THEN
                 microstate(i)=j
                 filtro(i)=.false.
              END IF
           END IF
        END DO
     END IF
  END DO

  ms_short2=microstate


END SUBROUTINE STATE_SHORT_NOSIMETRIC

SUBROUTINE REMOVE_PERMUT_SHORT_NOSIMETRIC(mol)

  IMPLICIT NONE

  INTEGER::i,j,g,h,ii,vigilo
  INTEGER,INTENT(IN)::mol
  INTEGER,DIMENSION(17)::microstate,bb,key,key_aux
  LOGICAL,DIMENSION(17)::filtro
  INTEGER,DIMENSION(4,2)::aux_permut
  LOGICAL::interruptor
  INTEGER::ceros4,ceros5,ceros,x_ellos
  INTEGER::x_primera4,x_primera5,x_primera
  INTEGER::x_core,x_core4,x_core5
  INTEGER::x_segunda4,x_segunda5,x_segunda

  key=mss_ind_wat
  key_aux=mss_ind_wat
  microstate=ms_short2
  aux_permut(1,1)=6
  aux_permut(1,2)=7
  aux_permut(2,1)=9
  aux_permut(2,2)=10
  aux_permut(3,1)=13
  aux_permut(3,2)=14
  aux_permut(4,1)=16
  aux_permut(4,2)=17


  
  DO h=1,2
     i=aux_permut(h,1)
     j=aux_permut(h,2)
     IF (ms_short2(i)>0) THEN
        IF (ms_short2(i)>ms_short2(j)) THEN
           microstate(i)=ms_short2(j)
           microstate(j)=ms_short2(i)
           key_aux(i)=key(j)
           key_aux(j)=key(i)
           filtro=.true.
           DO g=1,17
              IF ((microstate(g)==i).and.(filtro(g)==.true.)) THEN
                 microstate(g)=j
                 filtro(g)=.false.
              END IF
           END DO
           DO g=1,17
              IF ((microstate(g)==j).and.(filtro(g)==.true.)) THEN
                 microstate(g)=i
                 filtro(g)=.false.
              END IF
           END DO
           ms_short2=microstate
           key=key_aux
        END IF
     END IF
  END DO
  
  key=key_aux
  ms_short2=microstate



  ceros4=0
  ceros5=0
  ceros=0
  x_core4=0
  x_core5=0
  x_core=0
  x_primera4=0
  x_primera5=0
  x_primera=0
  x_segunda4=0
  x_segunda5=0
  x_segunda=0
  x_ellos=0

  DO i=12,14
     IF (ms_short2(i)==0) THEN
        ceros4=ceros4+1
     END IF
     IF (ms_short2(i)==1) THEN
        x_core4=x_core4+1
     END IF
     IF ((ms_short2(i)<=5).and.(ms_short2(i)>=2)) THEN
        x_primera4=x_primera4+1
     END IF
     IF ((ms_short2(i)<=11).and.(ms_short2(i)>=6)) THEN
        x_segunda4=x_segunda4+1
     END IF
     IF ((ms_short2(i)/=i).and.(ms_short2(i)>11)) THEN
        x_ellos=x_ellos+1
     END IF
  END DO
  DO i=15,17
     IF (ms_short2(i)==0) THEN
        ceros5=ceros5+1
     ELSE
        IF (ms_short2(i)==1) THEN
           x_core5=x_core5+1
        ELSE
           IF ((ms_short2(i)<=5).and.(ms_short2(i)>=2)) THEN
              x_primera5=x_primera5+1
           ELSE
              IF ((ms_short2(i)<=11).and.(ms_short2(i)>=6)) THEN
                 x_segunda5=x_segunda5+1
              ELSE
                 IF ((ms_short2(i)/=i).and.(ms_short2(i)>11)) THEN
                    x_ellos=x_ellos+1
                 END IF
              END IF
           END IF
        END IF
     END IF
  END DO
  ceros=ceros4+ceros5
  x_core=x_core4+x_core5
  x_primera=x_primera4+x_primera5
  x_segunda=x_segunda4+x_segunda5


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interruptor=.false.
  
  IF (ms_short2(5)==0) THEN
     interruptor=.true.
  END IF
  
  IF (interruptor==.false.) THEN
     IF (x_core/=0) THEN
        IF (x_core==1) THEN
           IF (x_core5==1) THEN
              CALL DOY_VUELTA()
              CALL DOY_VUELTA_KEY (key,key_aux)
              microstate=ms_short2
              interruptor=.true.
           END IF
        ELSE
           !print*,'problema x_core',frame,mol
           !print 117, ms_short2(:)
           !STOP
        END IF
     END IF
  END IF

  IF (interruptor==.false.) THEN
     IF (ceros4/=ceros5) interruptor=.true.
     IF (ceros5<ceros4) THEN
        CALL DOY_VUELTA()
        CALL DOY_VUELTA_KEY (key,key_aux)
     END IF
  END IF
  
  IF (interruptor==.false.) THEN
     IF (x_primera4/=x_primera5) interruptor=.true.
     IF (x_primera5<x_primera4) THEN
        CALL DOY_VUELTA()
        CALL DOY_VUELTA_KEY (key,key_aux)
     END IF
  END IF
  
  IF (interruptor==.false.) THEN
     IF (x_segunda4/=x_segunda5) interruptor=.true.
     IF (x_segunda5<x_segunda4) THEN
        CALL DOY_VUELTA()
        CALL DOY_VUELTA_KEY (key,key_aux)
     END IF
  END IF
  
  !!
  microstate=ms_short2
  key_aux=key
  DO h=3,4
     i=aux_permut(h,1)
     j=aux_permut(h,2)
     IF (ms_short2(i)/=0) THEN
        IF (ms_short2(i)>ms_short2(j)) THEN
           key_aux(i)=key(j)
           key_aux(j)=key(i)
           microstate(i)=ms_short2(j)
           microstate(j)=ms_short2(i)
           filtro=.true.
           DO g=1,17
              IF ((microstate(g)==i).and.(filtro(g)==.true.)) THEN
                 microstate(g)=j
                 filtro(g)=.false.
              END IF
           END DO
           DO g=1,17
              IF ((microstate(g)==j).and.(filtro(g)==.true.)) THEN
                 microstate(g)=i
                 filtro(g)=.false.
              END IF
           END DO
           ms_short2=microstate
           key=key_aux
        END IF
     END IF
  END DO
  ms_short2=microstate
  key=key_aux
  !!


  !sigo
  
  IF (interruptor==.false.) THEN
     IF (ceros>0) THEN
        IF ((ms_short2(12)==0).and.(ms_short2(15)/=0)) THEN
           interruptor=.true.
        ELSE
           IF ((ms_short2(12)/=0).and.(ms_short2(15)==0)) THEN
              CALL DOY_VUELTA()
              CALL DOY_VUELTA_KEY (key,key_aux)
              interruptor=.true.
           END IF
        END IF
        IF (interruptor==.false.) THEN
           IF ((ms_short2(13)==0).and.(ms_short2(16)/=0)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(13)/=0).and.(ms_short2(16)==0)) THEN
                 CALL DOY_VUELTA()
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
        END IF
        IF (interruptor==.false.) THEN
           IF ((ms_short2(14)==0).and.(ms_short2(17)/=0)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(14)/=0).and.(ms_short2(17)==0)) THEN
                 CALL DOY_VUELTA()
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
        END IF
     END IF
  END IF


  IF (interruptor==.false.) THEN
     IF ((ms_short2(12)==12).and.(ms_short2(15)/=15)) THEN
        interruptor=.true.
     ELSE
        IF ((ms_short2(12)/=12).and.(ms_short2(15)==15)) THEN
           CALL DOY_VUELTA()
           CALL DOY_VUELTA_KEY (key,key_aux)
           interruptor=.true.
        END IF
     END IF
     IF (interruptor==.false.) THEN
        IF ((ms_short2(13)==13).and.(ms_short2(16)/=16)) THEN
           interruptor=.true.
        ELSE
           IF ((ms_short2(13)/=13).and.(ms_short2(16)==16)) THEN
              CALL DOY_VUELTA()
              CALL DOY_VUELTA_KEY (key,key_aux)
              interruptor=.true.
           END IF
        END IF
     END IF
     IF (interruptor==.false.) THEN
        IF ((ms_short2(14)==14).and.(ms_short2(17)/=17)) THEN
           interruptor=.true.
        ELSE
           IF ((ms_short2(14)/=14).and.(ms_short2(17)==17)) THEN
              CALL DOY_VUELTA()
              CALL DOY_VUELTA_KEY (key,key_aux)
              interruptor=.true.
           END IF
        END IF
     END IF
  END IF
  microstate=ms_short2
  key_aux=key

  !bb=0
  !! Cruces entre ellos:
  !  - Un cruce
  !  bb(12:17)=(/ 12,13,14,14,16,17 /)  !SII Elijo este como representante
  !  bb(12:17)=(/ 12,13,14,13,16,17 /)  !SII ---> 12,13,14,14,16,17
  !  bb(12:17)=(/ 12,13,14,15,16,12 /)  !NO
  !bb(12:17)=(/ 12,13,14,15,12,17 /)  !SII
  !  - Dos cruces
  !  bb(12:17)=(/ 12,13,14,14,16,12 /)  !NO
  !  bb(12:17)=(/ 12,13,14,14,12,17 /)  !SII Elijo este como representante
  !  bb(12:17)=(/ 12,13,14,13,16,12 /)  !NO
  !  bb(12:17)=(/ 12,13,14,13,12,17 /)  !SII ---> 12,13,14,14,12,17
  !  - Tres cruces da igual
  !! Cruces entre ellos con la primera capa tambien:
  !  - Un cruce entre ellos y uno con primera capa
  !  bb(12:17)=(/ 5,13,14,14,4,17 /)
  !  bb(12:17)=(/ 5,13,14,14,16,4 /)
  ! 
  ! Corrijo esto:

  IF (ms_short2(16)==12) THEN
     CALL DOY_VUELTA()
     CALL DOY_VUELTA_KEY (key,key_aux)
  END IF
  microstate=ms_short2
  key_aux=key
  IF (ms_short2(15)==13) THEN
     key(13)=key_aux(14)
     key(14)=key_aux(13)
     key_aux(13:14)=key(13:14)
     ms_short2(13)=microstate(14)
     ms_short2(14)=microstate(13)
     microstate(13:14)=ms_short(13:14)
     filtro=.true.
     DO i=12,17
        j=ms_short2(i)
        IF ((j>11).and.(filtro(i)==.true.)) THEN
           ms_short2(i)=i
           filtro(i)=.false.
           DO ii=i+1,17
              IF ((microstate(ii)==j).and.(filtro(ii)==.true.)) THEN
                 ms_short2(ii)=i
                 filtro(ii)=.false.
              END IF
           END DO
        END IF
     END DO
  END IF
  microstate=ms_short2
  key_aux=key

  IF (interruptor==.false.) THEN
     IF ((ms_short2(13)==3).and.(ms_short2(16==2))) THEN
        CALL DOY_VUELTA()
        CALL DOY_VUELTA_KEY (key,key_aux)
        microstate=ms_short2
        key_aux=key
     END IF
  END IF


  mss_ind_wat=key
  !!!!!!!

  !IF (0) THEN
  !   interruptor=.true.
  !   DO i=12,17
  !      IF (ms_short2(i)/=bb(i)) THEN
  !         interruptor=.false.
  !         EXIT
  !      END IF
  !   END DO
  !   IF (interruptor==.true.) THEN
  !      print*,'...>',frame,mol
  !      print 117,ms_short2(:)
  !   END IF
  !END IF


  !IF (interruptor==.false.) THEN
  !   print*,'TATE',frame,mol
  !   print 117,ms_short2
  !END IF


117 format (1I3," |",4I3," |",3I3," |",3I3," |",3I3," |",3I3)

END SUBROUTINE REMOVE_PERMUT_SHORT_NOSIMETRIC


SUBROUTINE DOY_VUELTA ()

  IMPLICIT NONE

  INTEGER::i,j,ii
  INTEGER,DIMENSION(17)::microstate
  LOGICAL,DIMENSION(17)::filtro

  filtro=.true.
  microstate=ms_short2

  ms_short2(4)=microstate(5)
  ms_short2(5)=microstate(4)
  ms_short2(12:14)=microstate(15:17)
  ms_short2(15:17)=microstate(12:14)
  microstate=ms_short2
  
  DO i=1,17
     j=ms_short2(i)
     IF ((j>1).and.(filtro(i)==.true.)) THEN
        ms_short2(i)=i
        filtro(i)=.false.
        DO ii=i+1,17
           IF ((microstate(ii)==j).and.(filtro(ii)==.true.)) THEN
              ms_short2(ii)=i
              filtro(ii)=.false.
           END IF
        END DO
     END IF
  END DO

END SUBROUTINE DOY_VUELTA

SUBROUTINE DOY_VUELTA_KEY (key,key_aux)

  IMPLICIT NONE


  INTEGER,DIMENSION(17),INTENT(INOUT)::key,key_aux

  key_aux=key
  key(4)=key_aux(5)
  key(5)=key_aux(4)
  key(12:14)=key_aux(15:17)
  key(15:17)=key_aux(12:14)
  key_aux=key


END SUBROUTINE DOY_VUELTA_KEY

!##################################################################################
!##################################################################################


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! SKINNER HBONDS
!!!!
!!!! Source: R. Kumar, J. R. Schmidt and J. L. Skinner. J. Chem. Phys. 126, 204107 (2007)


SUBROUTINE HBONDS_SKINNER()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::sk_val,aux_cos,aux_dist

  DO i=1,NW
     DO j=1,2
        gg=0
        DO jj=1,nparts
           g=iarr(i,j,jj)
           aux_cos=dot_product(wat_perp(g,:),vect_norm_htoo(i,j,jj,:))
           aux_dist=darr(i,j,jj)
           CALL SKINNER_PARAMETER_DEFINITION (aux_cos,aux_dist,sk_val)
           IF (sk_val>sk_param) THEN
              gg=gg+1
              num_h2o(i,j)=gg
              h2o(i,j,gg)=g
              strength_h2o(i,j,gg)=sk_val     ! Tomaré N como criterio para eliminar hbonds
           END IF
        END DO
     END DO
  END DO

END SUBROUTINE HBONDS_SKINNER

SUBROUTINE SKINNER_PARAMETER_DEFINITION (aux_cos,aux_dist,sk_val)

  IMPLICIT NONE

  DOUBLE PRECISION,INTENT(INOUT)::aux_cos,aux_dist
  DOUBLE PRECISION,INTENT(OUT)::sk_val

  if (aux_cos>=1.0d0) aux_cos=1.0d0
  if (aux_cos<=-1.0d0) aux_cos=-1.0d0
  aux_cos=acos(aux_cos)
  aux_cos=aux_cos*(90/pi)
  IF (aux_cos>90) THEN
     print*,'aquiii error 3.14',aux_cos
     STOP
  END IF
  IF (aux_cos<0) THEN
     print*,'aquiii error 3.14',aux_cos
     STOP
  END IF

  !! darr is in 10^{-10} m
  sk_val=exp(-aux_dist/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)

END SUBROUTINE SKINNER_PARAMETER_DEFINITION


SUBROUTINE SKINNER_PARAMETER (index_wat_o,index_wat_h,index_h,sk_val)

  INTEGER,INTENT(IN)::index_wat_o,index_wat_h,index_h
  DOUBLE PRECISION,    ALLOCATABLE, DIMENSION(:) :: norm_htoo    !! (index_water,Hi,num_neights)
  DOUBLE PRECISION,    ALLOCATABLE, DIMENSION(:) :: perp,aux_vect
  DOUBLE PRECISION::aux_dist,aux_cos
  DOUBLE PRECISION,INTENT(OUT)::sk_val

  Lbox2=Lbox/2.0d0

  ALLOCATE(norm_htoo(3),aux_vect(3))
  ALLOCATE(perp(3))

  norm_htoo=0.0d0
  perp=0.0d0
  aux_vect=0.0d0
  aux_dist=0.0d0
  aux_cos=0.0d0
  sk_val=0.0d0

  aux_vect=XARR(index_wat_o+1,1,:)-XARR(index_wat_h+1,index_h+1,:)

  CALL PBC (aux_vect)
  aux_dist=sqrt(dot_product(aux_vect,aux_vect))
  norm_htoo=aux_vect/aux_dist

  CALL PERPENDICULAR_WATER(index_wat_o+1,perp)

  aux_cos=dot_product(perp(:),norm_htoo(:))

  CALL SKINNER_PARAMETER_DEFINITION (aux_cos,aux_dist,sk_val)

  DEALLOCATE(norm_htoo,perp,aux_vect)

END SUBROUTINE skinner_parameter




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Distance OH; R(o,h)
!!!!
!!!! Source: V. J. Buch. J. Chem. Phys. 96, 3814-3823 (1992)


SUBROUTINE HBONDS_DIST_OH()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::roh_val,aux_cos,aux_dist

  DO i=1,NW
     DO j=1,2
        gg=0
        DO jj=1,nparts
           g=iarr(i,j,jj)
           roh_val=darr(i,j,jj)
           IF (roh_val<roh_param) THEN
              gg=gg+1
              num_h2o(i,j)=gg
              h2o(i,j,gg)=g
              strength_h2o(i,j,gg)=roh_val     ! Tomaré N como criterio para eliminar hbonds
           END IF
        END DO
     END DO
  END DO

END SUBROUTINE HBONDS_DIST_OH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Distance OO and Angle OOH; R(o,o)-Ang(o,o,h)
!!!!
!!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)


SUBROUTINE HBONDS_DIST_OO_ANG_OOH()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::angooh,aux_val,aux_cos

  DOUBLE PRECISION,DIMENSION(3)::aux,aux2
  LOGICAL,DIMENSION(:,:),ALLOCATABLE::filter
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::roo_val
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::roo_vect_norm

  ALLOCATE(filter(NW,NW),roo_val(NW,NW),roo_vect_norm(NW,NW,3))
  filter=.false.
  

  DO i=1,NW
     aux=XARR(i,1,:)
     DO j=1,2
        DO jj=1,nparts
           g=iarr(i,j,jj)
           IF (filter(i,g)==.false.) THEN
              filter(i,g)=.true.
              filter(g,i)=.true.
              aux2=XARR(g,1,:)-aux
              CALL PBC (aux2)
              aux_val=sqrt(dot_product(aux2,aux2))
              aux2=aux2/aux_val
              roo_val(i,g)=aux_val
              roo_val(g,i)=aux_val
              roo_vect_norm(i,g,:)=aux2
              roo_vect_norm(g,i,:)=-aux2
           END IF
        END DO
     END DO
  END DO

  DO i=1,NW
     aux=XARR(i,1,:)
     DO j=1,2
        gg=0
        aux2=XARR(i,j+1,:)-aux
        CALL PBC (aux2)
        aux_val=sqrt(dot_product(aux2,aux2))
        aux2=aux2/aux_val
        DO jj=1,nparts
           g=iarr(i,j,jj)
           IF (roo_val(i,g)<roo_param) THEN
              aux_cos=dot_product(aux2,roo_vect_norm(i,g,:))
              IF (aux_cos>cos_angooh_param) THEN
                 gg=gg+1
                 num_h2o(i,j)=gg
                 h2o(i,j,gg)=g
                 strength_h2o(i,j,gg)=roo_val(i,g)     ! Tomaré N como criterio para eliminar hbonds
                 !strength_h2o(i,j,gg)=acos(aux_cos)
              END IF
           END IF
        END DO
     END DO
  END DO

  DEALLOCATE(filter,roo_val,roo_vect_norm)

END SUBROUTINE HBONDS_DIST_OO_ANG_OOH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Donor-Acceptor-Number
!!!!
!!!! Source: A. D. Hammerich, V. J. Buch. J. Chem. Phys. 128, 111101 (2008)


SUBROUTINE HBONDS_DONOR_ACCEPTOR_NUMBER()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::roh_val
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::neighs_oh_dist
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE::neighs_oh
  INTEGER,DIMENSION(:),ALLOCATABLE::neighs_oh_num

  ALLOCATE(neighs_oh(NW,2,2),neighs_oh_dist(NW,2),neighs_oh_num(NW))
  neighs_oh_num=0

  DO i=1,NW
     DO j=1,2
        g=iarr(i,j,1)
        roh_val=darr(i,j,1)
        h=neighs_oh_num(g)
        IF (h<2) THEN
           h=h+1
           neighs_oh_num(g)=h
           neighs_oh_dist(g,h)=roh_val
           neighs_oh(g,h,1)=i
           neighs_oh(g,h,2)=j
        ELSE
           hh=MAXLOC(neighs_oh_dist(g,:),DIM=1)
           IF (roh_val<neighs_oh_dist(g,hh)) THEN
              neighs_oh_dist(g,hh)=roh_val
              neighs_oh(g,hh,1)=i
              neighs_oh(g,hh,2)=j
           END IF
        END IF
     END DO
  END DO

  DO i=1,NW
     IF (neighs_oh_num(i)==2) THEN
        IF (neighs_oh_dist(i,1)>neighs_oh_dist(i,2)) THEN
           roh_val=neighs_oh_dist(i,1)
           ii=neighs_oh(i,1,1)
           jj=neighs_oh(i,1,2)
           neighs_oh_dist(i,1)=neighs_oh_dist(i,2)
           neighs_oh(i,1,:)=neighs_oh(i,2,:)
           neighs_oh_dist(i,2)=roh_val
           neighs_oh(i,2,1)=ii
           neighs_oh(i,2,2)=jj
        END IF
     END IF
  END DO

  DO i=1,NW
     g=neighs_oh_num(i)
     num_o2h(i)=g
     DO j=1,g
        ii=neighs_oh(i,j,1)
        jj=neighs_oh(i,j,2)
        roh_val=neighs_oh_dist(i,j)
        num_h2o(ii,jj)=1
        h2o(ii,jj,1)=i
        strength_h2o(ii,jj,1)=roh_val
        o2h(i,j)=ii
        o2which(i,j)=jj
        strength_o2h(i,j)=roh_val
     END DO
  END DO

  DEALLOCATE(neighs_oh,neighs_oh_dist,neighs_oh_num)

END SUBROUTINE HBONDS_DONOR_ACCEPTOR_NUMBER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Topological
!!!!
!!!! Source: R. H. Henchman and S. J. Irudayam. J. Phys. Chem. B. 114, 16792-16810 (2010)



SUBROUTINE HBONDS_TOPOLOGICAL()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::roh_val
  INTEGER,DIMENSION(2)::reord
  LOGICAL::interr

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::back1,back2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::back3


  reord(1)=2
  reord(2)=1


  DO i=1,NW
     DO j=1,2
        g=iarr(i,j,1)
        roh_val=darr(i,j,1)
        interr=.true.
        DO jj=1,nparts
           IF (iarr(i,reord(j),jj)==g) THEN
              IF (roh_val>darr(i,reord(j),jj)) interr=.false.
           END IF
        END DO
        IF (interr==.true.) THEN
           DO h=1,2
              DO jj=1,nparts
                 IF (iarr(g,h,jj)==i) THEN
                    IF (roh_val>darr(g,h,jj)) interr=.false.
                 END IF
              END DO
           END DO
        END IF
        IF (interr==.true.) THEN
           num_h2o(i,j)=1
           strength_h2o(i,j,1)=roh_val
           h2o(i,j,1)=g
           h=num_o2h(g)+1
           num_o2h(g)=h
           o2h(g,h)=i
           o2which(g,h)=j
           strength_o2h(g,h)=roh_val
        END IF
     END DO
  END DO

  ALLOCATE(filtro(6),back1(6),back2(6),back3(6))
  filtro=.false.

  DO i=1,NW
     g=num_o2h(i)
     IF (g>1) THEN
        back1=o2h(i,:)
        back2=o2which(i,:)
        back3=strength_o2h(i,:)
        filtro(1:g)=.true.
        DO jj=1,g
           h=MINLOC(back3(:),DIM=1,MASK=filtro)
           o2h(i,jj)=back1(h)
           o2which(i,jj)=back2(h)
           strength_o2h(i,jj)=back3(h)
           filtro(h)=.false.
        END DO
     END IF
  END DO

  DEALLOCATE(filtro,back1,back2,back3)


END SUBROUTINE HBONDS_TOPOLOGICAL



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Donor number and Angle OOH; Donor-Number-Ang(o,o,h)
!!!!
!!!! Source: J. D. Smith, C. D. Cappa, et al. Proc. Natl. Acad. Sci. U.S.A. 102, 14171 (2005).


SUBROUTINE HBONDS_DONOR_ANG_OOH()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::angooh,aux_val,aux_cos,roh_val

  DOUBLE PRECISION,DIMENSION(3)::aux,aux2

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::roo_val
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::roo_vect_norm

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::back1,back2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::back3


  ALLOCATE(roo_val(NW,2),roo_vect_norm(NW,2,3))

  DO i=1,NW
     aux=XARR(i,1,:)
     DO j=1,2
        g=iarr(i,j,1)
        aux2=XARR(g,1,:)-aux
        CALL PBC (aux2)
        aux_val=sqrt(dot_product(aux2,aux2))
        aux2=aux2/aux_val
        roo_val(i,j)=aux_val
        roo_vect_norm(i,j,:)=aux2
     END DO
  END DO

  DO i=1,NW
     aux=XARR(i,1,:)
     DO j=1,2
        aux2=XARR(i,j+1,:)-aux
        CALL PBC (aux2)
        aux_val=sqrt(dot_product(aux2,aux2))
        aux2=aux2/aux_val
        g=iarr(i,j,1)
        aux_cos=dot_product(aux2,roo_vect_norm(i,j,:))
        IF (aux_cos>cos_angooh_param) THEN
           roh_val=darr(i,j,1)
           num_h2o(i,j)=1
           h2o(i,j,1)=g
           strength_h2o(i,j,1)=roh_val     ! Tomaré N como criterio para eliminar hbonds
           !strength_h2o(i,j,1)=acos(aux_cos)
           h=num_o2h(g)+1
           num_o2h(g)=h
           o2h(g,h)=i
           o2which(g,h)=j
           strength_o2h(g,h)=roh_val
        END IF
     END DO
  END DO


  ALLOCATE(filtro(6),back1(6),back2(6),back3(6))
  filtro=.false.

  DO i=1,NW
     g=num_o2h(i)
     IF (g>1) THEN
        back1=o2h(i,:)
        back2=o2which(i,:)
        back3=strength_o2h(i,:)
        filtro(1:g)=.true.
        DO jj=1,g
           h=MINLOC(back3(:),DIM=1,MASK=filtro)
           o2h(i,jj)=back1(h)
           o2which(i,jj)=back2(h)
           strength_o2h(i,jj)=back3(h)
           filtro(h)=.false.
        END DO
     END IF
  END DO

  DEALLOCATE(filtro,back1,back2,back3)
  DEALLOCATE(roo_val,roo_vect_norm)

END SUBROUTINE HBONDS_DONOR_ANG_OOH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Nearest Neighbour
!!!!
!!!! Source: THIS IS NOT A HYDROGEN BOND


SUBROUTINE HBONDS_NEAREST_NEIGHBOUR()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  DOUBLE PRECISION::roh_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::back1,back2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::back3

  DO i=1,NW
     DO j=1,2
        g=iarr(i,j,1)
        roh_val=darr(i,j,1)
        num_h2o(i,j)=1
        h2o(i,j,1)=g
        strength_h2o(i,j,1)=roh_val     ! Tomaré N como criterio para eliminar hbonds
        h=num_o2h(g)+1
        num_o2h(g)=h
        o2h(g,h)=i
        o2which(g,h)=j
        strength_o2h(g,h)=roh_val
     END DO
  END DO


  ALLOCATE(filtro(6),back1(6),back2(6),back3(6))
  filtro=.false.

  DO i=1,NW
     g=num_o2h(i)
     IF (g>1) THEN
        back1=o2h(i,:)
        back2=o2which(i,:)
        back3=strength_o2h(i,:)
        filtro(1:g)=.true.
        DO jj=1,g
           h=MINLOC(back3(:),DIM=1,MASK=filtro)
           o2h(i,jj)=back1(h)
           o2which(i,jj)=back2(h)
           strength_o2h(i,jj)=back3(h)
           filtro(h)=.false.
        END DO
     END IF
  END DO

  DEALLOCATE(filtro,back1,back2,back3)


END SUBROUTINE HBONDS_NEAREST_NEIGHBOUR





END MODULE WAT

!!!! f2py --f90flags='-fast' -c -m pyn_fort_water pyn_fort_water.f90
