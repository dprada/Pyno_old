MODULE WAT

INTEGER::switch
INTEGER::nw,natw
INTEGER::nparts,nparts2
REAL::Lbox,Lbox2,sk_param

REAL, DIMENSION(:,:,:), ALLOCATABLE :: XARR    ! posiciones (molecula,atomo,coordenada)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: DARR    !! (index_water,Hi,num_neights)
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IARR    !! (index_water,Hi,num_neights)
REAL,    ALLOCATABLE, DIMENSION(:,:,:,:) :: vect_norm_htoo    !! (index_water,Hi,num_neights)
REAL,    ALLOCATABLE, DIMENSION(:,:) :: wat_perp

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: num_h2o, o2h, o2which
INTEGER, ALLOCATABLE, DIMENSION(:) :: num_o2h
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: h2o
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: strength_h2o
REAL, ALLOCATABLE, DIMENSION(:,:) :: strength_o2h

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: hbsmol,shell_w
INTEGER, ALLOCATABLE, DIMENSION(:) :: ms_short,ms_short2,mss_ind_wat

REAL::pi

PARAMETER (pi=acos(-1.0d0))
PARAMETER (natw=3)
PARAMETER (nparts=8)
PARAMETER (nparts2=nparts*nparts*2*2)

CONTAINS

SUBROUTINE hbonds_skinner (num_waters,lenbox)

  nw=num_waters
  Lbox=lenbox
  Lbox2=Lbox/2.0d0

  ALLOCATE(vect_norm_htoo(nw,2,nparts,3))
  ALLOCATE(wat_perp(nw,3))

  vect_norm_htoo=0.0d0
  wat_perp=0.0d0

  !! Optimization for hbonds

  ALLOCATE(num_h2o(nw,2),num_o2h(nw),h2o(nw,2,6),o2h(nw,6),o2which(nw,6))
  ALLOCATE(strength_o2h(nw,6),strength_h2o(nw,2,6))

  num_h2o=0
  num_o2h=0
  o2h=0
  h2o=0
  o2which=0
  strength_h2o=0.0d0
  strength_o2h=0.0d0

  IF (switch==0) THEN
     CALL DMAT()
     switch=1
  ELSE
     CALL DMAT_EF()
  END IF

  CALL HBONDS_SKINNER_INTERNAL()

  DEALLOCATE(vect_norm_htoo,wat_perp)

END SUBROUTINE hbonds_skinner

SUBROUTINE skinner_parameter (index_wat_o,index_wat_h,index_h,lenbox,Nval)

!  IMPLICIT NONE

  INTEGER,INTENT(IN)::index_wat_o,index_wat_h,index_h
  REAL, INTENT(IN)::lenbox
  REAL,    ALLOCATABLE, DIMENSION(:) :: norm_htoo    !! (index_water,Hi,num_neights)
  REAL,    ALLOCATABLE, DIMENSION(:) :: perp,aux_vect
  REAL::dd,aux
  REAL,INTENT(OUT)::Nval


  Lbox=lenbox
  Lbox2=Lbox/2.0d0

!  ALLOCATE(vect_norm_htoo(nw,2,nparts,3))
!  ALLOCATE(wat_perp(nw,3))

  ALLOCATE(norm_htoo(3),aux_vect(3))
  ALLOCATE(perp(3))

  norm_htoo=0.0d0
  perp=0.0d0
  aux_vect=0.0d0
  dd=0.0d0
  aux=0.0d0
  Nval=0.0d0

  aux_vect=XARR(index_wat_o+1,1,:)-XARR(index_wat_h+1,index_h+1,:)

  CALL PBC (aux_vect)
  dd=sqrt(dot_product(aux_vect,aux_vect))
  norm_htoo=aux_vect/dd



  CALL PERPENDICULAR_WATER(index_wat_o+1,perp)


  aux=dot_product(perp(:),norm_htoo(:))
  aux=acos(aux)
  aux=aux*(90/pi)
  IF (aux>90) THEN
     print*,'here error 3.14',aux
     STOP
  END IF
  IF (aux<0) THEN
     print*,'here error 3.14',aux
     STOP
  END IF

  !! darr is in 10^{-10} m
  Nval=exp(-dd/0.3430d0)*(7.10d0-0.050d0*aux+0.000210d0*aux**2)


  DEALLOCATE(norm_htoo,perp,aux_vect)

END SUBROUTINE skinner_parameter



SUBROUTINE free_hbonds()
  DEALLOCATE(num_h2o,num_o2h,h2o,o2h,o2which)
  DEALLOCATE(strength_o2h,strength_h2o)
END SUBROUTINE free_hbonds

SUBROUTINE microstates (num_waters,lenbox,mss)

  IMPLICIT NONE 

  INTEGER,INTENT(IN)::num_waters
  REAL,INTENT(IN)::lenbox

  INTEGER,DIMENSION(num_waters,17),INTENT(OUT)::mss

  INTEGER::i,j,k
  
  nw=num_waters
  Lbox=lenbox
  Lbox2=Lbox/2.0d0



!  ALLOCATE(IARR(nw,2,nparts))    ! indice de los nparts atomos de O primeros vecinos de un H dado (indice molecula, H1o H2, orden primeros vecinos)
!  ALLOCATE(DARR(nw,2,nparts))    ! distancia de los pares O-H de acuerdo con iarr
  ALLOCATE(vect_norm_htoo(nw,2,nparts,3))
  ALLOCATE(wat_perp(nw,3))


!  DARR=0.0d0
!  IARR=0
  vect_norm_htoo=0.0d0
  wat_perp=0.0d0

  !! Optimization for hbonds

  ALLOCATE(num_h2o(nw,2),num_o2h(nw),h2o(nw,2,6),o2h(nw,6),o2which(nw,6))
  ALLOCATE(strength_o2h(nw,6),strength_h2o(nw,2,6))
  ALLOCATE(hbsmol(nw,6))

  num_h2o=0
  num_o2h=0
  o2h=0
  h2o=0
  o2which=0
  strength_h2o=0.0d0
  strength_o2h=0.0d0
  hbsmol=0

  IF (switch==0) THEN
     CALL DMAT()
     switch=1
  ELSE
     CALL DMAT_EF ()
  END IF


  ALLOCATE(shell_w(nw,17),ms_short(17),ms_short2(17),mss_ind_wat(17))     ! numero maximo de moleculas vecinas en las dos capas

  shell_w=0
  ms_short=0
  ms_short2=0
  mss_ind_wat=0

  CALL HBONDS_SKINNER_INTERNAL()

  CALL BUILD_HBONDS_LIMIT_NOSIMETRIC ()
  DO j=1,NW      
     CALL STATE_SHORT_NOSIMETRIC(j)   
     CALL REMOVE_PERMUT_SHORT_NOSIMETRIC(j)
     mss(j,:)=ms_short2(:)
  END DO


  DEALLOCATE(vect_norm_htoo,wat_perp)
  DEALLOCATE(num_h2o,num_o2h,h2o,o2h,o2which)
  DEALLOCATE(strength_o2h,strength_h2o)
  DEALLOCATE(shell_w,ms_short2,ms_short,mss_ind_wat)
  DEALLOCATE(hbsmol)

END SUBROUTINE microstates


SUBROUTINE microstates_ind_wat (num_waters,lenbox,mss)

  IMPLICIT NONE 

  INTEGER,INTENT(IN)::num_waters
  REAL,INTENT(IN)::lenbox

  INTEGER,DIMENSION(num_waters,17),INTENT(OUT)::mss

  INTEGER::i,j,k
  
  nw=num_waters
  Lbox=lenbox
  Lbox2=Lbox/2.0d0



!  ALLOCATE(IARR(nw,2,nparts))    ! indice de los nparts atomos de O primeros vecinos de un H dado (indice molecula, H1o H2, orden primeros vecinos)
!  ALLOCATE(DARR(nw,2,nparts))    ! distancia de los pares O-H de acuerdo con iarr
  ALLOCATE(vect_norm_htoo(nw,2,nparts,3))
  ALLOCATE(wat_perp(nw,3))


!  DARR=0.0d0
!  IARR=0
  vect_norm_htoo=0.0d0
  wat_perp=0.0d0

  !! Optimization for hbonds

  ALLOCATE(num_h2o(nw,2),num_o2h(nw),h2o(nw,2,6),o2h(nw,6),o2which(nw,6))
  ALLOCATE(strength_o2h(nw,6),strength_h2o(nw,2,6))
  ALLOCATE(hbsmol(nw,6))

  num_h2o=0
  num_o2h=0
  o2h=0
  h2o=0
  o2which=0
  strength_h2o=0.0d0
  strength_o2h=0.0d0
  hbsmol=0

  IF (switch==0) THEN
     CALL DMAT()
     switch=1
  ELSE
     CALL DMAT_EF ()
  END IF


  ALLOCATE(shell_w(nw,17),ms_short(17),ms_short2(17),mss_ind_wat(17))     ! numero maximo de moleculas vecinas en las dos capas

  shell_w=0
  ms_short=0
  ms_short2=0
  mss_ind_wat=0

  CALL HBONDS_SKINNER_INTERNAL()

  CALL BUILD_HBONDS_LIMIT_NOSIMETRIC ()
  DO j=1,NW      
     CALL STATE_SHORT_NOSIMETRIC(j)   
     CALL REMOVE_PERMUT_SHORT_NOSIMETRIC(j)
     mss(j,:)=mss_ind_wat(:)
     DO i=1,17
        mss(j,i)=mss(j,i)-1
     END DO
  END DO


  DEALLOCATE(vect_norm_htoo,wat_perp)
  DEALLOCATE(num_h2o,num_o2h,h2o,o2h,o2which)
  DEALLOCATE(strength_o2h,strength_h2o)
  DEALLOCATE(shell_w,ms_short2,ms_short,mss_ind_wat)
  DEALLOCATE(hbsmol)

END SUBROUTINE microstates_ind_wat

!!!!!!!!!!! INTERNAL SUBROUTINES:

SUBROUTINE DMAT()

  IMPLICIT NONE
  INTEGER::i,j,jj,g
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter_Hs
  INTEGER,DIMENSION(:),ALLOCATABLE::lista
  REAL,DIMENSION(:),ALLOCATABLE::dHs
  REAL,DIMENSION(:,:),ALLOCATABLE::vect_aux
  REAL,DIMENSION(3)::aux,aux2
  REAL::norm,val

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
  REAL,DIMENSION(:),ALLOCATABLE::dHs
  REAL,DIMENSION(3)::vect_aux,aux2
  REAL,DIMENSION(:,:),ALLOCATABLE::vect_aux2
  INTEGER,DIMENSION(:),ALLOCATABLE::list,list2
  INTEGER,DIMENSION(:,:,:),ALLOCATABLE::iarr2
  REAL::norm,val

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


SUBROUTINE HBONDS_SKINNER_INTERNAL()

  IMPLICIT NONE
  INTEGER::i,ii,j,jj,g,gg,h,hh
  REAL::aux1,aux2
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro
  INTEGER,DIMENSION(:),ALLOCATABLE::back1,back2
  REAL,DIMENSION(:),ALLOCATABLE::back3

  ALLOCATE(filtro(6),back1(6),back2(6),back3(6))


  CALL ALL_NORM_WATER ()

  num_h2o=0
  num_o2h=0
  h2o=0
  o2h=0
  o2which=0
  strength_h2o=0.0d0
  strength_o2h=0.0d0

  !! Si r<2.13570 A entra para cualquier angulo
  !! Si r>2.3077 A sale para cualquier angulo


  DO i=1,NW
     DO j=1,2
        gg=0
        DO jj=1,nparts
           g=iarr(i,j,jj)
           CALL SKINNER_PARAMETER_INTERNAL(i,j,jj,g,aux1) ! mol H, Hi, vecino, mol O, N val

           IF (aux1>sk_param) THEN
              gg=gg+1
              num_h2o(i,j)=gg
              h2o(i,j,gg)=g
              strength_h2o(i,j,gg)=aux1     ! Tomaré N como criterio para eliminar hbonds
              !lo meto en los oxigenos
              h=num_o2h(g)+1
              num_o2h(g)=h
              o2h(g,h)=i
              o2which(g,h)=j
              strength_o2h(g,h)=aux1
           END IF

        END DO
        !ordeno lo que sale de los hidrogenos
        IF (gg>1) THEN
           back1=h2o(i,j,:)
           back3=strength_h2o(i,j,:)
           filtro=.false.
           filtro(1:gg)=.true.
           DO jj=1,gg
              h=MAXLOC(back3(:),DIM=1,MASK=filtro)
              h2o(i,j,jj)=back1(h)
              strength_h2o(i,j,jj)=back3(h)
              filtro(h)=.false.
           END DO
        END IF

     END DO
  END DO


  !ordeno lo que sale de los oxigenos
  DO i=1,NW
     IF (num_o2h(i)>1) THEN
        gg=num_o2h(i)
        back1=o2h(i,:)
        back2=o2which(i,:)
        back3=strength_o2h(i,:)
        filtro=.false.
        filtro(1:gg)=.true.
        DO jj=1,gg
           h=MAXLOC(back3(:),DIM=1,MASK=filtro)
           o2h(i,jj)=back1(h)
           o2which(i,jj)=back2(h)
           strength_o2h(i,jj)=back3(h)
           filtro(h)=.false.
        END DO
     END IF

  END DO

  DEALLOCATE(filtro,back1,back2,back3)

END SUBROUTINE HBONDS_SKINNER_INTERNAL


SUBROUTINE SKINNER_PARAMETER_INTERNAL (molh,hi,vecino,molo,Nval)

  IMPLICIT NONE

  INTEGER::i
  INTEGER,INTENT(IN)::molh,hi,vecino,molo
  REAL,INTENT(OUT)::Nval
  REAL::aux

  aux=dot_product(wat_perp(molo,:),vect_norm_htoo(molh,hi,vecino,:))
  aux=acos(aux)
  aux=aux*(90/pi)
  IF (aux>90) THEN
     print*,'aquiii error 3.14',aux
     STOP
  END IF
  IF (aux<0) THEN
     print*,'aquiii error 3.14',aux
     STOP
  END IF

  !! darr is in 10^{-10} m
  Nval=exp(-darr(molh,hi,vecino)/0.3430d0)*(7.10d0-0.050d0*aux+0.000210d0*aux**2)


END SUBROUTINE SKINNER_PARAMETER_INTERNAL


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
  
  DO i=4,17
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




SUBROUTINE ALL_NORM_WATER()

  IMPLICIT NONE
  INTEGER::i
  REAL,DIMENSION(3)::norm

  DO i=1,NW
     CALL PERPENDICULAR_WATER(i,norm)
     wat_perp(i,:)=norm
  END DO

END SUBROUTINE ALL_NORM_WATER



SUBROUTINE PBC(vector)

  implicit none
  integer::i
  real,dimension(3),intent(INOUT)::vector
  
  DO i=1,3
     IF (abs(vector(i))>Lbox2) THEN
        IF (vector(i)>Lbox2) THEN
           vector(i)=vector(i)-Lbox
        ELSE
           vector(i)=vector(i)+Lbox
        END IF
     END IF
  END DO
  
END SUBROUTINE PBC

SUBROUTINE PERPENDICULAR_WATER (a,vect)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::a
  REAL,DIMENSION(3),INTENT(OUT)::vect
  REAL,DIMENSION(3)::posh1,posh2

  vect=0.0d0

  posh1=xarr(a,1,:)-xarr(a,2,:)
  posh2=xarr(a,1,:)-xarr(a,3,:)

  CALL PRODUCT_VECT(posh1,posh2,vect)
  CALL NORMALIZE_VECT (vect)

END SUBROUTINE PERPENDICULAR_WATER

SUBROUTINE PRODUCT_VECT(a,b,normal)

  REAL,DIMENSION(3),INTENT(IN)::a,b
  REAL,DIMENSION(3),INTENT(OUT)::normal
  REAL::norm
  
  normal(1)=a(2)*b(3)-a(3)*b(2)
  normal(2)=-a(1)*b(3)+a(3)*b(1)
  normal(3)=a(1)*b(2)-a(2)*b(1)

END SUBROUTINE PRODUCT_VECT

SUBROUTINE NORMALIZE_VECT (a)

  REAL,DIMENSION(3),INTENT(INOUT)::a
  REAL::norm

  norm=sqrt(dot_product(a,a))
  a=a/norm

END SUBROUTINE NORMALIZE_VECT



END MODULE WAT

!!!! f2py --f90flags='-fast' -c -m pyn_fort_water pyn_fort_water.f90
