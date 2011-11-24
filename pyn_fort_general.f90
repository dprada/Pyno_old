MODULE aux_funcs_general

CONTAINS

SUBROUTINE dist (pbc_opt,coors1,box1,coors2,n1,n2,matrix)

IMPLICIT NONE

INTEGER,INTENT(IN)::pbc_opt
integer,intent(in)::n1,n2
real,dimension(n1,3),intent(in)::coors1
REAL,DIMENSION(3,3),INTENT(IN)::box1
real,dimension(n2,3),intent(in)::coors2
real,dimension(n1,n2),intent(out)::matrix
integer::i,j
real,dimension(3)::vect


IF (pbc_opt==1) THEN

   vect=0.0d0
   do i=1,n1
      do j=1,n2
         
         vect=(coors1(i,:)-coors2(j,:))
         CALL PBC (vect,box1)
         matrix(i,j)=sqrt(dot_product(vect,vect))
         
      end do
   end do
   
ELSE

   vect=0.0d0
   do i=1,n1
      do j=1,n2
         vect=(coors1(i,:)-coors2(j,:))
         matrix(i,j)=sqrt(dot_product(vect,vect))
      end do
   end do
   
END IF

END SUBROUTINE dist


SUBROUTINE neighbs_limit(pbc_opt,ident,limit,coors1,box1,coors2,n_atoms1,n_atoms2,neighb_list,neighb_dist,neighb_uvect)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::pbc_opt,ident
  INTEGER,INTENT(IN)::limit
  INTEGER,INTENT(IN)::n_atoms1,n_atoms2
  REAL,DIMENSION(n_atoms1,3),INTENT(IN)::coors1
  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
  REAL,DIMENSION(3,3),INTENT(IN)::box1
  INTEGER,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_list
  REAL,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_dist
  REAL,DIMENSION(n_atoms1,limit,3),INTENT(OUT)::neighb_uvect

  LOGICAL::lpbc,lident
  INTEGER::i,j,g
  INTEGER,DIMENSION(:),ALLOCATABLE::list
  REAL,DIMENSION(:),ALLOCATABLE::list_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
  REAL::norm

  lident=.FALSE.
  lpbc=.FALSE.
  IF (ident>0) lident=.TRUE.
  IF (pbc_opt>0) lpbc=.TRUE.


  ALLOCATE(list(limit),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))

  list=0
  list_dists=0.0d0
  list_vects=0.0d0
  aux=0.0d0
  aux2=0.0d0
  filter=.true.

  DO i=1,n_atoms1
     aux=coors1(i,:)
     DO j=1,n_atoms2
        aux2=coors2(j,:)-aux
        IF (lpbc.eqv..true.) CALL PBC (aux2,box1)
        list_dists(j)=sqrt(dot_product(aux2,aux2))
        list_vects(j,:)=aux2
     END DO

     IF (lident.eqv..true.) filter(i)=.false.

     DO j=1,limit
        g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
        list(j)=g
        norm=list_dists(g)
        neighb_dist(i,j)=norm
        neighb_uvect(i,j,:)=list_vects(g,:)/norm
        neighb_list(i,j)=g
        filter(g)=.false.
     END DO

     DO j=1,limit
        g=list(j)
        filter(g)=.true.
     END DO
     filter(i)=.true.
  END DO

  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)

  neighb_list=neighb_list-1

END SUBROUTINE NEIGHBS_LIMIT

SUBROUTINE neighbs_dist2(pbc_opt,ident,ii,dist,coors1,box1,coors2,n_atoms2,neighb_list,neighb_dist,neighb_uvect)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::pbc_opt,ident
  REAL,INTENT(IN)::dist
  INTEGER,INTENT(IN)::n_atoms2,ii
  REAL,DIMENSION(3),INTENT(IN)::coors1
  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
  REAL,DIMENSION(3,3),INTENT(IN)::box1

  INTEGER,DIMENSION(n_atoms2),INTENT(OUT)::neighb_list
  REAL,DIMENSION(n_atoms2),INTENT(OUT)::neighb_dist
  REAL,DIMENSION(n_atoms2,3),INTENT(OUT)::neighb_uvect

  LOGICAL::lpbc,lident
  INTEGER::j,g,limit
  INTEGER,DIMENSION(:),ALLOCATABLE::list
  REAL,DIMENSION(:),ALLOCATABLE::list_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
  REAL::norm

  lident=.FALSE.
  lpbc=.FALSE.
  IF (ident>0) lident=.TRUE.
  IF (pbc_opt>0) lpbc=.TRUE.

  

  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))

  list=0
  list_dists=0.0d0
  list_vects=0.0d0
  aux=0.0d0
  aux2=0.0d0
  filter=.false.


  aux=coors1(:)
  limit=0

  DO j=1,n_atoms2
     aux2=coors2(j,:)-aux
     IF (lpbc.eqv..true.) CALL PBC (aux2,box1)
     norm=sqrt(dot_product(aux2,aux2))
     IF (norm<=dist) THEN
        limit=limit+1
        filter(j)=.true.
        list_dists(j)=norm
        list_vects(j,:)=aux2
     END IF
  END DO
  
  IF (lident.eqv..true.) THEN 
     filter(ii)=.false.
     limit=limit-1
  END IF
  
  
  DO j=1,limit
     g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
     list(j)=g
     norm=list_dists(g)
     neighb_dist(j)=norm
     neighb_uvect(j,:)=list_vects(g,:)/norm
     neighb_list(j)=g
     filter(g)=.false.
  END DO
  
  DO j=1,limit
     g=list(j)
     filter(g)=.false.
  END DO
  filter(ii)=.false.


  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)

  neighb_list=neighb_list-1


END SUBROUTINE NEIGHBS_DIST2


SUBROUTINE neighbs_dist(pbc_opt,ident,dist,coors1,box1,coors2,n_atoms1,n_atoms2,neighb_list,neighb_dist,neighb_uvect)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::pbc_opt,ident
  REAL,INTENT(IN)::dist
  INTEGER,INTENT(IN)::n_atoms1,n_atoms2
  REAL,DIMENSION(n_atoms1,3),INTENT(IN)::coors1
  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
  REAL,DIMENSION(3,3),INTENT(IN)::box1
  INTEGER,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_list
  REAL,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_dist
  REAL,DIMENSION(n_atoms1,n_atoms2,3),INTENT(OUT)::neighb_uvect

  LOGICAL::lpbc,lident
  INTEGER::i,j,g,limit
  INTEGER,DIMENSION(:),ALLOCATABLE::list
  REAL,DIMENSION(:),ALLOCATABLE::list_dists
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
  REAL::norm

  lident=.FALSE.
  lpbc=.FALSE.
  IF (ident>0) lident=.TRUE.
  IF (pbc_opt>0) lpbc=.TRUE.


  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))

  list=0
  list_dists=0.0d0
  list_vects=0.0d0
  aux=0.0d0
  aux2=0.0d0
  filter=.false.

  DO i=1,n_atoms1
     aux=coors1(i,:)
     limit=0
     DO j=1,n_atoms2
        aux2=coors2(j,:)-aux
        IF (lpbc.eqv..true.) CALL PBC (aux2,box1)
        norm=sqrt(dot_product(aux2,aux2))
        IF (norm<=dist) THEN
           limit=limit+1
           filter(j)=.true.
           list_dists(j)=norm
           list_vects(j,:)=aux2
        END IF
     END DO

     IF (lident.eqv..true.) THEN 
        filter(i)=.false.
        limit=limit-1
     END IF

     DO j=1,limit
        g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
        list(j)=g
        norm=list_dists(g)
        neighb_dist(i,j)=norm
        neighb_uvect(i,j,:)=list_vects(g,:)/norm
        neighb_list(i,j)=g
        filter(g)=.false.
     END DO

     DO j=1,limit
        g=list(j)
        filter(g)=.false.
     END DO
     filter(i)=.false.
  END DO

  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)

  neighb_list=neighb_list-1

END SUBROUTINE NEIGHBS_DIST


SUBROUTINE PBC(vector,box)

  IMPLICIT NONE

  REAL,DIMENSION(3),INTENT(INOUT)::vector
  REAL,DIMENSION(3,3),INTENT(IN)::box
  INTEGER::i
  REAL::x,L,Lhalf

  DO i=1,3
     L=box(i,i)
     Lhalf=0.50d0*L
     x=vector(i)
     IF (abs(x)>Lhalf) THEN
        IF (x>Lhalf) THEN
           x=x-L
        ELSE
           x=x+L
        END IF
        vector(i)=x
     END IF
  END DO
  
END SUBROUTINE PBC


END MODULE aux_funcs_general
