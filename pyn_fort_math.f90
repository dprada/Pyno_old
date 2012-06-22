MODULE stats

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::histo_x,histo_y

  CONTAINS

  SUBROUTINE free_mem ()
    IF (ALLOCATED(histo_x)) DEALLOCATE(histo_x)
    IF (ALLOCATED(histo_y)) DEALLOCATE(histo_y)
  END SUBROUTINE free_mem

  SUBROUTINE average (idatos,l,aa,sigma)
    
    implicit none
    INTEGER,INTENT(IN)::l
    DOUBLE PRECISION,DIMENSION(l),INTENT(IN)::idatos
    DOUBLE PRECISION,INTENT(OUT)::aa,sigma

    INTEGER::i,j
    DOUBLE PRECISION::aux_ave,aux_ave2

    aux_ave=0.0d0
    aux_ave2=0.0d0

    DO i=1,l
       aux_ave=aux_ave+idatos(i)
       aux_ave2=aux_ave2+idatos(i)**2
    END DO

    aux_ave=aux_ave/(l*1.0d0)
    aux_ave2=aux_ave2/(l*1.0d0)

    aa=aux_ave
    sigma=sqrt(aux_ave2-aux_ave**2)

  END SUBROUTINE average
  


  SUBROUTINE histograma (opt_norm,opt_prec,opt_range,opt,idatos,ibins,imin_x,imax_x,idelta_x,iprec,l)
    
    implicit none
    INTEGER,INTENT(IN)::opt_norm,opt_prec,opt_range,opt,l,ibins
    DOUBLE PRECISION,DIMENSION(l),INTENT(IN)::idatos
    DOUBLE PRECISION,INTENT(IN)::idelta_x,imax_x,imin_x,iprec

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::frecuencias
    INTEGER::i,j,k,prec,aux,bins
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::datos
    DOUBLE PRECISION::max,min,delta_x,total,sobra
    
    bins=ibins

    ALLOCATE(datos(l))
    datos=0.0d0

    !! Por un problema que hay con python 2.6 corrijo la precision
    IF (opt_prec==1) THEN
       prec=nint(1.0/iprec)
       datos=0.0d0
       DO i=1,l
          aux=nint(idatos(i)*prec)
          datos(i)=(aux*1.0d0)/(prec*1.0d0)
       END DO
       max=(nint(imax_x*prec)*1.0d0)/(prec*1.0d0)
       min=(nint(imin_x*prec)*1.0d0)/(prec*1.0d0)
       delta_x=idelta_x
    ELSE
       datos=idatos
       max=imax_x
       min=imin_x
       delta_x=idelta_x
    END IF
    
    IF (opt_range==0) THEN
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
          sobra=(bins*delta_x-(max-min))/2.0d0
          bins=bins+1
          min=min-sobra
          max=max+sobra
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
          sobra=delta_x/2.0d0
          min=min-sobra
          max=max+sobra
          bins=bins+1
       END IF
    ELSE
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    END IF
    
    ALLOCATE(frecuencias(bins))
    frecuencias=0.0d0
    
    DO k=1,l
       DO j=bins,1,-1
          IF (datos(k)<(min+j*delta_x)) THEN
             frecuencias(j)=frecuencias(j)+1.0d0
          ELSE
             EXIT
          END IF
       END DO
    END DO
    
    DO j=bins,2,-1
       
       frecuencias(j)=frecuencias(j)-frecuencias(j-1)
       
    END DO

    
    IF (opt_norm==1) THEN
       total=SUM(frecuencias)*delta_x
       frecuencias=frecuencias/total
    END IF
    
    
    !! Output
    ALLOCATE(histo_y(bins),histo_x(bins))
    histo_y=frecuencias
    DO i=1,bins
       histo_x(i)=min+(i*1.0d0-0.50d0)*delta_x
    END DO
    
    
    DEALLOCATE(datos,frecuencias)

  END SUBROUTINE histograma
  

  SUBROUTINE binning (opt_prec,opt_range,opt,idatos,ibins,imin_x,imax_x,idelta_x,iprec,l,tray_bins)
    
    implicit none
    INTEGER,INTENT(IN)::opt_prec,opt_range,opt,l,ibins
    DOUBLE PRECISION,DIMENSION(l),INTENT(IN)::idatos
    DOUBLE PRECISION,INTENT(IN)::idelta_x,imax_x,imin_x,iprec
    INTEGER,DIMENSION(l),INTENT(OUT)::tray_bins

    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::frecuencias
    INTEGER::i,j,k,prec,aux,bins
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::datos
    DOUBLE PRECISION::max,min,delta_x,total,sobra
    
    tray_bins=0
    bins=ibins

    ALLOCATE(datos(l))
    datos=0.0d0

    !! Por un problema que hay con python 2.6 corrijo la precision
    IF (opt_prec==1) THEN
       prec=nint(1.0/iprec)
       datos=0.0d0
       DO i=1,l
          aux=nint(idatos(i)*prec)
          datos(i)=(aux*1.0d0)/(prec*1.0d0)
       END DO
       max=(nint(imax_x*prec)*1.0d0)/(prec*1.0d0)
       min=(nint(imin_x*prec)*1.0d0)/(prec*1.0d0)
       delta_x=idelta_x
    ELSE
       datos=idatos
       max=imax_x
       min=imin_x
       delta_x=idelta_x
    END IF
    
    IF (opt_range==0) THEN
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
          sobra=(bins*delta_x-(max-min))/2.0d0
          bins=bins+1
          min=min-sobra
          max=max+sobra
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
          sobra=delta_x/2.0d0
          min=min-sobra
          max=max+sobra
          bins=bins+1
       END IF
    ELSE
       IF (opt==1) THEN
          bins=CEILING((max-min)/delta_x)
       ELSE
          delta_x=(max-min)/(bins*1.0d0)
       END IF
    END IF
    
    DO k=1,l
       DO j=bins,1,-1
          IF (((min+(j-1)*delta_x)<=datos(k)).and.(datos(k)<(min+j*delta_x))) THEN
             tray_bins(k)=j-1
             EXIT
          END IF
       END DO
    END DO
    
    
    !! Output
    ALLOCATE(histo_x(bins))

    DO i=1,bins
       histo_x(i)=min+(i*1.0d0-0.50d0)*delta_x
    END DO
    
    
    DEALLOCATE(datos)

  END SUBROUTINE binning
  


END MODULE stats

! f2py -c -m pyn_fort_math pyn_fort_math.f90

