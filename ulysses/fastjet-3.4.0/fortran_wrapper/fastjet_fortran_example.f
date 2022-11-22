C     example program to run siscone and/or pp sequential recombination
C     algorithms from f77
C     
C     To compile, first make sure that the installation bin directory is
C     in your path (so as to have access to fastjet-config) and then
C     type make -f Makefile.alt fastjet_fortran_example
C     
C     Given the complications inherent in mixing C++ and fortran, your
C     mileage may vary...
C     
C     To use, type: ./fastjet_fortran_example < ../example/data/single-event.dat
C     
C     $Id$
C     
      program fastjet_fortran_example
      implicit none
c ... maximum number of particles
      integer n
      parameter (n = 1000)
      integer i,j
c ... momenta: first index is Lorentz index (1=px,2=py,3=pz,4=E),
c ... second index indicates which particle it is 
c ... [note, indices are inverted relative to convention in Pythia]
      double precision p(4,n)
c ... parameters of the jet algorithm
      double precision  R, f, palg    
c ... array to store the returned jets
      double precision jets(4,n)
      double precision fastjetdmerge,fastjetarea
      integer constituents(n)
      integer npart, njets, nconst ! <= n
      double precision ghost_maxrap, ghost_area
      double precision rapmin,rapmax,phimin,phimax,rho,sigma,meanarea
      integer nrepeat
c ... fill in p (NB, energy is p(4,i))
      do i=1,n
         read(*,*,end=500) p(1,i),p(2,i),p(3,i),p(4,i)
      enddo
      
 500  npart = i-1

      R = 0.6d0
      f = 0.75d0
c.....run the clustering with SISCone
c      call fastjetsiscone(p,npart,R,f,jets,njets)   ! ... now you have the jets
c.....or with a pp generalised-kt sequential recombination alg
      palg = 1d0 ! 1.0d0 = kt, 0.0d0 = Cam/Aachen, -1.0d0 = anti-kt
c      call fastjetppgenkt(p,npart,R,palg,jets,njets)   ! ... now you have the jets

c.....the same, but calculating area information too 
c.....(uselessy slower if you do not need areas)
      ghost_maxrap = 6.0d0 ! make sure you define this as a double precision (with the d0)
      nrepeat = 1
      ghost_area = 0.01d0 ! make sure you define this as a double precision (with the d0)
      call fastjetppgenktwitharea(p,npart,R,palg,
     #                            ghost_maxrap,nrepeat,ghost_area,
     #                            jets,njets)   ! ... now you have the jets


c.....write out all inclusive jets, in order of decreasing pt
      write(*,*) '      px         py          pz         E         pT  
     #        area'
      do i=1,njets
         write(*,*) i,(jets(j,i),j=1,4), sqrt(jets(1,i)**2+jets(2,i)**2)
     #      , fastjetarea(i)
      enddo
      
c.....write out indices of constituents of first jet
      write(*,*)
      write(*,*) 'Indices of constituents of first jet'
      i = 1;
      call fastjetconstituents(i, constituents, nconst)
      write(*,*) (constituents(i),i=1,nconst)

c.....write out the last 5 dmerge values
      write(*,*)
      write(*,*) "dmerge values from last 5 steps"
      do i=0,4
         write(*,*) " dmerge from ",i+1," to ",i," = ", fastjetdmerge(i)
      end do


c.....write out the values of rho, sigma and mean_area in the event
      write(*,*)
      write(*,*) "Background determination"
      rapmin = -3d0
      rapmax = 3d0
      phimin = 0d0
      phimax = 8d0*datan(1d0) ! 2pi
      call fastjetglobalrhoandsigma(rapmin,rapmax,phimin,phimax,
     #                              rho,sigma,meanarea)
      write(*,*) " rho       = ", rho
      write(*,*) " sigma     = ", sigma
      write(*,*) " mean area = ", meanarea
      end
      
