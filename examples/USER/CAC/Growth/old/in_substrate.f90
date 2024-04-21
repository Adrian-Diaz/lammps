
!****************************************************************************
!FORTRAN code to create a CAC model containing a PbTe substrate
!****************************************************************************

      program lammps

      implicit real*8 (a-h,o-z)

      parameter (max=20000000)

      real,dimension(20) :: xb1,yb1,zb1,xb2,yb2,zb2
	  real,dimension(9) :: cell1,cell2
      integer,dimension(10) :: id
	  integer,dimension(8,3) :: node	  
	  integer,dimension(8)  :: ijk
	  integer,dimension(20)  :: iprint,irest
      dimension xx(max),yy(max),zz(max),x(max),y(max),z(max) 
      character*30 title
      character*10 name(20)
	  
	  real    :: xmax,xmin,ymax,ymin,zmax,zmin,bb,na
	  integer :: i,j,k,ii,jj,kk,l,types,p
      integer :: ntot,naa,ne
	  integer :: nnx=8,nny=12
	  !Thickness of different regions in the substrate
	  !nnz has to be multiples of 3, as every 3 unit cells in z form a periodicity
	  integer :: nnz1=6,nnz2=6,nnz3=3,nnz4=9
	  integer :: iul=4 ! Smallest element size
	  real    :: ra,alloy=0 ! For PbTe-PbSe epialyer
	  real    :: midx,midy,midz,boundaryz=0.5
	  real    :: axx,ayy,azz,cxx,cyy,czz
	  integer,dimension(20):: name1,name2
	  
      open(11,file='Basis.PbTe') !PbTe unit cell	  
	  open(12,file='Basis.PbSe') !PbSe unit cell
	  open(4,file='ovito.txt') !Use the file to visulize the model in OVITO
	  open(2,file='Pb.txt') !CAC data file

!--------------------------------------------------------------
!     read in the information of the cube unit cell for LiF
!--------------------------------------------------------------

      natms=2
	  !PbTe: grain1
	  read(11,*)
	  read(11,*)
      read(11,*)cell1(1),cell1(2),cell1(3)
      read(11,*)cell1(4),cell1(5),cell1(6)
      read(11,*)cell1(7),cell1(8),cell1(9)
	  read(11,*)
      do i=1,natms
      read(11,*)name1(i),xb1(i),yb1(i),zb1(i) 
      enddo 
	  
      ! PbSe: grain2
	  read(12,*)
	  read(12,*)
      read(12,*)cell2(1),cell2(2),cell2(3)
      read(12,*)cell2(4),cell2(5),cell2(6)
      read(12,*)cell2(7),cell2(8),cell2(9)
	  read(12,*)
      do i=1,natms
      read(12,*)name2(i),xb2(i),yb2(i),zb2(i) 
      enddo

	  ! Thermal expansion causes lattice constant increase
	  do i=1,8
	  cell1(i)=cell1(i)*1.0319 !1100K:1.0426; 1000K:1.03677 950K:1.0343 900K:1.0319 800K:1.0274 650K:1.02114 
	  enddo
	  cell1(9)=cell1(9)*1.0319
	  do i=1,natms
	  xb1(i)=xb1(i)*1.0319
	  yb1(i)=yb1(i)*1.0319
	  zb1(i)=zb1(i)*1.0319
	  enddo	   	  	 
!--------------------------------------------------------------
!    node order
!--------------------------------------------------------------
      	 node(1,1)=0  !node 1
         node(1,2)=0
         node(1,3)=0

      	 node(2,1)=1  !node 2
         node(2,2)=0
         node(2,3)=0

      	 node(3,1)=1  !node 3
         node(3,2)=1
         node(3,3)=0

      	 node(4,1)=0  !node 4
         node(4,2)=1
         node(4,3)=0

		 node(5,1)=0  !node 5
         node(5,2)=0
         node(5,3)=1

      	 node(6,1)=1  !node 6
         node(6,2)=0
         node(6,3)=1

      	 node(7,1)=1  !node 7
         node(7,2)=1
         node(7,3)=1

      	 node(8,1)=0  !node 8
         node(8,2)=1
         node(8,3)=1
	 
 
!--------------------------------------------------------------
!    Compute number of nodes and atoms 
!--------------------------------------------------------------

        xmax=0d0
		xmin=0d0
		ymax=0d0
		ymin=0d0
        zmax=0d0
		zmin=0d0

		midx=-(cell1(1)+cell1(4)+cell1(7))/2d0
		midy=-(cell1(2)+cell1(5)+cell1(8))/2d0
		midz=-(cell1(3)+cell1(6)+cell1(9))/2d0 
	   
		!write(*,*)'midx,midy:',midx,midy
		!write(4,*)nd
		!write(4,*)'Atoms.'
		

        !element number        
		 n=0
		! DOF number
		 ntot=0
		! ovito DOF
		 nd=0
		!equivalent atom number
		 na=0
		!atom number 
		 naa=0
		!CG element number
		 ne=0

		!Bttom of PbTe substrate
		iulx=iul*4
		iuly=iul*4
		iulz=3
	     do  k=0,1-1
           do  j=0,nny-1    
             do  i=0,nnx-1 

			 n=n+1
		     ne=ne+1
			 xx(n)=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             yy(n)=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             zz(n)=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 
			 !write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=xx(n)+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=yy(n)+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=zz(n)+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 !write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 !write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
		 
			 end do
		  end do
		
            end do
	       end do
	      end do
		  
		  midz=midz+1.0*iulz*cell1(9)

		!PbTe substrate: large CG elements
		iulx=iul*4
		iuly=iul*4
		iulz=iul*4
	     do  k=0,nnz1-1
           do  j=0,nny-1    
             do  i=0,nnx-1  
			 n=n+1
			 ne=ne+1
	 
			 xx(n)=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             yy(n)=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             zz(n)=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 
			 !write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=xx(n)+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=yy(n)+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=zz(n)+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 !write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 !write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
		 
			 end do
		  end do
		
            end do
	       end do
	      end do
		  
		  midz=midz+nnz1*iulz*cell1(9)
		  
		!PbTe substrate: Mid CG elements
		iulx=iul*2
		iuly=iul*2
		iulz=iul*2
	     do  k=0,nnz2-1
           do  j=0,nny*2-1    !nny=24
             do  i=0,nnx*2-1  !nnx=24

			 n=n+1
			 ne=ne+1
	 
			 xx(n)=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             yy(n)=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             zz(n)=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 
			 !write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=xx(n)+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=yy(n)+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=zz(n)+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 !write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 !write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
			 
			 end do
		  end do
		
            end do
	       end do
	      end do
		  
		  midz=midz+nnz2*iulz*cell1(9)

		!PbTe substrate: Small CG elements
		iulx=iul
		iuly=iul
		iulz=iul
	     do  k=0,nnz3-1
           do  j=0,nny*4-1   
             do  i=0,nnx*4-1  
			
			 cxx=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             cyy=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             czz=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 			 
			 n=n+1
			 ne=ne+1

			 !write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=cxx+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=cyy+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=czz+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 !write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 !write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
			 
			 end do
		  end do
		  
            end do
	       end do
	      end do
		  
		  midz=midz+nnz3*iulz*cell1(9)+(cell1(3)+cell1(6)+cell1(9))/2d0
		  
		!PbTe substrate: atomic region
	     do  k=0,nnz4-1
           do  j=0,nny*iul*4-1    !nny=24
             do  i=0,nnx*iul*4-1  !nnx=24
               do  l=1,natms

			n=n+1
            axx=xb1(l)+cell1(1)*i+cell1(4)*mod(j,2)+cell1(7)*k
            ayy=yb1(l)+cell1(2)*i+cell1(5)*j+cell1(8)*mod(k,3)
            azz=zb1(l)+cell1(3)*i+cell1(6)*j+cell1(9)*k+midz
			na=na+1
			naa=naa+1
			!write(2,'(i8,a15,4i8)')n,'Atom',1,1,1,1

             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if

			dis=dz*1.47-sqrt(axx*axx+ayy*ayy)
			if (dis>0) azz=azz-strain*dis/dz/1.47
            !write(2,'(3i8,f10.2,3f24.8)')1,1,types,charge,axx,ayy,azz
		  
			nd=nd+1
			!write(4,'(2i8,3f24.8)')nd,types,axx,ayy,azz

			 if(azz .ge. zmax) zmax=azz
			 if(azz .le. zmin) zmin=azz			 
				enddo
		     enddo
	       enddo
		 enddo	
		 
		 midz=midz+nnz4*cell1(9)

		!PbSe epilayer: atomic region
		!Set alloy>0, you will have PbSe-PbTe epilayer
		n_Se=0
		n_Te=0
		call random_seed()
		!PbSe atoms
	     do  k=0,0-1
           do  j=0,nny*iul*4-1    !nny=24
             do  i=0,nnx*iul*4-1  !nnx=24
               do  l=1,natms
			
			if(l==1 .or. k==0)then			
			n=n+1
            axx=xb1(l)+cell1(1)*i+cell1(4)*mod(j,2)+cell1(7)*k
            ayy=yb1(l)+cell1(2)*i+cell1(5)*j+cell1(8)*mod(k,3)
            azz=zb2(l)+cell2(3)*i+cell2(6)*j+cell2(9)*k+midz
			na=na+1
			naa=naa+1
			!write(2,'(i8,a15,4i8)')n,'Atom',1,1,1,1

             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=2
             charge=-0.8d0
             end if

			 call random_number(ra)
			 if (ra<alloy .and. types==2)then
				types=3
				n_Te=n_Te+1
			 elseif (ra>=alloy .and. types==2)then
			    n_Se=n_Se+1
			 end if
			 
            !write(2,'(3i8,f10.2,3f24.3)')1,1,types,charge,axx,ayy,azz

			nd=nd+1
			!write(4,'(2i8,3f24.3)')nd,types,axx,ayy,azz
			
			 if(azz .ge. zmax) zmax=azz
			 if(azz .le. zmin) zmin=azz

			end if
				enddo
		     enddo
	       enddo
		 enddo			  
		  		
		  		 		 
		  write(*,*)'substrate size z',midz
		  write(*,*)'substrate size x:',(nnx*iul*4)*cell1(1)
		  write(*,*)'substrate size y:',(nny*iul*4)*cell1(5)
		  write(*,*) 'equivalent atom number:',na
		  write(*,*) 'atom number:',naa
		  write(*,*) 'element number:',ne	
  
!--------------------------------------------------------------
!    title
!--------------------------------------------------------------		  
      write(2,'(a10,/)') 'PbTe/PbSe SL'
	  
	  write(2,'(i10,a30)') n,'cac elements'
	  write(2,'(i10,a30,/)') 3,'atom types'
	  
      write(2,'(2f15.6,2a4)')0,(nnx*iul*4)*cell1(1),'xlo','xhi'
      write(2,'(2f15.6,2a4)')0,(nny*iul*4)*cell1(5),'ylo','yhi'
      write(2,'(2f15.6,2a4)')zmin-boundaryz,zmax+300,'zlo','zhi' 
      write(2,'(a20,/)')'Masses'
	  
	  write(*,*) 'Z boundary:',zmax+600

      write(2,'(i10,f20.6)')1,207.2
      write(2,'(i10,f20.6)')2,78.971
	  write(2,'(i10,f20.6)')3,127.60

      write(2,'(a20,/)')'CAC Elements'		  
		  
!--------------------------------------------------------------
!    coordinates of element nodes and atoms
!--------------------------------------------------------------
		
        xmax=0d0
		xmin=0d0
		ymax=0d0
		ymin=0d0
        zmax=0d0
		zmin=0d0

		midx=-(cell1(1)+cell1(4)+cell1(7))/2d0
		midy=-(cell1(2)+cell1(5)+cell1(8))/2d0
		midz=-(cell1(3)+cell1(6)+cell1(9))/2d0 
	   
		!write(*,*)'midx,midy:',midx,midy
		write(4,*)nd
		write(4,*)'Atoms.'
		

        !element number        
		 n=0
		! DOF number
		 ntot=0
		! ovito DOF
		 nd=0
		!equivalent atom number
		 na=0
		!atom number 
		 naa=0
		!CG element number
		 ne=0

		!Bttom of PbTe substrate
		iulx=iul*4
		iuly=iul*4
		iulz=3
	     do  k=0,1-1
           do  j=0,nny-1    
             do  i=0,nnx-1 

			 n=n+1
		     ne=ne+1
			 xx(n)=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             yy(n)=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             zz(n)=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 
			 write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=xx(n)+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=yy(n)+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=zz(n)+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
		 
			 end do
		  end do
		
            end do
	       end do
	      end do
		  
		  midz=midz+1.0*iulz*cell1(9)

		!PbTe substrate: large CG elements
		iulx=iul*4
		iuly=iul*4
		iulz=iul*4
	     do  k=0,nnz1-1
           do  j=0,nny-1    
             do  i=0,nnx-1  
			 n=n+1
			 ne=ne+1
	 
			 xx(n)=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             yy(n)=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             zz(n)=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 
			 write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=xx(n)+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=yy(n)+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=zz(n)+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
		 
			 end do
		  end do
		
            end do
	       end do
	      end do
		  
		  midz=midz+nnz1*iulz*cell1(9)
		  
		!PbTe substrate: Mid CG elements
		iulx=iul*2
		iuly=iul*2
		iulz=iul*2
	     do  k=0,nnz2-1
           do  j=0,nny*2-1    !nny=24
             do  i=0,nnx*2-1  !nnx=24

			 n=n+1
			 ne=ne+1
	 
			 xx(n)=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             yy(n)=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             zz(n)=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 
			 write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=xx(n)+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=yy(n)+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=zz(n)+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
			 
			 end do
		  end do
		
            end do
	       end do
	      end do
		  
		  midz=midz+nnz2*iulz*cell1(9)

		!PbTe substrate: Small CG elements
		iulx=iul
		iuly=iul
		iulz=iul
	     do  k=0,nnz3-1
           do  j=0,nny*4-1   
             do  i=0,nnx*4-1  
			
			 cxx=iulx*i*cell1(1)+iuly*mod(j,2)*cell1(4)+iulz*k*cell1(7)+midx
             cyy=iulx*i*cell1(2)+iuly*j*cell1(5)+iulz*mod(k,3)*cell1(8)+midy
             czz=iulx*i*cell1(3)+iuly*j*cell1(6)+iulz*k*cell1(9)+midz
			 			 
			 n=n+1
			 ne=ne+1

			 write(2,'(i8,a15,4i8)')n,'Eight_Node',natms,iulx,iuly,iulz
             na=na+iulx*iuly*iulz*natms
			 
        !node circle
		  do l=1,natms
		   	do  m=1,8
			 ntot=ntot+1
			 ii=node(m,1)
			 jj=node(m,2)
			 kk=node(m,3)             
			 x(ntot)=cxx+iulx*ii*cell1(1)+iuly*jj*cell1(4)+iulz*kk*cell1(7)
             y(ntot)=cyy+iulx*ii*cell1(2)+iuly*jj*cell1(5)+iulz*kk*cell1(8)
			 z(ntot)=czz+iulx*ii*cell1(3)+iuly*jj*cell1(6)+iulz*kk*cell1(9)
			 
			 nd=nd+1
			 write(4,'(2i8,3f24.3)') nd,4,x(ntot),y(ntot),z(ntot)
			 
             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if
			 
			 write(2,'(3i8,f10.2,3f24.3)')m,l,types,charge,x(ntot)+xb1(l),y(ntot)+yb1(l),z(ntot)+zb1(l)
			 
			 end do
		  end do
		  
            end do
	       end do
	      end do
		  
		  midz=midz+nnz3*iulz*cell1(9)+(cell1(3)+cell1(6)+cell1(9))/2d0
		  
		!PbTe substrate: atomic region
	     do  k=0,nnz4-1
           do  j=0,nny*iul*4-1    !nny=24
             do  i=0,nnx*iul*4-1  !nnx=24
               do  l=1,natms

			n=n+1
            axx=xb1(l)+cell1(1)*i+cell1(4)*mod(j,2)+cell1(7)*k
            ayy=yb1(l)+cell1(2)*i+cell1(5)*j+cell1(8)*mod(k,3)
            azz=zb1(l)+cell1(3)*i+cell1(6)*j+cell1(9)*k+midz
			na=na+1
			naa=naa+1
			write(2,'(i8,a15,4i8)')n,'Atom',1,1,1,1

             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=3 
             charge=-0.8d0
             end if

			dis=dz*1.47-sqrt(axx*axx+ayy*ayy)
			if (dis>0) azz=azz-strain*dis/dz/1.47
            write(2,'(3i8,f10.2,3f24.8)')1,1,types,charge,axx,ayy,azz
		  
			nd=nd+1
			write(4,'(2i8,3f24.8)')nd,types,axx,ayy,azz

			 if(azz .ge. zmax) zmax=azz
			 if(azz .le. zmin) zmin=azz			 
				enddo
		     enddo
	       enddo
		 enddo	
		 
		 midz=midz+nnz4*cell1(9)

		!PbSe epilayer: atomic region
		!Set alloy>0, you will have PbSe-PbTe epilayer
		n_Se=0
		n_Te=0
		call random_seed()
		!PbSe atoms
	     do  k=0,0-1
           do  j=0,nny*iul*4-1    !nny=24
             do  i=0,nnx*iul*4-1  !nnx=24
               do  l=1,natms
			
			if(l==1 .or. k==0)then			
			n=n+1
            axx=xb1(l)+cell1(1)*i+cell1(4)*mod(j,2)+cell1(7)*k
            ayy=yb1(l)+cell1(2)*i+cell1(5)*j+cell1(8)*mod(k,3)
            azz=zb2(l)+cell2(3)*i+cell2(6)*j+cell2(9)*k+midz
			na=na+1
			naa=naa+1
			write(2,'(i8,a15,4i8)')n,'Atom',1,1,1,1

             if(l<=natms/2) then
			 types=1
			 charge=0.8d0
			 end if
             if(l>natms/2) then
			 types=2
             charge=-0.8d0
             end if

			 call random_number(ra)
			 if (ra<alloy .and. types==2)then
				types=3
				n_Te=n_Te+1
			 elseif (ra>=alloy .and. types==2)then
			    n_Se=n_Se+1
			 end if
			 
            write(2,'(3i8,f10.2,3f24.3)')1,1,types,charge,axx,ayy,azz

			nd=nd+1
			write(4,'(2i8,3f24.3)')nd,types,axx,ayy,azz
			
			 if(azz .ge. zmax) zmax=azz
			 if(azz .le. zmin) zmin=azz

			end if
				enddo
		     enddo
	       enddo
		 enddo			  
		  
		  write(*,*) 'equivalent atom number:',na
		  write(*,*) 'atom number:',naa
		  write(*,*) 'element number:',ne

!--------------------------------------------------------------	  		  
      end		  
