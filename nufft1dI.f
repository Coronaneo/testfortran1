	subroutine nufft1dIapp(nj,c,U,V,xsub,ns,iflag,r,S)
	implicit none
	integer  r,i,j,k,nj,ns,iflag
	integer xsub(nj)
	complex*16 M(nj,r),N(ns,r),S(ns),c(nj),U(ns,r),V(nj,r)
	double complex in1, out1
	dimension in1(nj), out1(nj)
	integer*8 plan


	M=0


	do i = 1,nj
	   do k = 1,r
              j=xsub(i)
	      M(j,k) = M(j,k)+conjg(V(i,k))*c(i)
	   enddo
	enddo

	!if (iflag < 0) 


	call dfftw_plan_dft_1d(plan,nj,in1,out1,-1,0)
	  do i = 1,r
	     in1 = conjg(M(:,i))
	     call dfftw_execute_dft(plan, in1, out1)
	     N(:,i) = conjg(out1)
	  enddo


	!end
        N=U*N
	S = sum(N,2)


	end
