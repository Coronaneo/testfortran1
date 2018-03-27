	program test
	implicit none
        integer i,iflag,xsub(128),ier,num
        integer nj,ns,r
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(12,128),V1(12,128),U2(12,128),V2(12,128)
        real*16 re1(128),re2(128),x(128),xsubb(128),pi
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(128,12),V(128,12),c(128),S(128),re(128),M(128,12)
        complex*16 fk(-64:63)
        real*8 x1(128),eps
        nj=128
        ns=128
        r=12
        iflag=1
        eps=1E-12
        num=1000
        open(unit = 10,file = 'Ur.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui.txt')
        read(10,*) U2
        open(unit = 10,file = 'Vi.txt')
        read(10,*) V2


        re=dcmplx(re1,re2)        
        U=dcmplx(transpose(U1),transpose(U2))
        V=dcmplx(transpose(V1),transpose(V2))
        !print *,V(2,:)
        !print *,U(1,:)
        do i = 1,128
           x(i) = i*pi/8
        enddo
        xsub=mod(floor(x+0.5),ns)+1
        do i = 1,128
           c(i) = exp(dcmplx(0,1)*i/ns)
        enddo
        !print *,c(1:12)
        do i = 1,128
           x1(i) = i*pi*2*pi/(8*nj)
        enddo

        call system_clock(time_begin,countrage,countmax)
        do i=1,num
        call nufft1dIapp(nj,c,U,V,xsub,ns,iflag,r,S)
        enddo
        call system_clock(time_end,countrage,countmax)
        print *,' T_our         = ',(time_end-time_begin)/num

        call system_clock(time_begin,countrage,countmax)
        do i=1,num
        call nufft1d1f90(nj,x1,c,iflag,eps,ns,fk,ier)
        enddo
        call system_clock(time_end,countrage,countmax)
        print *,' T_nyu         = ',(time_end-time_begin)/num
        print *,sum((S-conjg(fk)*nj)*conjg(S-conjg(fk)*nj))
         


	end program


