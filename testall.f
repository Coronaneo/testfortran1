	program test
       ! use  Random_Mod
	implicit none
        integer i,iflag
        integer nj,ns,r
        real*16 U1(128,12),V1(128,12),U2(128,12),V2(128,12)
        real*16 re1(128),re2(128),xsub(128)
        complex*16 U(128,12),V(128,12),c(128),S(128),re(128)
        nj=128
        ns=128
        r=12
        iflag=-1
        open(unit = 10,file = 'Ur.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui.txti')
        read(10,*) U2
        open(unit = 10,file = 'Vi.txt')
        read(10,*) V2
        open(unit = 10,file = 'xsub.txt')
        read(10,*) xsub
        open(unit = 10,file = 'Rer.txt')
        read(10,*) re1
        open(unit = 10,file = 'Rei.txt')
        read(10,*) re2

        re=dcmplx(re1,re2)        
        U=dcmplx(U1,U2)
        V=dcmplx(V1,V2)
        do i = 1,128
           c(i) = exp(dcmplx(0,1)*i/ns)
        enddo
        call nufft1dIapp(nj,c,U,V,xsub,ns,iflag,r,S)
        print *,abs(re-S)


	end program


