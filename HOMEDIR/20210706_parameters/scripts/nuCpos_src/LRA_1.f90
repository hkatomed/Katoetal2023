subroutine LRA_1(inseq_num, logascL, logascR, &
                        &freqL1, tranL1, TtranL2, TtranL3, TtranL4, &
                        &TfreqN4L, TfreqN4R, TtranN4)

  implicit none
  character(147)   inseq
  integer    c(147), w(147), i, j, k, l, n, inseq_num(147)
  real(8)    ascL, ascR
  real(8)    logascL, logascR
  real(8)    TtranL2(16,4),TtranL3(64,4),TtranL4(256,4)
  real(8)    TtranN4((147-4)*256,4), TfreqN4L(64,4), TfreqN4R(64,4)
  real(8)    freqL1(4),tranL1(4,4),tranL2(4,4,4),tranL3(4,4,4,4),tranL4(4,4,4,4,4)
  real(8)    freqN4L(4,4,4,4), freqN4R(4,4,4,4)
  real(8)    tranN4(5:147,4,4,4,4,4),freqL4(4,4,4,4)

  inseq=char(inseq_num(1))
  do i=2,147;
     inseq = inseq(1:(i-1)) // char(inseq_num(i))
  end do

  do i=1,147
    if(inseq(i:i)=='A'.or.inseq(i:i)=='a') then
      w(i) = 1
    elseif(inseq(i:i)=='C'.or.inseq(i:i)=='c') then
      w(i) = 2
    elseif(inseq(i:i)=='G'.or.inseq(i:i)=='g') then
      w(i) = 3
    elseif(inseq(i:i)=='T'.or.inseq(i:i)=='t') then
      w(i) = 4
    else
      return
    endif
  enddo

  do i=1,147
    c(i) = 5_1 - w(147 - i + 1)
  end do

  do i=1,4; do j=1,4; tranL2(i,j,:)=TtranL2((i-1)*4+j,:)
  do k=1,4; tranL3(i,j,k,:)=TtranL3((i-1)*16+(j-1)*4+k,:)
  do l=1,4; tranL4(i,j,k,l,:)=TtranL4((i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4
    freqN4L(i,j,k,:)=TfreqN4L((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4
    freqN4R(i,j,k,:)=TfreqN4R((i-1)*16+(j-1)*4+k,:)
  end do; end do; end do

  do n=5,147; do i=1,4; do j=1,4; do k=1,4; do l=1,4
    tranN4(n,i,j,k,l,:)=TtranN4((n-5)*256+(i-1)*64+(j-1)*16+(k-1)*4+l,:)
  end do; end do; end do; end do; end do

  do i=1,4; do j=1,4; do k=1,4; do l=1,4
    freqL4(i,j,k,l)=freqL1(i)*tranL1(i,j)*tranL2(i,j,k)*tranL3(i,j,k,l)
  end do; end do; end do; end do

! 1:73 (ascL)
  ascL=freqN4L(w(1),w(2),w(3),w(4))&
      &/freqL4(w(1),w(2),w(3),w(4))&
      &*freqN4R(c(75),c(76),c(77),c(78))&
      &/freqL4(c(75),c(76),c(77),c(78))
  do i=5,73
    ascL=ascL*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=79,147
    ascL=ascL*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

! 75:147 (ascR)
  ascR=freqN4R(w(75),w(76),w(77),w(78))&
      &/freqL4(w(75),w(76),w(77),w(78))&
      &*freqN4L(c(1),c(2),c(3),c(4))&
      &/freqL4(c(1),c(2),c(3),c(4))
  do i=79,147
    ascR=ascR*tranN4(i,w(i-4),w(i-3),w(i-2),w(i-1),w(i))&
      &/tranL4(w(i-4),w(i-3),w(i-2),w(i-1),w(i))
  end do
  do i=5,73
    ascR=ascR*tranN4(i,c(i-4),c(i-3),c(i-2),c(i-1),c(i))&
      &/tranL4(c(i-4),c(i-3),c(i-2),c(i-1),c(i))
  end do

  logascL = log(ascL)
  logascR = log(ascR)

end subroutine LRA_1

