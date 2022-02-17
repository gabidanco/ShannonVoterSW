! This program computes the Shannon entropy measure of the
! stationary probability distribution of the voter model
! dynamics with external influencers in small world networks.
! N1 and N0 influence all network nodes homogeneously.
! The relation between t and Temp may need to be adjusted for
! finding the entropy maxima with better precision.
!
! Marcus A.M. de Aguiar and Gabriella D. Franco - Last version - 10/02/2022

program shannon_voter_sw2kd

IMPLICIT REAL*8 (A-H,O-Z)
REAL*8 aux,aux2,anr,annb,prba,auxi,annbi
REAL*8 ai,p, n0, n1, Temp, alog2
INTEGER, ALLOCATABLE :: x(:),a(:,:),nn(:),b(:,:),pr(:)
INTEGER  i,j,k,h,u,ik,imk,tot,kn,l, newneigh, rew, knj
INTEGER  n,nt,nb,nnb,ind,t,nmax, icheck
INTEGER :: iseed(12)
CHARACTER*70 filename
CHARACTER*3 stu, stprew, sttemp

! n0 = number of frozen nodes 0
! n1 = number of frozen nodes 1

OPEN(UNIT=50,FILE='seed.in',STATUS='OLD')
READ(50,*) iseed
CLOSE(50)
CALL RANDOM_SEED(put=iseed)


OPEN(UNIT=7,FILE='input.in',STATUS='OLD',POSITION='REWIND')
	!model parameters
	READ(7,*)n       !population size
	READ(7,*)kn      !network degree = 2*kn
	READ(7,*)nt      !equilibrium time
	READ(7,*)nb      !number of measures to acquire (lines in the output file)
CLOSE(7)

an = dfloat(n)
anb = dfloat(nb)
nmax = nt + 10000*nb

ALLOCATE (x(n),pr(0:n),a(n,n),b(n,n),nn(n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Generate a regular ring network
!     a(i,j) = adjacency matrix
!     b(i,j) = j-th neighbor of node i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! rewire with probability prew
DO rew = 0,10
    prew  = dfloat(rew)*0.01d0
    print *,'prew= ' ,prew
    CALL NUMBSTR(3, rew, stprew)

    DO u=1,10
        CALL NUMBSTR(3, u, stu)
        a = 0
        do i=1,n
            do j=1,kn
                ik = i + j
                imk = i - j
                if(ik > n) ik = ik - n
                if(imk < 1) imk = imk + n
                a(i,ik) = 1
                a(ik,i) = 1
                a(i,imk) = 1
                a(imk,i) = 1
            end do
        end do

        do i=1,n
            do j=i+1,n
                if(a(i,j) /= 0) then
                    call random_number(aux)
                    if(aux < prew) then
                        icheck = 1
                        do while (icheck == 1)
                            call random_number(aux)
                            newneigh = int(aux*n) + 1
                            call random_number(aux)
                            if(aux < 0.5) then
                                ! replace j
                                ! check if newneigh is already a neigbor of i
                                if (newneigh /= i .and. a(i,newneigh) /= 1) then
                                    a(i,j) = 0
                                    a(j,i) = 0
                                    a(i,newneigh) = 1
                                    a(newneigh,i) = 1
                                    icheck = 0
                                else
                                    icheck = 1
                                end if
                            else
                                ! replace i
                                ! check if newneigh is already a neigbor of j
                                if (newneigh /= j .and. a(j,newneigh) /= 1) then
                                    a(i,j) = 0
                                    a(j,i) = 0
                                    a(newneigh,j) = 1
                                    a(j,newneigh) = 1
                                    icheck = 0
                                else
                                    icheck = 1
                                end if
                            end if
                        end do
                    end if
                end if
            end do
        end do

        ! write b(i,j)
        b = 0
        nn = 0
        do i=1,n
            do j=1,n
                if(a(i,j) /= 0) then
                    nn(i) = nn(i) + 1
                    b(i,nn(i)) = j
                end if
            end do
        end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! Realizations
        DO t=0,100
            Temp = dfloat(t)*0.026d0/100.0d0 + 0.024d0
            n0 = Temp
            n1 = Temp

            ! INITIAL CONFIGURATION
            S = 0.0d0
            pr = 0
            prba = 0.0d0
            x = 0 ! network state
            ! dynamical probability
            p = 0 ! Implies that there is dynamics in every time step

            DO i=1,n
                call random_number(aux)
                if(aux < 0.5) x(i) = 1
            END DO

            nf = 1
            ! Dynamics
            DO i=1,nmax
                if(i > nt) nf = mod(i,10000)
                call random_number(aux)
                j = int(aux*n) + 1
                call random_number(aux)

                knj = nn(j)

                annb = dfloat(knj)
                if(aux > p) then
                    call random_number(aux)
                    aux2 = aux*(annb+n0+n1)
                    if(aux2 > annb + n0) then
                        x(j) = 1
                    else if(aux2 > annb) then
                        x(j) = 0
                    else
                        k = int(aux2) + 1
                        x(j) = x(b(j,k))
                    end if
                end if
                ! Testing if equilibrium is reached (given parameter) and calculating the magnetization
                if(nf == 0) then
                    tot = sum(x)  ! count number of states 1
                    pr(tot) = pr(tot) + 1
                end if
            END DO

            alog2 = DLOG(2.0d0)

            filename = 'SW_Shannon_'//trim(stprew)//'_'//trim(stu)//'.dat'
            OPEN(UNIT=58,FILE=filename,STATUS='UNKNOWN')
            DO i=0,n
                prba = dfloat(pr(i))/anb
                if(prba/=0) then
                    S = S + prba*DLOG(prba)
                end if
            END DO
            S = S/alog2

            write(58,*) Temp, -S
            print *, 'S=', -S
            print *,'t= ',t
            print *, 'u=', u


        END DO


        CLOSE(58)

    END DO


END DO

CALL RANDOM_SEED(get=iseed)
OPEN(UNIT=50,FILE='seed.in',STATUS='OLD', POSITION='REWIND')
WRITE (50,*) iseed
close(50)

end program shannon_voter_sw2kd


FUNCTION seq(i,j,l) result(k)
!  INTEGER, intent(in) :: i,j
!  INTEGER :: k
  k = l*(i-1)+j
END FUNCTION seq


FUNCTION line(k,l) result(i)
!  INTEGER, intent(in) :: k
!  INTEGER :: i
  i = INT((k-1)/l)+1
END FUNCTION line

FUNCTION column(k,l) result(j)
!  INTEGER, intent(in) :: k
!  INTEGER :: j
  j = k-l*INT((k-1)/l)
END FUNCTION column


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE NUMBSTR(ID,NUMBER,STR)
CHARACTER*(*) STR
INTEGER*4 ID,NUMBER
CHARACTER*1 B
INTEGER*4 IA0,NN,II,IT
IA0 = ICHAR('0')
NN = NUMBER
DO II=1,ID
J = ID + 1 - II
IT = MOD(NN,10)
B = CHAR(IA0 + IT)
STR(J:J) = B
NN = (NN - IT)/10
END DO
RETURN
END


