!=====================================================
!   Matrix determinant computing for double complex matrixes.
!=====================================================
Subroutine Det(N, Nfill, A, Det1)
    Implicit None
    Double Complex::S, R, Det1, P
    Integer :: Rank, I, K, K1, J
    integer, intent(in) :: N
    integer, intent(in) :: Nfill
    Double Complex, intent(in) :: A(N, N)
    Double Complex, Allocatable, Dimension (:,:) :: Tmp
    Rank= Nfill
    Allocate(Tmp(Rank,Rank))
    Tmp(:,:)=A(1:Rank,1:Rank)
    P=1._8
    Do K=1, Rank-1
        K1=K+1
        S=Tmp(K,K)
        J=K
        Do I=K1,Rank
            R=Tmp(I, K)
            If (Abs(R).GT.Abs(S)) Then
                S=R; J=I
            EndIf
        EndDo
        If (S.EQ.(0._8,0._8)) Return
        If (J.NE.K) Then
            Do I=K, Rank
                R=Tmp(K,I); Tmp(K,I)=Tmp(J,I); Tmp(J,I)=R
            EndDo
            P=-P
        EndIf
        Tmp(K,K1:Rank)=Tmp(K,K1:Rank)/S
        Do I=K1, Rank
            Tmp(I,K1:Rank)=Tmp(I,K1:Rank)-Tmp(K,K1:Rank)*Tmp(I,K)
        EndDo
        P=P*S
    EndDo
    S=P*Tmp(Rank,Rank)
    Det1=S
    Deallocate(Tmp)
Return
End Subroutine Det

