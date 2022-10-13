      subroutine ENDFL(I)

      implicit none
      integer, intent(in) :: I

      do
         READ(I,*,END=2)
      enddo
2     continue
      BACKSPACE I

      return
      end subroutine ENDFL
