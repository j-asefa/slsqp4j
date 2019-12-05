subroutine square_cube(i, isquare, icube)
    integer, intent(in)  :: i              ! input
    integer, intent(out) :: isquare, icube ! output
    isquare = i**2
    icube   = i**3
end subroutine square_cube