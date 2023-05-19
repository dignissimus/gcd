program main
    print *, naive_gcd(14, 56)
    print *, naive_gcd_half(14, 56)
    print *, erastosthenes(14, 56)
    print *, euclid(14, 56)
    print *, euclidean(14, 56)
    print *, binary(14, 56)
contains

integer function binary(left, right)
    implicit none
    integer, intent(in) :: left, right
    integer :: smaller, larger
    integer :: next
    integer :: k = 0

    smaller = min(left, right)
    larger = max(left, right)

    do while ((mod(smaller, 2) .eq. 0) .and. (mod(larger, 2) .eq. 0))
        smaller = smaller / 2
        larger = larger / 2
        k = k + 1
    end do

    do while (smaller .ne. 0)
        do while (mod(larger, 2) .eq. 0)
            larger = larger / 2
        end do
        do while (mod(smaller, 2) .eq. 0)
            smaller = smaller / 2
        end do

       next = larger - smaller
       larger = max(smaller, next)
       smaller = min(smaller, next)
    end do
    binary = lshift(larger, k)
end function

integer function euclid(left, right)
    implicit none
    integer, intent(in) :: left, right
    integer :: smaller, larger
    integer :: next

    smaller = min(left, right)
    larger = max(left, right)
    do while ((smaller .ne. 0) .and. (smaller .ne. 1))
        next = larger - smaller
        larger = max(smaller, next)
        smaller = min(smaller, next)
    end do
    euclid = larger 
end function

integer function euclidean(left, right)
    implicit none
    integer, intent(in) :: left, right
    integer :: smaller, larger
    integer :: next

    smaller = min(left, right)
    larger = max(left, right)
    do while ((smaller .ne. 0) .and. (smaller .ne. 1))
        next = mod(larger, smaller)
        larger = max(smaller, next)
        smaller = min(smaller, next)
    end do
    euclidean = larger 
end function

integer function erastosthenes(left, right)
    implicit none
    integer, intent(in) :: left, right
    integer :: divisor
    integer :: i
    logical, dimension(:), allocatable :: isnondivisor
    allocate(isnondivisor(left), source = .false.)

    erastosthenes = 0
    if (mod(right, left) .eq. 0) then
        erastosthenes = left
        return
    end if
    do divisor = 1, left / 2
        if (isnondivisor(divisor)) then
            continue
        end if
        if ((mod(left, divisor) .eq. 0) .and. (mod(right, divisor) .eq. 0)) then
            erastosthenes = divisor
        else
            do i = 1, left / divisor
                isnondivisor(i * divisor) = .true.
            end do
        end if
    end do
end function

integer function naive_gcd(left, right)
    integer, intent(in) :: left, right
    integer :: smaller
    integer :: divisor
    logical :: divides_left, divides_right

    naive_gcd = 1
    smaller = min(left, right)

    do divisor = 1, smaller
        divides_left = mod(left, divisor) .eq. 0
        divides_right = mod(right, divisor) .eq. 0
        if (divides_left .and. divides_right) then
            naive_gcd = divisor
        end if
    end do
end function

integer function naive_gcd_half(left, right)
    integer, intent(in) :: left, right
    integer :: smaller, larger
    integer :: divisor
    logical :: smaller_divides_larger
    logical :: divides_left, divides_right

    naive_gcd_half = 1
    smaller = min(left, right)
    larger = max(left, right)

    smaller_divides_larger = mod(larger, smaller) .eq. 0
    if (smaller_divides_larger) then
        naive_gcd_half = smaller
        return
    end if

    do divisor = 1, smaller / 2
        divides_left = mod(left, divisor) .eq. 0
        divides_right = mod(right, divisor) .eq. 0
        if (divides_left .and. divides_right) then
            naive_gcd_half = divisor
        end if
    end do
end function
end program
