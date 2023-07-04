subroutine calc_dist_matrix(cell, cell_inverse, num_of_nodes, cartesian_coordinates, dist_matrix)
    implicit none
    real(8), intent(in) :: cell(3,3), cell_inverse(3,3)
    integer, intent(in) :: num_of_nodes
    real(8), intent(in) :: cartesian_coordinates(num_of_nodes, 3)
    real(8), intent(out) :: dist_matrix(num_of_nodes, num_of_nodes)


    integer :: n1, n2
    real(8) :: frac_coords(3, num_of_nodes), delta(3), dist2(3), dist

    call dgemm('N', 'T', 3, num_of_nodes, 3, 1d0, cell_inverse, 3, cartesian_coordinates, num_of_nodes, 0d0, frac_coords, 3)
    frac_coords = mod(frac_coords, 1d0)

    !$OMP PARALLEL DO PRIVATE(n1, n2, delta, dist, dist2) SCHEDULE(DYNAMIC)
    do n1 = 2, num_of_nodes
        do n2 = 1, n1 - 1
            delta = frac_coords(:, n1) - frac_coords(:, n2)
            delta = delta - nint(delta)
            call dgemm('N', 'N', 1, 3, 3, 1d0, delta, 1, cell, 3, 0d0, dist2, 1)
            dist = sqrt(sum(dist2**2))
            if (dist < 0.1) then
                write(6,'(A,I,A,I,A)') "ERROR: distance between atom ", n1, " and ", n2, " are less than 0.1 Angstrom in the unit cell! Please check your input file for overlapping atoms."
            end if
            dist_matrix(n1, n2) = dist
            dist_matrix(n2, n1) = dist
        end do
    end do
    !$OMP END PARALLEL DO
end subroutine calc_dist_matrix
