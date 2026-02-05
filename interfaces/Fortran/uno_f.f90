!---------------------------------------------
! uno_create_model_f
!---------------------------------------------
function uno_create_model_f(problem_type, number_variables, variables_lower_bounds, &
                            variables_upper_bounds, base_indexing) result(model)
   character(len=*) :: problem_type
   integer(uno_int), value :: number_variables, base_indexing
   real(c_double), dimension(*) :: variables_lower_bounds, variables_upper_bounds
   type(c_ptr) :: model
   character(c_char), allocatable :: problem_type_c(:)
   integer :: i, n
   n = len_trim(problem_type)
   allocate(problem_type_c(n+1))
   do i = 1, n
      problem_type_c(i) = problem_type(i:i)
   end do
   problem_type_c(n+1) = c_null_char
   model = uno_create_model(problem_type_c, number_variables, variables_lower_bounds, variables_upper_bounds, base_indexing)
   deallocate(problem_type_c)
end function uno_create_model_f

!---------------------------------------------
! uno_set_solver_integer_option_f
!---------------------------------------------
function uno_set_solver_integer_option_f(solver, option_name, option_value) result(success)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   integer(uno_int), value :: option_value
   logical(c_bool) :: success
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n
   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   success = uno_set_solver_integer_option(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_set_solver_integer_option_f

!---------------------------------------------
! uno_set_solver_double_option_f
!---------------------------------------------
function uno_set_solver_double_option_f(solver, option_name, option_value) result(success)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   real(c_double), value :: option_value
   logical(c_bool) :: success
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n
   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   success = uno_set_solver_double_option(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_set_solver_double_option_f

!---------------------------------------------
! uno_set_solver_bool_option_f
!---------------------------------------------
function uno_set_solver_bool_option_f(solver, option_name, option_value) result(success)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   logical(c_bool), value :: option_value
   logical(c_bool) :: success
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n
   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   success = uno_set_solver_bool_option(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_set_solver_bool_option_f

!---------------------------------------------
! uno_set_solver_string_option_f
!---------------------------------------------
function uno_set_solver_string_option_f(solver, option_name, option_value) result(success)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   character(len=*) :: option_value
   logical(c_bool) :: success
   character(c_char), allocatable :: option_name_c(:)
   character(c_char), allocatable :: option_value_c(:)
   integer :: i, n
   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   n = len_trim(option_value)
   allocate(option_value_c(n+1))
   do i = 1, n
      option_value_c(i) = option_value(i:i)
   end do
   option_value_c(n+1) = c_null_char
   success = uno_set_solver_string_option(solver, option_name_c, option_value_c)
   deallocate(option_name_c)
   deallocate(option_value_c)
end function uno_set_solver_string_option_f

!---------------------------------------------
! uno_get_solver_option_type_f
!---------------------------------------------
function uno_get_solver_option_type_f(solver, option_name) result(option_type)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   integer(uno_int) :: option_type
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   option_type = uno_get_solver_option_type(solver, option_name_c)
   deallocate(option_name_c)
end function uno_get_solver_option_type_f

!---------------------------------------------
! uno_load_solver_option_file_f
!---------------------------------------------
function uno_load_solver_option_file_f(solver, file_name) result(success)
   type(c_ptr), value :: solver
   character(len=*) :: file_name
   logical(c_bool) :: success
   character(c_char), allocatable :: file_name_c(:)
   integer :: i, n

   n = len_trim(file_name)
   allocate(file_name_c(n+1))
   do i = 1, n
      file_name_c(i) = file_name(i:i)
   end do
   file_name_c(n+1) = c_null_char
   success = uno_load_solver_option_file(solver, file_name_c)
   deallocate(file_name_c)
end function uno_load_solver_option_file_f

!---------------------------------------------
! uno_set_solver_preset_f
!---------------------------------------------
function uno_set_solver_preset_f(solver, preset_name) result(success)
   type(c_ptr), value :: solver
   character(len=*) :: preset_name
   logical(c_bool) :: success
   character(c_char), allocatable :: preset_name_c(:)
   integer :: i, n

   n = len_trim(preset_name)
   allocate(preset_name_c(n+1))
   do i = 1, n
      preset_name_c(i) = preset_name(i:i)
   end do
   preset_name_c(n+1) = c_null_char
   success = uno_set_solver_preset(solver, preset_name_c)
   deallocate(preset_name_c)
end function uno_set_solver_preset_f

!---------------------------------------------
! uno_get_solver_integer_option_f
!---------------------------------------------
function uno_get_solver_integer_option_f(solver, option_name) result(solver_integer_option)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   integer(uno_int) :: solver_integer_option
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   solver_integer_option = uno_get_solver_integer_option(solver, option_name_c)
   deallocate(option_name_c)
end function uno_get_solver_integer_option_f

!---------------------------------------------
! uno_get_solver_double_option_f
!---------------------------------------------
function uno_get_solver_double_option_f(solver, option_name) result(solver_double_option)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   real(c_double) :: solver_double_option
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n
   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   solver_double_option = uno_get_solver_double_option(solver, option_name_c)
   deallocate(option_name_c)
end function uno_get_solver_double_option_f

!---------------------------------------------
! uno_get_solver_bool_option_f
!---------------------------------------------
function uno_get_solver_bool_option_f(solver, option_name) result(solver_bool_option)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   logical(c_bool) :: solver_bool_option
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   solver_bool_option = uno_get_solver_bool_option(solver, option_name_c)
   deallocate(option_name_c)
end function uno_get_solver_bool_option_f

!---------------------------------------------
! uno_get_solver_string_option_f
!---------------------------------------------
function uno_get_solver_string_option_f(solver, option_name) result(solver_string_option)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   character(:), allocatable :: solver_string_option
   character(c_char), allocatable :: option_name_c(:)
   type(c_ptr) :: ptr_solver_string_option_c
   integer :: i, n
   character(c_char), pointer :: solver_string_option_c(:)

   n = len_trim(option_name)
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   ptr_solver_string_option_c = uno_get_solver_string_option(solver, option_name_c)
   call c_f_pointer(ptr_solver_string_option_c, solver_string_option_c, [0])
   n = 0
   do while (solver_string_option_c(n+1) /= c_null_char)
      n = n + 1
   end do
   allocate(character(len=n)::solver_string_option)
   do i = 1, n
      solver_string_option(i:i) = solver_string_option_c(i)
   end do
   deallocate(option_name_c)
end function uno_get_solver_string_option_f
