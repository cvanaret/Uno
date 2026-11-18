! Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
! Licensed under the MIT license. See LICENSE file in the project directory for details.

!==============================================================
! Fortran interfaces -- uno_fortran.f90
!==============================================================

!---------------------------------------------
! uno_create_model
!---------------------------------------------
function uno_create_model(problem_type, number_variables, variables_lower_bounds, &
                          variables_upper_bounds, base_indexing) &
   result(model)
   character(len=*) :: problem_type
   integer(uno_int), value :: number_variables
   real(c_double) :: variables_lower_bounds(*)
   real(c_double) :: variables_upper_bounds(*)
   integer(uno_int), value :: base_indexing
   type(c_ptr) :: model
   character(c_char), allocatable :: problem_type_c(:)
   integer :: i, n

   interface
      function uno_create_model_c(problem_type, number_variables, &
                                  variables_lower_bounds, variables_upper_bounds, &
                                  base_indexing) &
         result(model) &
         bind(C, name="uno_create_model")
         import :: c_char, uno_int, c_double, c_ptr
         character(c_char) :: problem_type(*)
         integer(uno_int), value :: number_variables
         real(c_double) :: variables_lower_bounds(*)
         real(c_double) :: variables_upper_bounds(*)
         integer(uno_int), value :: base_indexing
         type(c_ptr) :: model
      end function uno_create_model_c
   end interface

   n = len_trim(problem_type)
   if (index(problem_type(1:n), c_null_char) > 0) then
      model = c_null_ptr
      return
   end if
   allocate(problem_type_c(n+1))
   do i = 1, n
      problem_type_c(i) = problem_type(i:i)
   end do
   problem_type_c(n+1) = c_null_char
   model = uno_create_model_c(problem_type_c, number_variables, variables_lower_bounds, &
                              variables_upper_bounds, base_indexing)
   deallocate(problem_type_c)
end function uno_create_model

!---------------------------------------------
! uno_create_unconstrained_model
!---------------------------------------------
function uno_create_unconstrained_model(problem_type, number_variables, base_indexing) &
   result(model)
   character(len=*) :: problem_type
   integer(uno_int), value :: number_variables
   integer(uno_int), value :: base_indexing
   type(c_ptr) :: model
   character(c_char), allocatable :: problem_type_c(:)
   integer :: i, n

   interface
      function uno_create_unconstrained_model_c(problem_type, number_variables, &
                                                base_indexing) &
         result(model) &
         bind(C, name="uno_create_unconstrained_model")
         import :: c_char, uno_int, c_ptr
         character(c_char) :: problem_type(*)
         integer(uno_int), value :: number_variables
         integer(uno_int), value :: base_indexing
         type(c_ptr) :: model
      end function uno_create_unconstrained_model_c
   end interface

   n = len_trim(problem_type)
   if (index(problem_type(1:n), c_null_char) > 0) then
      model = c_null_ptr
      return
   end if
   allocate(problem_type_c(n+1))
   do i = 1, n
      problem_type_c(i) = problem_type(i:i)
   end do
   problem_type_c(n+1) = c_null_char
   model = uno_create_unconstrained_model_c(problem_type_c, number_variables, &
                                            base_indexing)
   deallocate(problem_type_c)
end function uno_create_unconstrained_model

!---------------------------------------------
! uno_set_model_name
!---------------------------------------------
function uno_set_model_name(model, name) result(success)
   type(c_ptr), value :: model
   character(len=*) :: name
   logical(c_bool) :: success
   character(c_char), allocatable :: name_c(:)
   integer :: i, n

   interface
      function uno_set_model_name_c(model, name) &
         result(success) &
         bind(C, name="uno_set_model_name")
         import :: c_ptr, c_char, c_bool
         type(c_ptr), value :: model
         character(c_char) :: name(*)
         logical(c_bool) :: success
      end function uno_set_model_name_c
   end interface

   n = len_trim(name)
   if (index(name(1:n), c_null_char) > 0) then
      success = .false.
      return
   end if
   allocate(name_c(n+1))
   do i = 1, n
      name_c(i) = name(i:i)
   end do
   name_c(n+1) = c_null_char
   success = uno_set_model_name_c(model, name_c)
   deallocate(name_c)
end function uno_set_model_name

!---------------------------------------------
! uno_set_solver_integer_option
!---------------------------------------------
function uno_set_solver_integer_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   integer(uno_int), value :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_set_solver_integer_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_set_solver_integer_option")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         integer(uno_int), value :: option_value
         integer(uno_int) :: status
      end function uno_set_solver_integer_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   status = uno_set_solver_integer_option_c(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_set_solver_integer_option

!---------------------------------------------
! uno_set_solver_double_option
!---------------------------------------------
function uno_set_solver_double_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   real(c_double), value :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_set_solver_double_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_set_solver_double_option")
         import :: c_ptr, c_char, c_double, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         real(c_double), value :: option_value
         integer(uno_int) :: status
      end function uno_set_solver_double_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   status = uno_set_solver_double_option_c(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_set_solver_double_option

!---------------------------------------------
! uno_set_solver_bool_option
!---------------------------------------------
function uno_set_solver_bool_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   logical(c_bool), value :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_set_solver_bool_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_set_solver_bool_option")
         import :: c_ptr, c_char, c_bool, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         logical(c_bool), value :: option_value
         integer(uno_int) :: status
      end function uno_set_solver_bool_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   status = uno_set_solver_bool_option_c(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_set_solver_bool_option

!---------------------------------------------
! uno_set_solver_string_option
!---------------------------------------------
function uno_set_solver_string_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   character(len=*) :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   character(c_char), allocatable :: option_value_c(:)
   integer :: i, n

   interface
      function uno_set_solver_string_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_set_solver_string_option")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         character(c_char) :: option_value(*)
         integer(uno_int) :: status
      end function uno_set_solver_string_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   n = len_trim(option_value)
   if (index(option_value(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_INVALID_VALUE
      return
   end if
   allocate(option_value_c(n+1))
   do i = 1, n
      option_value_c(i) = option_value(i:i)
   end do
   option_value_c(n+1) = c_null_char
   status = uno_set_solver_string_option_c(solver, option_name_c, option_value_c)
   deallocate(option_name_c)
   deallocate(option_value_c)
end function uno_set_solver_string_option

!---------------------------------------------
! uno_get_solver_option_type
!---------------------------------------------
function uno_get_solver_option_type(solver, option_name) result(solver_option_type)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   integer(uno_int) :: solver_option_type
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_get_solver_option_type_c(solver, option_name) &
         result(solver_option_type) &
         bind(C, name="uno_get_solver_option_type")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         integer(uno_int) :: solver_option_type
      end function uno_get_solver_option_type_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      solver_option_type = UNO_OPTION_TYPE_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   solver_option_type = uno_get_solver_option_type_c(solver, option_name_c)
   deallocate(option_name_c)
end function uno_get_solver_option_type

!---------------------------------------------
! uno_options_iterator_name
!---------------------------------------------
function uno_options_iterator_name(iterator) &
   result(options_iterator_name)
   type(c_ptr), value :: iterator
   character(:), allocatable :: options_iterator_name
   type(c_ptr) :: ptr_options_iterator_name_c
   integer :: i, n
   character(c_char), pointer :: options_iterator_name_c(:)

   interface
      function uno_options_iterator_name_c(iterator) &
         result(options_iterator_name) &
         bind(C, name="uno_options_iterator_name")
         import :: c_ptr
         type(c_ptr), value :: iterator
         type(c_ptr) :: options_iterator_name
      end function uno_options_iterator_name_c
   end interface

   ptr_options_iterator_name_c = uno_options_iterator_name_c(iterator)
   if (.not. c_associated(ptr_options_iterator_name_c)) then
      options_iterator_name = ""
      return
   end if
   call c_f_pointer(ptr_options_iterator_name_c, options_iterator_name_c, [UNO_MAX_STRING_LENGTH + 1])
   n = 0
   do while (n < UNO_MAX_STRING_LENGTH)
      if (options_iterator_name_c(n+1) == c_null_char) exit
      n = n + 1
   end do
   allocate(character(len=n)::options_iterator_name)
   do i = 1, n
      options_iterator_name(i:i) = options_iterator_name_c(i)
   end do
end function uno_options_iterator_name

!---------------------------------------------
! uno_load_solver_option_file
!---------------------------------------------
function uno_load_solver_option_file(solver, file_name) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: file_name
   integer(uno_int) :: status
   character(c_char), allocatable :: file_name_c(:)
   integer :: i, n

   interface
      function uno_load_solver_option_file_c(solver, file_name) &
         result(status) &
         bind(C, name="uno_load_solver_option_file")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: file_name(*)
         integer(uno_int) :: status
      end function uno_load_solver_option_file_c
   end interface

   n = len_trim(file_name)
   if (index(file_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_INVALID_FILE
      return
   end if
   allocate(file_name_c(n+1))
   do i = 1, n
      file_name_c(i) = file_name(i:i)
   end do
   file_name_c(n+1) = c_null_char
   status = uno_load_solver_option_file_c(solver, file_name_c)
   deallocate(file_name_c)
end function uno_load_solver_option_file

!---------------------------------------------
! uno_set_solver_preset
!---------------------------------------------
function uno_set_solver_preset(solver, preset_name) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: preset_name
   integer(uno_int) :: status
   character(c_char), allocatable :: preset_name_c(:)
   integer :: i, n

   interface
      function uno_set_solver_preset_c(solver, preset_name) &
         result(status) &
         bind(C, name="uno_set_solver_preset")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: preset_name(*)
         integer(uno_int) :: status
      end function uno_set_solver_preset_c
   end interface

   n = len_trim(preset_name)
   if (index(preset_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_INVALID_VALUE
      return
   end if
   allocate(preset_name_c(n+1))
   do i = 1, n
      preset_name_c(i) = preset_name(i:i)
   end do
   preset_name_c(n+1) = c_null_char
   status = uno_set_solver_preset_c(solver, preset_name_c)
   deallocate(preset_name_c)
end function uno_set_solver_preset

!---------------------------------------------
! uno_get_method_description
!---------------------------------------------
function uno_get_method_description(solver, method_description) result(status)
   type(c_ptr), value :: solver
   character(:), allocatable, intent(inout) :: method_description
   integer(uno_int) :: status
   type(c_ptr) :: ptr_method_description_c
   integer :: i, n
   character(c_char), pointer :: method_description_c(:)

   interface
      function uno_get_method_description_c(solver, method_description) &
         result(status) &
         bind(C, name="uno_get_method_description")
         import :: c_ptr, uno_int
         type(c_ptr), value :: solver
         type(c_ptr) :: method_description
         integer(uno_int) :: status
      end function uno_get_method_description_c
   end interface

   ptr_method_description_c = c_null_ptr
   status = uno_get_method_description_c(solver, ptr_method_description_c)
   ! method_description is assigned only if the C function wrote it
   if (c_associated(ptr_method_description_c)) then
      call c_f_pointer(ptr_method_description_c, method_description_c, [UNO_MAX_STRING_LENGTH + 1])
      n = 0
      do while (n < UNO_MAX_STRING_LENGTH)
         if (method_description_c(n+1) == c_null_char) exit
         n = n + 1
      end do
      if (allocated(method_description)) deallocate(method_description)
      allocate(character(len=n)::method_description)
      do i = 1, n
         method_description(i:i) = method_description_c(i)
      end do
   end if
end function uno_get_method_description

!---------------------------------------------
! uno_get_solver_integer_option
!---------------------------------------------
function uno_get_solver_integer_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   integer(uno_int) :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_get_solver_integer_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_get_solver_integer_option")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         integer(uno_int) :: option_value
         integer(uno_int) :: status
      end function uno_get_solver_integer_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   status = uno_get_solver_integer_option_c(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_get_solver_integer_option

!---------------------------------------------
! uno_get_solver_double_option
!---------------------------------------------
function uno_get_solver_double_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   real(c_double) :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_get_solver_double_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_get_solver_double_option")
         import :: c_ptr, c_char, c_double, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         real(c_double) :: option_value
         integer(uno_int) :: status
      end function uno_get_solver_double_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   status = uno_get_solver_double_option_c(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_get_solver_double_option

!---------------------------------------------
! uno_get_solver_bool_option
!---------------------------------------------
function uno_get_solver_bool_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   logical(c_bool) :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   integer :: i, n

   interface
      function uno_get_solver_bool_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_get_solver_bool_option")
         import :: c_ptr, c_char, c_bool, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         logical(c_bool) :: option_value
         integer(uno_int) :: status
      end function uno_get_solver_bool_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   status = uno_get_solver_bool_option_c(solver, option_name_c, option_value)
   deallocate(option_name_c)
end function uno_get_solver_bool_option

!---------------------------------------------
! uno_get_solver_string_option
!---------------------------------------------
function uno_get_solver_string_option(solver, option_name, option_value) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: option_name
   character(:), allocatable, intent(inout) :: option_value
   integer(uno_int) :: status
   character(c_char), allocatable :: option_name_c(:)
   type(c_ptr) :: ptr_option_value_c
   integer :: i, n
   character(c_char), pointer :: option_value_c(:)

   interface
      function uno_get_solver_string_option_c(solver, option_name, option_value) &
         result(status) &
         bind(C, name="uno_get_solver_string_option")
         import :: c_ptr, c_char, uno_int
         type(c_ptr), value :: solver
         character(c_char) :: option_name(*)
         type(c_ptr) :: option_value
         integer(uno_int) :: status
      end function uno_get_solver_string_option_c
   end interface

   n = len_trim(option_name)
   if (index(option_name(1:n), c_null_char) > 0) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   allocate(option_name_c(n+1))
   do i = 1, n
      option_name_c(i) = option_name(i:i)
   end do
   option_name_c(n+1) = c_null_char
   ptr_option_value_c = c_null_ptr
   status = uno_get_solver_string_option_c(solver, option_name_c, ptr_option_value_c)
   deallocate(option_name_c)
   ! option_value is assigned only if the C function wrote it
   if (c_associated(ptr_option_value_c)) then
      call c_f_pointer(ptr_option_value_c, option_value_c, [UNO_MAX_STRING_LENGTH + 1])
      n = 0
      do while (n < UNO_MAX_STRING_LENGTH)
         if (option_value_c(n+1) == c_null_char) exit
         n = n + 1
      end do
      if (allocated(option_value)) deallocate(option_value)
      allocate(character(len=n)::option_value)
      do i = 1, n
         option_value(i:i) = option_value_c(i)
      end do
   end if
end function uno_get_solver_string_option
