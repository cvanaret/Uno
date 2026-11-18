! Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
! Licensed under the MIT license. See LICENSE file in the project directory for details.

!==============================================================
! Fortran interfaces -- uno_fortran.f90
!==============================================================

!---------------------------------------------
! uno_to_c_string
!---------------------------------------------
! converts a Fortran string (without its trailing blanks) into a null-terminated C string. Returns .false. (c_string is
! then not allocated) if the string has more than UNO_MAX_STRING_LENGTH characters or contains a null character, which
! would silently truncate it on the C side
function uno_to_c_string(f_string, c_string) result(valid)
   character(len=*), intent(in) :: f_string
   character(c_char), allocatable, intent(out) :: c_string(:)
   logical :: valid
   integer :: i, n

   valid = .false.
   n = len_trim(f_string)
   ! checked first, so that a very long string is neither scanned nor copied
   if (n > UNO_MAX_STRING_LENGTH) return
   if (index(f_string(1:n), c_null_char) /= 0) return
   allocate(c_string(n+1))
   do i = 1, n
      c_string(i) = f_string(i:i)
   end do
   c_string(n+1) = c_null_char
   valid = .true.
end function uno_to_c_string

!---------------------------------------------
! uno_from_c_string
!---------------------------------------------
! converts a null-terminated C string returned by Uno into a Fortran string. Uno guarantees that the string has at most
! UNO_MAX_STRING_LENGTH characters: at most UNO_MAX_STRING_LENGTH + 1 characters are read, so that the loop terminates
! even if the pointer is wrong
subroutine uno_from_c_string(c_pointer, f_string)
   type(c_ptr), intent(in) :: c_pointer
   character(:), allocatable, intent(inout) :: f_string
   character(c_char), pointer :: c_string(:)
   integer :: i, n

   call c_f_pointer(c_pointer, c_string, [UNO_MAX_STRING_LENGTH + 1])
   n = 0
   do while (n < UNO_MAX_STRING_LENGTH)
      if (c_string(n+1) == c_null_char) exit
      n = n + 1
   end do
   if (allocated(f_string)) deallocate(f_string)
   allocate(character(len=n)::f_string)
   do i = 1, n
      f_string(i:i) = c_string(i)
   end do
end subroutine uno_from_c_string

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

   if (.not. uno_to_c_string(problem_type, problem_type_c)) then
      model = c_null_ptr
      return
   end if
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

   if (.not. uno_to_c_string(problem_type, problem_type_c)) then
      model = c_null_ptr
      return
   end if
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

   if (.not. uno_to_c_string(name, name_c)) then
      success = .false.
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   if (.not. uno_to_c_string(option_value, option_value_c)) then
      status = UNO_OPTIONS_INVALID_VALUE
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      solver_option_type = UNO_OPTION_TYPE_NOT_FOUND
      return
   end if
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
   call uno_from_c_string(ptr_options_iterator_name_c, options_iterator_name)
end function uno_options_iterator_name

!---------------------------------------------
! uno_load_solver_option_file
!---------------------------------------------
function uno_load_solver_option_file(solver, file_name) result(status)
   type(c_ptr), value :: solver
   character(len=*) :: file_name
   integer(uno_int) :: status
   character(c_char), allocatable :: file_name_c(:)

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

   if (.not. uno_to_c_string(file_name, file_name_c)) then
      status = UNO_OPTIONS_INVALID_FILE
      return
   end if
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

   if (.not. uno_to_c_string(preset_name, preset_name_c)) then
      status = UNO_OPTIONS_INVALID_VALUE
      return
   end if
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
      call uno_from_c_string(ptr_method_description_c, method_description)
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
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

   if (.not. uno_to_c_string(option_name, option_name_c)) then
      status = UNO_OPTIONS_NOT_FOUND
      return
   end if
   ptr_option_value_c = c_null_ptr
   status = uno_get_solver_string_option_c(solver, option_name_c, ptr_option_value_c)
   deallocate(option_name_c)
   ! option_value is assigned only if the C function wrote it
   if (c_associated(ptr_option_value_c)) then
      call uno_from_c_string(ptr_option_value_c, option_value)
   end if
end function uno_get_solver_string_option
