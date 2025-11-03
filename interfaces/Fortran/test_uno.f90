program test_uno
   include 'uno.f90'
   integer(c_int32_t) :: major, minor, patch
   call uno_get_version(major, minor, patch)
   print *, 'Uno version:', major, '.', minor, '.', patch
end program test_uno