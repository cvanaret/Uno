# Copyright (c) 2022 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

# https://www.auctoris.co.uk/2017/04/29/calling-c-classes-from-python-with-ctypes/

import ctypes
libuno_C_wrapper = ctypes.cdll.LoadLibrary('./libuno_C_wrapper.so')

class SparseVector(object):
   # constructor
   def __init__(self, capacity):
      # list parameter and result types
      libuno_C_wrapper.SparseVector_new.argtypes = [ctypes.c_uint]
      libuno_C_wrapper.SparseVector_new.restype = ctypes.c_void_p # void*
      libuno_C_wrapper.SparseVector_delete.argtypes = [ctypes.c_void_p]
      libuno_C_wrapper.SparseVector_delete.restype = None # void
      libuno_C_wrapper.SparseVector_insert.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_double]
      libuno_C_wrapper.SparseVector_insert.restype = None # void
      libuno_C_wrapper.SparseVector_clear.argtypes = [ctypes.c_void_p]
      libuno_C_wrapper.SparseVector_clear.restype = None # void
      libuno_C_wrapper.SparseVector_display.argtypes = [ctypes.c_void_p]
      libuno_C_wrapper.SparseVector_display.restype = None # void
      
      self.obj = libuno_C_wrapper.SparseVector_new(capacity)
   
   def __del__(self):
      libuno_C_wrapper.SparseVector_delete(self.obj)
   
   def insert(self, key, value):
      libuno_C_wrapper.SparseVector_insert(self.obj, key, value)

   def clear(self):
      libuno_C_wrapper.SparseVector_clear(self.obj)

   def display(self):
      libuno_C_wrapper.SparseVector_display(self.obj)

