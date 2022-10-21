# Copyright (c) 2022 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

from UnoSparseVector import SparseVector

if __name__ == '__main__':
   v = SparseVector(10)
   v.insert(0, 1.5)
   v.insert(1, 4.)
   print "Vector:"
   v.display()
   v.clear()
   v.insert(5, 10.)
   print "Vector:"
   v.display()
