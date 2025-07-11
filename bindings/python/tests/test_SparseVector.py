# Copyright (c) 2022-2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

from unopy import SparseVector

if __name__ == '__main__':
   v = SparseVector(10)
   v.insert(0, 1.5)
   v.insert(1, 4.)
   print(v)
   v.clear()
   v.insert(5, 10.)
   print(v)