# Copyright (c) 2025 Charlie Vanaret
# Licensed under the MIT license. See LICENSE file in the project directory for details.

import unopy

if __name__ == '__main__':
	options = unopy.Options.get_default()
	unopy.Options.set_preset(options, "filtersqp")
	print(options)