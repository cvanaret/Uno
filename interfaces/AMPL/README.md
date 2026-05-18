## Uno's AMPL interface

To solve an AMPL model in the [.nl format](https://en.wikipedia.org/wiki/Nl_(format)), move to the `build` directory and:
- run `cmake` with the path to the ASL library: `-DAMPLSOLVER=path`
- compile the executable ```make uno_ampl```
- run the command ```./uno_ampl model.nl [-AMPL] [option=value ...]``` where ```[option=value ...]``` is a list of options separated by spaces. If the `-AMPL` flag is supplied, the solution is written to the AMPL solution file `model.sol`

Options can be set in three different ways (with decreasing precedence):
- setting a preset that mimics an existing solver (`preset=[filtersqp|ipopt]`).
- passing an option file (`option_file=file`) that contains an `option value` pair per line. Anything that follows a `#` is treated as a comment and is ignored. The attached `uno_default.opt` contains the default Uno options.
- setting individual options (see the [default options](https://github.com/cvanaret/Uno/blob/main/docs/options.md)).
For an overview of the available strategies, type: ```./uno_ampl --strategies```:

A couple of CUTEst instances are available in the `/examples` directory.
