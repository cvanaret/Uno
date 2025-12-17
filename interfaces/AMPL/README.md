## Uno's AMPL interface

To solve an AMPL model in the [.nl format](https://en.wikipedia.org/wiki/Nl_(format)), move to the `build` directory and:
- run `cmake` with the path to the ASL library: `-DAMPLSOLVER=path`
- compile the executable ```make uno_ampl```
- run the command ```./uno_ampl model.nl [-AMPL] [option=value ...]``` where ```[option=value ...]``` is a list of options separated by spaces. If the `-AMPL` flag is supplied, the solution is written to the AMPL solution file `model.sol`

Options can be set in three different ways (with decreasing precedence):
- passing an option file (`option_file=file`) that contains `option value` on each line;
- setting a preset that mimics an existing solver (`preset=[filtersqp|ipopt]`);
- setting individual options (see the [default options](https://github.com/cvanaret/Uno/blob/main/docs/options.md)).
For an overview of the available strategies, type: ```./uno_ampl --strategies```:

A couple of CUTEst instances are available in the `/examples` directory.