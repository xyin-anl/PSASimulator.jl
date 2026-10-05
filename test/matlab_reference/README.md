# MATLAB reference values

`reference.csv` holds the objectives and constraints of the original MATLAB simulator
([PEESEgroup/PSA](https://github.com/PEESEgroup/PSA), commit `3de0832320`) for the 48 demo cases,
computed with GNU Octave 10.3 (`N = 10`). The test suite compares `psacycle` against these values.

To regenerate:

```bash
julia export_cases.jl                      # writes cases.csv
git clone https://github.com/PEESEgroup/PSA && cd PSA
cp /path/to/test/matlab_reference/run_cases.m .
octave-cli run_cases.m /path/to/cases.csv 1 48 out.csv
```

Each case takes 2 to 20 minutes in Octave, so splitting the rows over several processes helps.
The columns of `out.csv` are scenario, material row, objective 1, objective 2, constraints 1 to 3 and run time;
`reference.csv` adds the scenario and material names.
