# Nonlinear time-dependent leak current in automated patch-clamp platforms
A study of the nonlinear and time-dependent leak current observed in automated patch-clamp platforms.
This repo contains all data and code for reproducing the results in the paper "*[A nonlinear and time-dependent leak current in the presence of calcium fluoride patch clamp seal enhancer](https://doi.org/10.12688/wellcomeopenres.15968.1)*" by Chon Lok Lei, Alan Fabbri, Dominic Whittaker, Michael Clerx, Monique Windley, Adam Hill, Gary Mirams, and Teun de Boer.

## Installing
This repository requires Python 3 with Numpy, Scipy, and Matplotlib.
An example installation:
```bash
python3 -m venv env                # Set up a virtual environment
source env/bin/activate            # Activate the environment
pip install --upgrade pip          # Update pip
pip install -r requirements.txt
```

## Results
- [Figure 2](fig/patch-auto.pdf): Run `patch-auto.py`.
- [Figure 3](fig/patch-manual.pdf): Run `patch-manual.py`.

## Folders
- [method](./method): Contains all modules/utility functions.
- [data](./data): Contains all experimental data and protocols, the individual data files are described in a readme in this folder and there is an associated licence.
- [fig](./fig): All result figures.
- [out](./out): Contains all fitting results.

## License
All the data provided in this repository within the [data](./data) directory are distributed under a [Creative Commons Attribution 4.0 International License](./data/LICENSE.txt).
The rest of the code in this repository are under a [BSD 3-Clause License](./LICENSE.txt).

## Acknowledging this work
Please cite this preprint:  
Lei CL, Fabbri A, Whittaker DG et al. A nonlinear and time-dependent leak current in the presence of calcium fluoride patch-clamp seal enhancer [version 1; peer review: awaiting peer review]. Wellcome Open Res 2020, 5:152 (https://doi.org/10.12688/wellcomeopenres.15968.1) 
