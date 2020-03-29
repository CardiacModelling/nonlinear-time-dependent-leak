# Nonlinear time-dependent leak current in automated patch-clamp platforms
A study of the nonlinear and time-dependent leak current observed in automated patch-clamp platforms.
This repo contains all data and code for reproducing the results in the paper "*A nonlinear and time-dependent leak current in the presence of calcium fluoride patch clamp seal enhancer*".

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
- [data](./data): Contains all experimental data and protocols.
- [fig](./fig): All result figures.
- [out](./out): Contains all fitting results.

## Acknowledging this work
[PLACEHOLDER]

