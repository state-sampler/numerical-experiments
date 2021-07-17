# Counting knowledge states

The code in this repository estimates the size of an ordinal knowledge space using a variation of the subset simulation algorithm introduced by Au and Beck (2001).

>Au, S.-K. and Beck, J. L. (2001). Estimation of small failure probabilities in high dimensions by subset simulation. Probabilistic Engineering Mechanics, 16(4):263 277.

# Requirements

A Python installation with NumPy is required.  The code has been tested on NumPy 1.19.1 and Python 3.7.9.

# Usage

Before running the Python code, the Fortran source in sampler.f95 must first be built into a Python module using NumPy's [`F2PY`](https://numpy.org/doc/stable/f2py/) tool.

```sh
python -m numpy.f2py -c -m sampler sampler.f95
```

The estimator can then be run in an interactive Python session as follows.

```python
import counting_states
partial_order = np.loadtxt('full_partial_order.csv', delimiter=',')
num_states, Ei_dict = counting_states.run_estimator(partial_order)
```

The file `full_partial_order.csv` contains the partial order for the full knowledge space of 300 items.  To run the estimator on the projected knowledge space of 100 items, use the file `projected_partial_order.csv` instead.
