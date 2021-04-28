'''
How to fit an IDEA model and compute intervention impact.
'''

# Imports.
import numpy as np
from scipy.optimize import curve_fit

# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0083622
def idea(t, R0, d):
    I = (R0 / ((1 + d) ** t)) ** t
    return I

# Some random data. Replace with real data here.
t = np.arange(1,101)  # Number of serial intervals.
y = np.random.randint(1,10,100)  # Case counts at each time step.

# Intervention time.
split = 80

# Fit the model.
fit_params = curve_fit(idea, t[:split], y[:split], bounds=([0, 0], [np.inf, 1]))
R0 = fit_params[0][0]
d = fit_params[0][1]

# Impact.
impact = np.sum(y[split:] - idea(t[split:], R0, d))

# Print results.
print('R0, d, impact')
print(R0, d, impact)
