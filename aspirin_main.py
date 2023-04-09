import numpy as np
from scipy.integrate import odeint


def aspirin_main(BW, dose_step):
    # Constants
    phi = 0.7
    PS = 0.5
    min_effective = 0.03
    max_safe = 0.1

    # Scaling of volumes with body weight
    Vp = 40 * (BW) ** 0.99
    Vg = 2.5 * (BW) ** 0.86
    ke = 2.1 * (BW) ** (0.75)

    # Initial conditions
    dose = dose_step * 50
    Cg0 = dose / Vg
    Cp0 = 0
    y0 = np.array([Cp0, Cg0])

    # Time array
    tmax = 24 * 60
    time = np.linspace(0, tmax, 1000)

    while dose > 0:
        # Solve the ODEs
        sol = odeint(aspirin_odes, y0, time, args=(Vp, ke, phi, PS, Vg))

        # Get the time when plasma concentration drops below MEC after peak
        peak_idx = np.argmax(sol[:, 0])
        idx = np.argmax(sol[peak_idx:, 0] < min_effective) + peak_idx
        time_next_dose = time[idx]

        # Add a second dose
        dose_2 = dose
        Cg02 = sol[idx, 1] + dose_2 / Vg
        Cp02 = sol[idx, 0]
        y02 = np.array([Cp02, Cg02])

        sol2 = odeint(aspirin_odes, y02, time, args=(Vp, ke, phi, PS, Vg))
        t2 = time + time_next_dose

        # If the concentration exceeds max_safe, decrease the dose and try again
        if np.max(sol2[:, 0]) > max_safe:
            dose -= dose_step
            Cg0 = dose / Vg
            y0 = np.array([0, Cg0])
        else:
            break

    # Return the data for plotting along with the dose and time_next_dose
    plot_data = {
        't1': time[:idx].tolist(),
        'y1': sol[:idx, 0].tolist(),
        't2': t2.tolist(),
        'y2': sol2[:, 0].tolist(),
        'min_effective': min_effective,
        'max_safe': max_safe
    }

    return dose, time_next_dose, plot_data


def aspirin_odes(y, t, Vp, ke, phi, PS, Vg):
    dydt = [
        (PS / Vp) * y[1] - (PS * phi / Vp) * y[0] - (ke / Vp) * y[0],  # mol balance in plasma
        -(PS / Vg) * y[1] + (PS * phi / Vg) * y[0]  # mol balance in gut
    ]

    return dydt
