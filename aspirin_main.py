import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


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

    # Plotting
    # plt.plot(time[:idx + 1], sol[:idx + 1, 0], 'black', linewidth=1.5, label='First Dose')
    # plt.plot(t2, sol2[:, 0], 'blue', linewidth=1.5, label='Subsequent Doses')
    # plt.axhline(min_effective, color='green', linewidth=1.5, label='Minimum Effective Concentration')
    # plt.axhline(max_safe, color='red', linewidth=1.5, label='Maximum Safe Concentration')
    # plt.xlabel('time (min)')
    # plt.ylabel('Plasma Concentration (mg/ml)')
    # plt.legend()
    # plt.show()
    #
    # print("Dose in mg is:", dose)
    # print("Time between doses in hours is:", time_next_dose / 60)
    print('Dose in mg is:')
    print(dose)
    print('Time between doses in hours is: ')
    print(time_next_dose / 60)

    return dose, time_next_dose


def aspirin_odes(y, t, Vp, ke, phi, PS, Vg):
    dydt = [
        (PS / Vp) * y[1] - (PS * phi / Vp) * y[0] - (ke / Vp) * y[0],  # mol balance in plasma
        -(PS / Vg) * y[1] + (PS * phi / Vg) * y[0]  # mol balance in gut
    ]

    return dydt
