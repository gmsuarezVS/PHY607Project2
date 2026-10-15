import numpy as np

from hyperbolic_dynamics import equations, rk4_step


def step_with_wall(state, dt, alpha, length, y0, Tw, Pd, rng):
    # Part II calculates one proposed step
    new_state = rk4_step(state, dt, equations, alpha, length)

    # Part III checks every particle
    for i in range(len(state)):
        if new_state[i, 1] < y0:
            new_state[i, 1] = y0

            if rng.random() >= Pd:
                # Specular reflection
                new_state[i, 3] = abs(new_state[i, 3])

            else:
                # Diffuse reflection
                new_state[i, 2] = y0 * np.sqrt(Tw) * rng.normal()
                new_state[i, 3] = y0 * np.sqrt(
                    -2.0 * Tw * np.log(1.0 - rng.random())
                )

    # Periodic boundary in x
    new_state[:, 0] %= length

    return new_state


def evolve(state, dt, steps, alpha, length, y0, Tw, Pd, rng):
    history = np.empty((steps + 1, len(state), 4))
    history[0] = state

    for step in range(steps):
        state = step_with_wall(
            state, dt, alpha, length, y0, Tw, Pd, rng
        )
        history[step + 1] = state

    return history