import numpy as np


# state[i] = [x, y, vx, vy]
def coulomb_force(state, alpha, length):
    number = len(state)
    force = np.zeros((number, 2))
    scale = 2.0 * np.pi / length
    factor = alpha * np.pi / length

    for i in range(number):
        xi, yi = state[i, 0], state[i, 1]

        for j in range(number):
            if i == j:
                continue

            xj, yj = state[j, 0], state[j, 1]
            x_phase = scale * (xi - xj)
            y_minus = scale * (yi - yj)
            y_plus = scale * (yi + yj)

            d_minus = np.cosh(y_minus) - np.cos(x_phase)
            d_plus = np.cosh(y_plus) - np.cos(x_phase)

            force[i, 0] += factor * np.sin(x_phase) * (
                1.0 / d_minus - 1.0 / d_plus
            )
            force[i, 1] += factor * (
                np.sinh(y_minus) / d_minus - np.sinh(y_plus) / d_plus
            )

    return force


def equations(state, alpha, length):
    y = state[:, 1]
    vx = state[:, 2]
    vy = state[:, 3]
    force = coulomb_force(state, alpha, length)

    ax = 2.0 * vx * vy / y + y**2 * force[:, 0]
    ay = (vy**2 - vx**2) / y + y**2 * force[:, 1]
    return np.column_stack((vx, vy, ax, ay))


def rk4_step(state, dt, derivative, *args):
    k1 = derivative(state, *args)
    k2 = derivative(state + dt * k1 / 2.0, *args)
    k3 = derivative(state + dt * k2 / 2.0, *args)
    k4 = derivative(state + dt * k3, *args)
    return state + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def evolve(state, dt, steps, alpha, length):
    history = np.empty((steps + 1, len(state), 4))
    history[0] = state

    for step in range(steps):
        state = rk4_step(state, dt, equations, alpha, length)
        state[:, 0] %= length
        history[step + 1] = state

    return history


def wall_crossings(state_before, state_after, y0):