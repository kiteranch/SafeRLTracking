import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os
from scipy.integrate import trapezoid

# -------------------- Load Data --------------------
# mat_file = '01_tracking_eight_rl.mat'
mat_file = '01_tracking_eight_saferl.mat'
# mat_file = 'compare_01_tracking_eight.mat'
mat_path = os.path.join('data', mat_file)   # Modify according to the actual path
data = loadmat(mat_path)
print(f"Load data {mat_file}.")

# -------------------- Environment and Obstacles --------------------
x_low, x_high, y_low, y_high = -1.6, 1.6, -1.6, 2.4
obstacles = [
    (np.array([0.0, 1.0]), 0.25),
    (np.array([-0.75, 0.0]), 0.2),
    (np.array([0.6, -0.4]), 0.2),
]
robot_radius = 0.24

# -------------------- Extract Variables --------------------
print(data['time'].flatten()[-1])
T_plot =53.0                      # Set the duration to be displayed (seconds) 
t = data['time'].flatten()
idx = t <= T_plot
t = t[idx]
state = data['state'][idx, :]                      # (N, 3) -> [x, y, psi]
ref = data['ref'][idx, :]                          # (N, 3) -> [x_d, y_d, psi_d]
# ud = data['ud'][idx, :]                            # (N, 2)
# mu = data['mu'][idx, :]                            # (N, 2)
cmd = data['cmd'][idx, :]                          # (N, 2)

# The actual position and heading angle from the state
x, y, psi = state[:, 0], state[:, 1], state[:, 2]
# Desired Trajectory
x_d, y_d, psi_d = ref[:, 0], ref[:, 1], ref[:, 2]

# Position Error (Euclidean Distance)
pos_err = np.sqrt((x - x_d)**2 + (y - y_d)**2)
# Heading Error
psi_err = np.arctan2(np.sin(psi - psi_d), np.cos(psi - psi_d))

# -------------------- Figure 1: State Trajectory --------------------
fig = plt.figure(1)
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim(-3, 3)
plt.ylim(-1.8, 2.6)
plt.xticks(np.arange(-2, 3, 1))

th = np.linspace(0, 2*np.pi, 100)

# 1) Robot Outline (alpha=0.3)
# α=0.3 Red on a white background → Pre-mixed Colors
col_robot = 0.3*np.array([1, 0, 0]) + 0.7*np.array([1, 1, 1])   # ≈ [1, 0.7, 0.7]

for i in range(0, len(x), 5):
    plt.plot(x[i] + robot_radius*np.cos(th),
            y[i] + robot_radius*np.sin(th),
            color=col_robot, linewidth=0.5)   # no alpha
# 2) Manually draw grid lines (overlapping the robot's outline)
xticks = plt.xticks()[0]
yticks = plt.yticks()[0]
xlims = plt.xlim(); ylims = plt.ylim()
dis = 0.00
for xt in xticks:
    plt.plot([xt, xt], [ylims[0]+dis, ylims[1]-dis],
            color=(0.6902, 0.6902, 0.6902), linewidth=0.8)
for yt in yticks:
    plt.plot([xlims[0]+dis, xlims[1]-dis], [yt, yt],
            color=(0.6902, 0.6902, 0.6902), linewidth=0.8)

# 3) Boundaries and Obstacles
for xy in [([x_low, x_low], [y_low, y_high]),
           ([x_high, x_high], [y_low, y_high]),
           ([x_low, x_high], [y_low, y_low]),
           ([x_low, x_high], [y_high, y_high])]:
    plt.plot(xy[0], xy[1], 'y', linewidth=2)

for c, r in obstacles:
    plt.plot(c[0] + r*np.cos(th), c[1] + r*np.sin(th),
            'k-.', linewidth=1.3)

# 4) Reference and Actual Trajectory
xp_plot,  = plt.plot(x,   y,   'r-',  linewidth=1.5, label='Actual')
xpd_plot, = plt.plot(x_d, y_d, 'b--', linewidth=1.5, label='Reference')

# 5) Add: Mark actual points at specific times
mark_times = [3.18, 18.63, 31.18, 44.00]
colors  = ['#4DBEEE', '#77AC30', '#7E2F8E', '#EDB120']
offsets = [(0.15, 0.0), (-0.75, -0.2), (-0.95, -0.2), (-0.5, -0.2)]
for i, t_mark in enumerate(mark_times):
    mi = np.argmin(np.abs(t - t_mark))
    plt.plot(x[mi], y[mi], 'x', markersize=8, markeredgewidth=1.5,
            color=colors[i])
    plt.text(x[mi] + offsets[i][0], y[mi] + offsets[i][1],
            f'$t={t_mark:.2f}$ s', color='k')

plt.xlabel(r'$X$ (m)')
plt.ylabel(r'$Y$ (m)')
# plt.axis('equal')
# plt.grid(True)
plt.legend([xp_plot, xpd_plot], ['Actual', 'Reference'],
          loc='upper right')
# plt.tight_layout()

# -------------------- Figure 2: Tracking Error, Control Input --------------------
plt.figure(2, figsize=[12,10])

plt.subplot(2, 2, 1)
plt.plot(t, pos_err, 'g-')
plt.xlabel('Time (s)')
plt.ylabel('Position error (m)')
plt.title('Position')
plt.grid(True)

plt.subplot(2, 2, 2)
# plt.plot(t, psi, 'r-', label='actual heading')
# plt.plot(t, psi_d, 'b-', label='desired heading')
plt.plot(t, psi_err, 'm-')
plt.xlabel('Time (s)')
plt.ylabel('Heading error (rad)')
plt.title('Heading')
plt.grid(True)

# plt.subplot(2, 2, 3)
# plt.plot(t, cmd[:, 0], label=r'$u_v$')
# plt.plot(t, cmd[:, 1], label=r'$u_\omega$', linestyle='-')
# plt.xlabel('Time (s)')
# plt.ylabel('Command')
# plt.title('Control Input')
# plt.legend()
# plt.grid(True)

# plt.subplot(2, 2, 4)
# plt.plot(t, ud, '-', label=r'$u_d$')
# plt.plot(t, mu, '--', label=r'$u_e$')
# plt.xlabel('Time (s)')
# plt.ylabel('ud & mu')
# plt.title('Tracking Controller')
# plt.legend()
# plt.grid(True)

# plt.tight_layout()

# -------------------- Figure 3: Control Barrier Function and safety constraints --------------------
if 'B' in data:
    u_safe = data['u_safe'][idx, :]                 # (N, 2)
    B = data['B'].flatten()[idx]                    # (N,)
    h_arr = data['h'][idx, :]      # shape (N, num_constraints)
    H_arr = data['H'][idx, :]      # shape (N, num_constraints)

    plt.figure(3, figsize=(12,10))
    plt.subplot(2, 2, 1)
    # Define a time interval
    # t_low, t_high = 50.0, 60.0
    # idx_interval = np.where((t >= t_low) & (t <= t_high))[0]  # Get the indices within the range
    # if len(idx_interval) > 0:
    #     local_max_idx = np.argmax(B[idx_interval]) # Get the index of the local maximum within the interval
    #     global_idx = idx_interval[local_max_idx]  # Get the global index
    #     t_max = t[global_idx]                 # Get the time of the maximum value
    #     print(f"In the time interval [{t_low}, {t_high}], the maximum value of B is {B[global_idx]:.4f} at t = {t_max:.4f} s")
    # else:
    #     print(f"No data points within the time interval [{t_low}, {t_high}]")
    print('The minimum value of B:',min(B))
    plt.plot(t, B, 'r-')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$B(x,v)$')
    plt.title('Lyapunov-like CBF')
    plt.grid()

    num_constraints = h_arr.shape[1]

    # Compute the minimum value of each constraint
    h_min = np.min(h_arr, axis=0)
    H_min = np.min(H_arr, axis=0)
    print("The minimum value of each constraint h_i:")
    for i in range(num_constraints):
        print(f"  h_{i+1}: {h_min[i]:.6f}")
    print("The minimum value of each constraint H_i:")
    for i in range(num_constraints):
        print(f"  H_{i+1}: {H_min[i]:.6f}")

    plt.subplot(2, 2, 2)
    plt.plot(t, u_safe[:, 0], label=r'$u_{s1}$')
    plt.plot(t, u_safe[:, 1], label=r'$u_{s2}$')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$u_s(x,v)$')
    plt.title('CBF-based Safeguarding Controller')
    plt.legend()
    plt.grid()
    
    plt.subplot(2, 2, 3)
    for i in range(num_constraints):
        plt.plot(t, h_arr[:, i], label=f'$h_{i+1}$')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$h_i(x)$')
    plt.title('Safety constraints h')
    plt.grid(True)
    plt.legend(loc='upper right', ncol=2)

    plt.subplot(2, 2, 4)
    for i in range(num_constraints):
        plt.plot(t, H_arr[:, i], label=f'$H_{i+1}$')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$H_i(x)$')
    plt.title('Safety constraints H')
    plt.grid(True)
    plt.legend(loc='upper right', ncol=2)
    # plt.tight_layout()
else:
    print(f"No u_safe/B/h/H data found, safeguarding may be disabled.")

# -------------------- Figure 4: Learning Weights and cost --------------------
if 'Wc' in data:
    Wc = data['Wc'][idx, :]          # (N, L)
    Wa = data['Wa'][idx, :]          # (N, L)
    cbar = data['cbar'].flatten()[idx]   # (N,)
    cost = data['cum_cost'].flatten()[idx]   # (N,)
    print('cost:',cost[-1])
    L = Wc.shape[1]

    plt.figure(4, figsize=(12,10))
    plt.subplot(2, 2, 1)
    plt.plot(t, Wc, linewidth=1.5)
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\hat{w}_c(t)$')
    plt.grid(True)
    plt.legend([f'$w_{{c{i+1}}}$' for i in range(L)], ncol=2)

    plt.subplot(2, 2, 2)
    plt.plot(t, Wa, linewidth=1.5)
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\hat{w}_a(t)$')
    plt.grid(True)
    plt.legend([f'$w_{{a{i+1}}}$' for i in range(L)], ncol=2)

    plt.subplot(2, 2, 3)
    plt.plot(t, cbar, linewidth=1.5, color='green')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$c_{bar}(t)$')
    plt.grid(True)

    plt.subplot(2, 2, 4)
    plt.plot(t, cost, label='Integral Cost')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\int_0^t (e^T Q e + \mu^T R \mu) d\tau$')
    plt.grid(True)
    plt.legend()
    # plt.tight_layout()

# ---------- Quantitative Metrics ----------
T = t[-1] - t[0]
error = np.sqrt(pos_err**2 + psi_err**2)
RMSE = np.sqrt(trapezoid(error**2, t) / T)

u_abs = abs(cmd[:, 0]) + abs(cmd[:, 1])
IAU = trapezoid(y=u_abs, x=t)  # ∫ |u| dt

print(f"RMSE: {RMSE:.4f}")
print(f"IAU: {IAU:.4f}")

# -------------------- Save as EPS --------------------
# fig.savefig('exp01_state.eps', bbox_inches='tight', pad_inches=0)
plt.show()

