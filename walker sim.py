import matplotlib.pyplot as plt
import numpy as np

#defining certain fixed quantities such as limits
time = 10
g = 9.81
a_slope = 0.5
dt = 0.01
n_steps = int(time/dt)
t_arr = np.zeros(n_steps + 1)
t_arr[0] = 0
I_a = 1
I_b = 1
I_ab = 2 #for phase 1
m_a = 1
m_b = 1
m_ab = m_a + m_b
d_step = 3
d_arr = np.zeros(n_steps + 1)
d_arr[0] = 0
a_1 = 1
#phase 1
a_phase1_change_arr = np.zeros(n_steps + 1)
a_phase1_change_arr[0] = 0
a_phase1 = 1 #how much the thing rotates in phase 1
va1_arr = np.zeros(n_steps + 1)
va1_arr[0] = 0
L1 = 1 #perpendicular distance between total c.g.and pivot on ground
#phase 2
a_phase2_change_arr = np.zeros(n_steps + 1)
a_phase2_change_arr[0] = 0
a_phase2 = 1
va2_arr = np.zeros(n_steps + 1)
va2_arr[0] = 0
a_2 = 2
L2 = 1 #between joint of a and b and c.g. of a
#phase 3
a_phase3 = 1
a_phase3_change_arr = np.zeros(n_steps + 1)
a_phase3_change_arr[0] = 0
va3_arr = np.zeros(n_steps + 1)
va3_arr[0] = 0
a_3 = 2
L3 = 1 #between c.g. of b and joint of a and b


print (n_steps)
n = 0
phase1_init = 1
phase2_init = 0
for i in range(1, n_steps + 1):
    if a_phase1_change_arr[i-1] < a_phase1:
        A = (np.cos(a_1 - a_phase1_change_arr[i-1] - a_slope))*(9.81)*(m_ab)*L1/(I_ab)
        va1_arr[i] = va1_arr[i-1] + A*dt
        a_phase1_change_arr[i] = a_phase1_change_arr[i-1] + va1_arr[i]*dt
        va2_arr[i] = va1_arr[i]
    elif a_phase2_change_arr[i-1] < a_phase2:
        A = np.sin(a_slope + np.pi/2 - a_2 - a_phase2_change_arr[i-1])*g*m_a*L2/I_a
        va2_arr[i] = va2_arr[i-1] + A*dt
        a_phase2_change_arr[i] = a_phase2_change_arr[i-1] + va2_arr[i]*dt
        a_phase1_change_arr[i] = a_phase1_change_arr[i-1]
    elif a_phase3_change_arr[i-1] < a_phase3:
        A = np.cos(a_3 - a_slope - a_phase3_change_arr[i-1])*L3*g*m_b/I_b
        va3_arr[i] = va3_arr[i-1] + A*dt
        a_phase3_change_arr[i] = a_phase3_change_arr[i-1] + va3_arr[i]*dt
        a_phase1_change_arr[i] = a_phase1_change_arr[i-1]
        a_phase2_change_arr[i] = a_phase2_change_arr[i-1]
    else:
        a_phase1_change_arr[i] = 0
        a_phase2_change_arr[i] = 0
        a_phase3_change_arr[i] = 0
        va1_arr[i] = 0
        va2_arr[i] = 0
        va3_arr[i] = 0
        n += 1
    d_arr[i] = n*d_step
    print(a_phase1_change_arr[i], a_phase2_change_arr[i])


    t_arr[i] = t_arr[i-1] + dt
    

fig = plt.figure()                                  # create figure
plt.plot(t_arr, d_arr, linewidth = 4, label = 'distance')    # plot P to t
plt.title('Title', fontsize = 12)    # add some title to your plot
plt.xlabel('t (in seconds)', fontsize = 12)
plt.ylabel('d(t)', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)                        # show grid
plt.axis([0, time, -10, 10])     # show axes measures
plt.legend()
plt.show()