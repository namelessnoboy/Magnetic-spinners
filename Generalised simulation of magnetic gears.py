import numpy as np
import matplotlib.pyplot as plt

def distances(X1, Y1, X2, Y2, A1, A2, m1, m2):
   d = ((X1 - X2)**2 + (Y1 - Y2)**2)**(1/2)
   r11 = ((R*(np.cos(A1)-np.cos(A2))-d)**2 + (R*(np.sin(A1)-np.sin(A2)))**2)**(0.5)
   r12 = ((R*(np.cos(A1)-np.cos(A2 + 2*np.pi))-d)**2 + (R*(np.sin(A1)-np.sin(A2 + 2*np.pi)))**2)**(0.5)
   r13 = ((R*(np.cos(A1)-np.cos(A2 + 4*np.pi))-d)**2 + (R*(np.sin(A1)-np.sin(A2 + 4*np.pi)))**2)**(0.5)
   r21 = ((R*(np.cos(A1 + 2*np.pi)-np.cos(A2))-d)**2 + (R*(np.sin(A1 + 2*np.pi)-np.sin(A2)))**2)**(0.5)
   r22 = ((R*(np.cos(A1 + 2*np.pi)-np.cos(A2 + 2*np.pi))-d)**2 + (R*(np.sin(A1 + 2*np.pi)-np.sin(A2 + 2*np.pi)))**2)**(0.5)
   r23 = ((R*(np.cos(A1 + 2*np.pi)-np.cos(A2 + 4*np.pi))-d)**2 + (R*(np.sin(A1 + 2*np.pi)-np.sin(A2 + 4*np.pi)))**2)**(0.5)
   r31 = ((R*(np.cos(A1 + 4*np.pi)-np.cos(A2))-d)**2 + (R*(np.sin(A1 + 4*np.pi)-np.sin(A2)))**2)**(0.5)
   r32 = ((R*(np.cos(A1 + 4*np.pi)-np.cos(A2 + 2*np.pi))-d)**2 + (R*(np.sin(A1 + 4*np.pi)-np.sin(A2 + 2*np.pi)))**2)**(0.5)
   r33 = ((R*(np.cos(A1 + 4*np.pi)-np.cos(A2 + 4*np.pi))-d)**2 + (R*(np.sin(A1 + 4*np.pi)-np.sin(A2 + 4*np.pi)))**2)**(0.5)
   #now we find dr/dA2 / r^4
   np.seterr(divide='ignore', invalid='ignore')
   dr11da = (((R*(np.cos(A1)-np.cos(A2)))-d)*(R*np.sin(A2))-R*(np.sin(A1)-np.sin(A2))*(R*np.cos(A2)))/r11**5
   dr12da = (((R*(np.cos(A1)-np.cos(A2 + 2*np.pi)))-d)*(R*np.sin(A2 + 2*np.pi))-R*(np.sin(A1)-np.sin(A2 + 2*np.pi))*(R*np.cos(A2 + 2*np.pi)))/r12**5
   dr13da = (((R*(np.cos(A1)-np.cos(A2 + 4*np.pi)))-d)*(R*np.sin(A2 + 4*np.pi))-R*(np.sin(A1)-np.sin(A2 + 4*np.pi))*(R*np.cos(A2 + 4*np.pi)))/r13**5
   dr21da = (((R*(np.cos(A1 + 2*np.pi)-np.cos(A2)))-d)*(R*np.sin(A2))-R*(np.sin(A1 + 2*np.pi)-np.sin(A2))*(R*np.cos(A2)))/r21**5
   dr22da = (((R*(np.cos(A1 + 2*np.pi)-np.cos(A2 + 2*np.pi)))-d)*(R*np.sin(A2 + 2*np.pi))-R*(np.sin(A1 + 2*np.pi)-np.sin(A2 + 2*np.pi))*(R*np.cos(A2 + 2*np.pi)))/r22**5
   dr23da = (((R*(np.cos(A1 + 2*np.pi)-np.cos(A2 + 4*np.pi)))-d)*(R*np.sin(A2 + 4*np.pi))-R*(np.sin(A1 + 2*np.pi)-np.sin(A2 + 4*np.pi))*(R*np.cos(A2 + 4*np.pi)))/r23**5
   dr31da = (((R*(np.cos(A1 + 4*np.pi)-np.cos(A2)))-d)*(R*np.sin(A2))-R*(np.sin(A1 + 4*np.pi)-np.sin(A2))*(R*np.cos(A2)))/r11**5
   dr32da = (((R*(np.cos(A1 + 4*np.pi)-np.cos(A2 + 2*np.pi)))-d)*(R*np.sin(A2 + 2*np.pi))-R*(np.sin(A1 + 4*np.pi)-np.sin(A2 + 2*np.pi))*(R*np.cos(A2 + 2*np.pi)))/r11**5
   dr33da = (((R*(np.cos(A1 + 4*np.pi)-np.cos(A2 + 4*np.pi)))-d)*(R*np.sin(A2 + 4*np.pi))-R*(np.sin(A1 + 4*np.pi)-np.sin(A2 + 4*np.pi))*(R*np.cos(A2 + 4*np.pi)))/r11**5
   vars()["dUda" + str(m1) + str(m2)] = -3*k*(dr11da + dr12da + dr13da + dr21da + dr22da + dr23da + dr31da + dr32da + dr33da)
   return vars()["dUda" + str(m1) + str(m2)]
#radius of spinners R
R = 0.033
#moment of inertia of spinners I
I = 66.8*10**-7
#magnetic moment m
m = 0.33
#k constant in potential energy calculation U = k/r^3
k = m**2 * 10**(-7)

#To generalise the whole simulation for a lot of spinners
#a represents angle, a1 is angle of first spinner, a2 is second etc.
n0 = int(input("Number of spinners: "))
n = n0

t_start = 0
t_end = 100
dt = 0.01
n_steps = int(round((t_end-t_start)/dt))
t_arr = np.zeros(n_steps + 1)



t_arr[0] = t_start

for n in range (1, n0 + 1):
   x = 'x' + str(n) 
   y = 'y' + str(n)
   a = 'a' + str(n) + "_arr"
   vars()[x] = float(input(x + "="))
   vars()[y] = float(input(y + "="))
   vars()[a] = np.zeros(n_steps + 1)
   vars()['v' + str(n) + '_arr'] = np.zeros(n_steps + 1)
   vars()['v' + str(n) + '_arr'][0] = float(input('speed='))
   vars()['a' + str(n) + '_arr'][0] = float(input('a' + str(n) + '='))

for i in range (1, n_steps + 1): 
    if n0 == 2:
        dUda12 = distances(x1, x2, y1, y2, a1_arr[i-1], a2_arr[i-1], 1, 2) 

        dUda2 = dUda12
        dUda1 = -dUda12

        A1 = -dUda1/I
        A2 = -dUda2/I

        v1_arr[i] = v1_arr[i-1] + A1*dt
        v2_arr[i] = v2_arr[i-1] + A2*dt

        a1_arr[i] = a1_arr[i-1] + (v1_arr[i-1] + A1*dt)*dt
        a2_arr[i] = a2_arr[i-1] + (v2_arr[i-1] + A2*dt)*dt

        if np.abs(a1_arr[i]) >= 2*np.pi:
            a1_arr[i] = 0
        if np.abs(a2_arr[i]) >= 2*np.pi:
            a2_arr[i] = 0
        
        
    if n0 == 3:
        dUda12 = distances(x1, x2, y1, y2, a1_arr[i-1], a2_arr[i-1], 1, 2)
        dUda23 = distances(x2, x3, y2, y3, a2_arr[i-1], a3_arr[i-1], 2, 3)
        dUda13 = distances(x1, x3, y1, y3, a1_arr[i-1], a3_arr[i-1], 1, 3)

        dUda1 = -dUda12 - dUda13
        dUda2 = dUda12 - dUda23
        dUda3 = dUda13 + dUda23

        A1 = -dUda1/I
        A2 = -dUda2/I
        A3 = -dUda3/I

        v1_arr[i] = v1_arr[i-1] + A1*dt
        v2_arr[i] = v2_arr[i-1] + A2*dt
        v3_arr[i] = v3_arr[i-1] + A3*dt

        a1_arr[i] = a1_arr[i-1] + (v1_arr[i-1] + A1*dt)*dt
        a2_arr[i] = a2_arr[i-1] + (v2_arr[i-1] + A2*dt)*dt
        a3_arr[i] = a3_arr[i-1] + (v3_arr[i-1] + A3*dt)*dt

        if np.abs(a1_arr[i]) >= 2*np.pi:
            a1_arr[i] = 0
        if np.abs(a2_arr[i]) >= 2*np.pi:
            a2_arr[i] = 0
        if np.abs(a3_arr[i]) >= 2*np.pi:
            a3_arr[i] = 0
    if n0 == 4:
        dUda12 = distances(x1, x2, y1, y2, a1_arr[i-1], a2_arr[i-1], 1, 2)
        dUda13 = distances(x1, x3, y1, y3, a1_arr[i-1], a3_arr[i-1], 1, 3)
        dUda14 = distances(x1, x4, y1, y4, a1_arr[i-1], a4_arr[i-1], 1, 4)
        dUda23 = distances(x2, x3, y2, y3, a2_arr[i-1], a3_arr[i-1], 2, 3)
        dUda24 = distances(x2, x4, y2, y4, a2_arr[i-1], a4_arr[i-1], 2, 4)
        dUda34 = distances(x3, x4, y3, y4, a3_arr[i-1], a4_arr[i-1], 3, 4)

        dUda1 = -dUda12 - dUda13 - dUda14
        dUda2 = dUda12 - dUda23 - dUda24
        dUda3 = dUda13 + dUda23 - dUda34
        dUda4 = dUda14 + dUda24 + dUda34

        A1 = -dUda1/I
        A2 = -dUda2/I
        A3 = -dUda3/I
        A4 = -dUda4/I

        v1_arr[i] = v1_arr[i-1] + A1*dt
        v2_arr[i] = v2_arr[i-1] + A2*dt
        v3_arr[i] = v3_arr[i-1] + A3*dt
        v4_arr[i] = v4_arr[i-1] + A4*dt
    
        a1_arr[i] = a1_arr[i-1] + (v1_arr[i-1] + A1*dt)*dt
        a2_arr[i] = a2_arr[i-1] + (v2_arr[i-1] + A2*dt)*dt
        a3_arr[i] = a3_arr[i-1] + (v3_arr[i-1] + A3*dt)*dt
        a4_arr[i] = a4_arr[i-1] + (v4_arr[i-1] + A4*dt)*dt

        if np.abs(a1_arr[i]) >= 2*np.pi:
            a1_arr[i] = 0
        if np.abs(a2_arr[i]) >= 2*np.pi:
            a2_arr[i] = 0
        if np.abs(a3_arr[i]) >= 2*np.pi:
            a3_arr[i] = 0
        if np.abs(a4_arr[i]) >= 2*np.pi:
            a4_arr[i] = 0
    if n0 == 5:
        dUda12 = distances(x1, x2, y1, y2, a1_arr[i-1], a2_arr[i-1], 1, 2)
        dUda13 = distances(x1, x3, y1, y3, a1_arr[i-1], a3_arr[i-1], 1, 3)
        dUda14 = distances(x1, x4, y1, y4, a1_arr[i-1], a4_arr[i-1], 1, 4)
        dUda15 = distances(x1, x5, y1, y5, a1_arr[i-1], a5_arr[i-1], 1, 5)
        dUda23 = distances(x2, x3, y2, y3, a2_arr[i-1], a3_arr[i-1], 2, 3)
        dUda24 = distances(x2, x4, y2, y4, a2_arr[i-1], a4_arr[i-1], 2, 4)
        dUda25 = distances(x2, x5, y2, y5, a2_arr[i-1], a5_arr[i-1], 2, 5)
        dUda34 = distances(x3, x4, y3, y4, a3_arr[i-1], a4_arr[i-1], 3, 4)
        dUda35 = distances(x3, x5, y3, y5, a3_arr[i-1], a5_arr[i-1], 3, 5)
        dUda45 = distances(x4, x5, y4, y5, a4_arr[i-1], a5_arr[i-1], 4, 5)

        dUda1 = -dUda12 - dUda13 - dUda14 - dUda15
        dUda2 = dUda12 - dUda23 - dUda24 - dUda25
        dUda3 = dUda13 + dUda23 - dUda34 - dUda35
        dUda4 = dUda14 + dUda24 + dUda34 - dUda45
        dUda5 = dUda15 + dUda25 + dUda35 + dUda45

        A1 = -dUda1/I
        A2 = -dUda2/I
        A3 = -dUda3/I
        A4 = -dUda4/I
        A5 = -dUda5/I

        v1_arr[i] = v1_arr[i-1] + A1*dt
        v2_arr[i] = v2_arr[i-1] + A2*dt
        v3_arr[i] = v3_arr[i-1] + A3*dt
        v4_arr[i] = v4_arr[i-1] + A4*dt
        v5_arr[i] = v5_arr[i-1] + A5*dt

        a1_arr[i] = a1_arr[i-1] + (v1_arr[i-1] + A1*dt)*dt
        a2_arr[i] = a2_arr[i-1] + (v2_arr[i-1] + A2*dt)*dt
        a3_arr[i] = a3_arr[i-1] + (v3_arr[i-1] + A3*dt)*dt
        a4_arr[i] = a4_arr[i-1] + (v4_arr[i-1] + A4*dt)*dt
        a5_arr[i] = a5_arr[i-1] + (v5_arr[i-1] + A5*dt)*dt

        if np.abs(a1_arr[i]) >= 2*np.pi:
            a1_arr[i] = 0
        if np.abs(a2_arr[i]) >= 2*np.pi:
            a2_arr[i] = 0
        if np.abs(a3_arr[i]) >= 2*np.pi:
            a3_arr[i] = 0
        if np.abs(a4_arr[i]) >= 2*np.pi:
            a4_arr[i] = 0
        if np.abs(a5_arr[i]) >= 2*np.pi:
            a5_arr[i] = 0
    if n0== 6:
        dUda12 = distances(x1, x2, y1, y2, a1_arr[i-1], a2_arr[i-1], 1, 2)
        dUda13 = distances(x1, x3, y1, y3, a1_arr[i-1], a3_arr[i-1], 1, 3)
        dUda14 = distances(x1, x4, y1, y4, a1_arr[i-1], a4_arr[i-1], 1, 4)
        dUda15 = distances(x1, x5, y1, y5, a1_arr[i-1], a5_arr[i-1], 1, 5)
        dUda16 = distances(x1, x6, y1, y6, a1_arr[i-1], a6_arr[i-1], 1, 6)
        dUda23 = distances(x2, x3, y2, y3, a2_arr[i-1], a3_arr[i-1], 2, 3)
        dUda24 = distances(x2, x4, y2, y4, a2_arr[i-1], a4_arr[i-1], 2, 4)
        dUda25 = distances(x2, x5, y2, y5, a2_arr[i-1], a5_arr[i-1], 2, 5)
        dUda26 = distances(x2, x6, y2, y6, a2_arr[i-1], a6_arr[i-1], 2, 6)
        dUda34 = distances(x3, x4, y3, y4, a3_arr[i-1], a4_arr[i-1], 3, 4)
        dUda35 = distances(x3, x5, y3, y5, a3_arr[i-1], a5_arr[i-1], 3, 5)
        dUda36 = distances(x3, x6, y3, y6, a3_arr[i-1], a6_arr[i-1], 3, 6)
        dUda45 = distances(x4, x5, y4, y5, a4_arr[i-1], a5_arr[i-1], 4, 5)
        dUda46 = distances(x4, x6, y4, y6, a4_arr[i-1], a6_arr[i-1], 4, 6)
        dUda56 = distances(x5, x6, y5, y6, a5_arr[i-1], a6_arr[i-1], 5, 6)

        dUda1 = -dUda12 - dUda13 - dUda14 - dUda15 - dUda16
        dUda2 = dUda12 - dUda23 - dUda24 - dUda25 - dUda26
        dUda3 = dUda13 + dUda23 - dUda34 - dUda35 - dUda36
        dUda4 = dUda14 + dUda24 + dUda34 - dUda45 - dUda46
        dUda5 = dUda15 + dUda25 + dUda35 + dUda45 - dUda46
        dUda6 = dUda16 + dUda26 + dUda36 + dUda46 + dUda56

        A1 = -dUda1/I
        A2 = -dUda2/I
        A3 = -dUda3/I
        A4 = -dUda4/I
        A5 = -dUda5/I
        A6 = -dUda6/I

        v1_arr[i] = v1_arr[i-1] + A1*dt
        v2_arr[i] = v2_arr[i-1] + A2*dt
        v3_arr[i] = v3_arr[i-1] + A3*dt
        v4_arr[i] = v4_arr[i-1] + A4*dt
        v5_arr[i] = v5_arr[i-1] + A5*dt
        v5_arr[i] = v5_arr[i-1] + A5*dt
        v6_arr[i] = v6_arr[i-1] + A6*dt

        a1_arr[i] = a1_arr[i-1] + (v1_arr[i-1] + A1*dt)*dt
        a2_arr[i] = a2_arr[i-1] + (v2_arr[i-1] + A2*dt)*dt
        a3_arr[i] = a3_arr[i-1] + (v3_arr[i-1] + A3*dt)*dt
        a4_arr[i] = a4_arr[i-1] + (v4_arr[i-1] + A4*dt)*dt
        a5_arr[i] = a5_arr[i-1] + (v5_arr[i-1] + A5*dt)*dt
        a6_arr[i] = a6_arr[i-1] + (v6_arr[i-1] + A6*dt)*dt

        if np.abs(a1_arr[i]) >= 2*np.pi:
            a1_arr[i] = 0
        if np.abs(a2_arr[i]) >= 2*np.pi:
            a2_arr[i] = 0
        if np.abs(a3_arr[i]) >= 2*np.pi:
            a3_arr[i] = 0
        if np.abs(a4_arr[i]) >= 2*np.pi:
            a4_arr[i] = 0
        if np.abs(a5_arr[i]) >= 2*np.pi:
            a5_arr[i] = 0
        if np.abs(a6_arr[i]) >= 2*np.pi:
            a6_arr[i] = 0
    t_arr[i] = t_arr[i-1] + dt
print(v2_arr[0])
fig = plt.figure()                                  # create figure
plt.plot(t_arr, a1_arr, linewidth = 4, label = 'angle spinner 1')    # plot P to t
plt.plot(t_arr, a2_arr, linewidth = 4, label = 'angle spinner 2') 
plt.plot(t_arr, v1_arr, linewidth = 4, label = 'speed spinner 1') 
plt.plot(t_arr, v2_arr, linewidth = 4, label = 'speed spinner 2')  
plt.plot(t_arr, v3_arr, linewidth = 4, label = 'speed spinner 3')
plt.plot(t_arr, v4_arr, linewidth = 4, label = 'speed spinner 4')
plt.plot(t_arr, v5_arr, linewidth = 4, label = 'speed spinner 5')
plt.plot(t_arr, v6_arr, linewidth = 4, label = 'speed spinner 6')
plt.title('Title', fontsize = 12)    # add some title to your plot
plt.xlabel('t (in seconds)', fontsize = 12)
plt.ylabel('E(t), a(t)', fontsize = 12)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.grid(True)                        # show grid
plt.axis([t_start, t_end, -10, 10])     # show axes measures
plt.legend()
plt.show()
    