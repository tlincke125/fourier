from scipy.io import loadmat
import matplotlib.pyplot as plt
from scipy import fftpack, stats
from scipy.signal import find_peaks
import numpy as np
import os

def reject_outliers(data, m=0.1):
    return abs(data - np.mean(data)) < m * np.std(data)

def peaks(data_file, n, h, r):
    x = loadmat(data_file)

    z_arr = np.array([])
    peak_arr = np.array([])
    
    # nth array is the nth soliton
    # nth element of the nth array is a pair of points for the position at time t
    soliton_positions = [[[], [], []] for i in range(n)]
    max_counter = 1
    counter = 0
    for i in range(len(x['Amat'])):
        peaks, _ = find_peaks(x['Amat'][i], height = h, prominence = 1)
        t = x['t_vec'][0][i]
       
        if len(peaks) > 0 and peaks[-1] > (len(x['z_vec'][0]) - 200):
            break

        if len(peaks) > n:
            peaks = peaks[-n:]

            # Distribute peaks into soliton position array
            for j in peaks:

                z = x['z_vec'][0][j]
                a = x['Amat'][i][j]
                k = n - 1
                while(k >= 0 and len(soliton_positions[k][1]) > 0 and z < soliton_positions[k][1][-1]):
                    k -= 1

                soliton_positions[k][0].append(t)
                soliton_positions[k][1].append(z)
                soliton_positions[k][2].append(a)
                #  N = len(soliton_positions[k][0])
                #  soliton_positions[k][2][0] = (soliton_positions[k][2][0] * (N - 1) + a) / N



    s = np.array([])
    a = np.array([])
    sp = np.array([])

    a0 = x["A0"][0][0]
    u0 = x["U0"][0][0]

    for i in range(n):
        if(len(soliton_positions[i][0]) < 4):
            continue
        soliton_positions[i][0] = np.array(soliton_positions[i][0][4:])
        soliton_positions[i][1] = np.array(soliton_positions[i][1][4:])
        soliton_positions[i][2] = np.array(soliton_positions[i][2][4:])

        slope, intercept, r_value, p_value, std_err = stats.linregress(soliton_positions[i][0], soliton_positions[i][1])

        if(r_value > r):
            s = np.append(s, slope)
            a = np.append(a, np.mean(soliton_positions[i][2]))
            sp = np.append(sp, u0 * theoretical_speed(np.mean(soliton_positions[i][2])))


    return s, a, sp

def theoretical_speed(a0):
    return (2 * a0**2 * np.log(a0) - a0**2 + 1) / (a0 - 1)**2



def display(data_file, n, h):
    x = loadmat(data_file)
    print(len(x['z_vec'][0]))

    for i in range(len(x['Amat'])):
        peaks, _ = find_peaks(x['Amat'][i], height = h, prominence = .5)

        if len(peaks) > 0 and peaks[-1] > (len(x['z_vec'][0]) - 200):
            break

        plt.plot(x['z_vec'][0], x['Amat'][i])
        plt.plot(x['z_vec'][0][peaks], x['Amat'][i][peaks], "vb")
        plt.pause(0.005)
        plt.clf()

#  display('./expData33.mat', 5, 0)

#  4.122208279577089 0.9996846168798904
#  4.445619825701175 0.9996624612682581
#  4.739776501941086 0.9998139109938673
#  5.018645732706719 0.9999047166422371
#  5.2722581905183175 0.9998685179897363



s = np.array([])
sp = np.array([])
a = np.array([])
for i in os.listdir('./'):

    if i.endswith('.mat'):
        print(i)
        s_, a_, sp_ = peaks(i, 5, 0, 0.99)
        #  plt.show()
        #  plt.clf()
        s = np.append(s, s_)
        a = np.append(a, a_)
        sp = np.append(sp, sp_)

plt.clf()
plt.scatter(a, s, s=1, label="Measured Speed")
plt.scatter(a, sp, s=1, label="Theoretical speed")
plt.ylabel("Speed")
plt.xlabel("Amplitude")
plt.title("Speed Amplitude Relationship")
plt.legend()
plt.show()



#  x = loadmat('./expData07.mat')
#
#
#  time_step = 0.0000000001
#
#  s = np.zeros(len(x['Amat']))
#  # first pass
#  for i in range(len(x['Amat'])):
#
#      speed = 0
#      if i > 0:
#          # get the maximum value amplitude
#          j = x['Amat'][i].argmax()
#
#          # Finite differences
#          speed = (x['z_vec'][0][j] - x['z_vec'][0][j - 1]) / (x['t_vec'][0][i] - x['t_vec'][0][i - 1])
#
#          # Very simple filter
#          if speed < 0:
#              speed = 0
#
#      s[i] = speed
#
#
#      plt.plot(x['z_vec'][0], x['Amat'][i])
#      plt.pause(0.005)
#      plt.clf()
#
#
#  #  for i in range(len(x['Amat'])):
#  #      speed = 0
#  #      if i > 0:
#  #          # get the maximum value amplitude
#  #          j = x['Amat'][i].argmax()
#  #
#  #          # Finite differences
#  #          speed = (x['z_vec'][0][j] - x['z_vec'][0][j - 1]) / (x['t_vec'][0][i] - x['t_vec'][0][i - 1])
#  #
#  #          # Very simple filter
#  #          if speed < 0:
#  #              speed = 0
#  #
#  #      s[i] = speed
#  #
#
#
#  plt.show()
