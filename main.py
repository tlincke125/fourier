from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt


# Gets the speed of the leading edge of the wave
def speed_vec(f):
    x = loadmat(f)

    # Get the speed of the wave
    a = x['A0'][0][0]
    u = x['U0'][0][0]
    c = (2 * a**2 * np.log(a) - a**2 + 1) / (a - 1)**2

    print(c)
    print(u)

    s = np.zeros(len(x['Amat']))
    for i in range(len(x['Amat'])):
        speed = 0.00001
        if i > 0:
            j = x['Amat'][i].argmax()
            speed = (x['z_vec'][0][j] - x['z_vec'][0][j - 1]) / (x['t_vec'][0][i] - x['t_vec'][0][i - 1])
            if speed <= 0:
                speed = 0.00001
        s[i] = speed
    return s


##
# @param u SCALAR
# @param a VECTOR (227 in length)
#
# @return 
def func(u, c):
    B = 1 - c
    D = (c + 1)/2
    return u * np.sqrt((2*D)/c - B / (c*u**2) - 2/u - (2*np.log(u))/c)

delu = 0.001


def numerical_sltn(u0, n0, data_file, delu, numerical_type):
    # Load input data
    x = loadmat(data_file)
    
    # Get the speed of the wave
    a = x['A0'][0][0]
    c = (2 * a**2 * np.log(a) - a**2 + 1) / (a - 1)**2
   
    # Iterator
    i = 0
    # eta value
    n = n0
    # Function value
    u = u0

    # Discrete function values
    u_vec = np.array([u])
    # discrete eta values
    n_vec = np.array([n])

    while(n < 0):
        if numerical_type == "rk4":
            k1 = delu / func(u, c)
            k2 = delu / func(u + delu / 2, c)
            k3 = delu / func(u + delu/2, c)
            k4 = delu / func(u + delu, c)

            n = n + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u = u + delu
        elif numerical_type == "euler":
            n = n + delu/func(u, c)
            u = u + delu

        # Update values
        u_vec = np.append(u_vec, u)
        n_vec = np.append(n_vec, n)

    # Truncate the last nan value and shift so centered at 0
    n_vec = n_vec[0:-1] - n_vec[-2]
    u_vec = u_vec[0:-1] 

    return n_vec, u_vec

def numerical_sltn2(u0, n0, data_file, delu, numerical_type):
    # Load input data
    x = loadmat(data_file)
    
    # Get the speed of the wave
    #  a = x['A0'][0][0]
    #  c = (2 * a**2 * np.log(a) - a**2 + 1) / (a - 1)**2

    speed = speed_vec("./expData07.mat")
   
    # Iterator
    i = 0
    # eta value
    n = n0
    # Function value
    u = u0

    # Discrete function values
    u_vec = np.array([u])
    # discrete eta values
    n_vec = np.array([n])

    while(n < 0):
        c = speed[i]
        if numerical_type == "rk4":
            k1 = delu / func(u, c)
            k2 = delu / func(u + delu / 2, c)
            k3 = delu / func(u + delu/2, c)
            k4 = delu / func(u + delu, c)

            n = n + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u = u + delu
        elif numerical_type == "euler":
            n = n + delu/func(u, c)
            u = u + delu

        # Update values
        u_vec = np.append(u_vec, u)
        n_vec = np.append(n_vec, n)

    # Truncate the last nan value and shift so centered at 0
    n_vec = n_vec[0:-1] - n_vec[-2]
    u_vec = u_vec[0:-1] 

    return n_vec, u_vec




if __name__ == '__main__':

    speed = speed_vec("./expData07.mat")
    plt.plot(speed)
    plt.show()

    #  f0 = 1.0001
    #  n0 = -100
    #  df = 0.001
    #  #  n_vec, u_vec = numerical_sltn(f0, n0, "./expData07.mat", df, "rk4")
    #  #  n_vec2, u_vec2 = numerical_sltn(f0, n0, "./expData07.mat", df, "euler")
    #  #
    #  #  # Plot euler vs rk4 side by side
    #  #  fig, axs = plt.subplots(1, 2)
    #  #  axs[0].plot(n_vec, u_vec, label = "f' = +sqrt")
    #  #  axs[0].plot(-n_vec, u_vec, label = "f' = -sqrt")
    #  #  axs[0].set(xlabel = "eta", ylabel = "u(eta)")
    #  #  axs[0].set_title(f"RK4 - (f0, eta0) = ({f0}, {n0}), df = {delu}")
    #  #  axs[0].legend()
    #  #
    #  #  axs[1].plot(n_vec2, u_vec2, label = "f' = +sqrt")
    #  #  axs[1].plot(-n_vec2, u_vec2, label = "f' = -sqrt")
    #  #  axs[1].set(xlabel = "eta", ylabel = "u(eta)")
    #  #  axs[1].set_title(f"Euler - (f0, eta0) = ({f0}, {n0}), df = {df}")
    #  #  axs[1].legend()
    #  #
    #  #  plt.show()
    #
    #  # Plot only 1
    #  method = "euler"
    #  n_vec, u_vec = numerical_sltn(f0, n0, "./expData07.mat", df, method)
    #  plt.plot(n_vec, u_vec, label = "f' = +sqrt")
    #  plt.plot(-n_vec, u_vec, label = "f' = -sqrt")
    #  plt.xlabel("eta")
    #  plt.ylabel("u(eta)")
    #  plt.title(f"RK4 - (f0, eta0) = ({f0}, {n0}), df = {delu}")
    #  plt.legend()
    #
    #  plt.show()
    #
    #
    #
    #
