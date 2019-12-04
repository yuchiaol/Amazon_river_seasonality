import numpy as np
import scipy.fftpack as fftp
import matplotlib.pyplot as plt    

def fourier_series(x, y, wn, n=None):
    # get FFT
    myfft = fftp.fft(y, n)
    # kill higher freqs above wavenumber wn
    myfft[wn:-wn] = 0
    # make new series
    y2 = fftp.ifft(myfft)

    plt.figure(num=None)
    plt.plot(x, y)
    plt.plot(x, y2)
    plt.show()

if __name__=='__main__':
    x = np.array([float(i) for i in range(0,360)])
    y = np.sin(2*np.pi/360*x) + np.sin(2*2*np.pi/360*x) + 5

#    fourier_series(x, y, 3, 360)

    n = 360
    wn = 3
    myfft = fftp.fft(y)
# kill higher freqs above wavenumber wn
#    myfft[wn:-wn] = 0
# make new series
    y2 = fftp.ifft(myfft*2)
    y_new = (y2-y2.mean()+y.mean()).real

    plt.figure(num=None)
    plt.plot(x, y, 'b-', label='original')
    plt.plot(x, y_new, 'r-', label='amplified')
    plt.xlabel('x')
    plt.legend()

    plt.savefig('fft_seasonality_tmp_plot.jpg', format='jpeg', dpi=200)


    plt.show()

