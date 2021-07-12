import matplotlib.pyplot as plt
import numpy as np

def main():
    
    ncv = 100

    filename = "lightcurves\lcvold%03d.dat" % ncv

    lcv_old = np.loadtxt(filename,delimiter=',')
    tv  = lcv_old[:,0]
    lum = lcv_old[:,1]
    fig = plt.scatter(tv,lum,label='old',s=10)

    lcv_new = np.loadtxt('lightcurves\lcvnew001.dat',delimiter=',')
    tv  = lcv_new[:,0]
    lum = lcv_new[:,1]

    fig = plt.scatter(tv,lum,label='new',s=10)
    fig.axes.set_xlabel('Time (min)')
    fig.axes.set_ylabel('Luminosity')

    plt.ylim((0.9,1))
    plt.legend()
    plt.show()

if __name__ == '__main__':
	main()