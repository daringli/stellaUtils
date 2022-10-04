#!/usr/bin/env python

import numpy as np
import os

def parse(filename):
    with open(filename) as f:
        lines = f.readlines()

    ts = []
    phi2s = []
    for l in lines:
        ts.append(float(l.split('=')[2].strip().split()[0]))
        phi2s.append(float(l.split('=')[3].strip().split()[0]))
        

    return np.array(ts),np.array(phi2s)

if __name__=="__main__":
    import sys
    import matplotlib.pyplot as plt


    lines = []
    def onclick(event):
        gammas  = []
        clicked_t = event.xdata
        for ts,phi2s in zip(tss,phi2ss):
            i  = (np.abs(ts-clicked_t)).argmin()
            y  = np.log(phi2s[i:])
            x = ts[i:]
            coefs = np.polyfit(x,y,1)
            gammas.append(coefs[0]/2)
        text1.set_text(gammas)
        
        line = plt.axvline(clicked_t)
        lines.append(line)
        if len(lines) > 1:
            lines[0].remove()
            del lines[0]
        plt.draw()
        
    

    argv  = sys.argv
    argc  = len(argv)
    if argc > 1:
        filenames = argv[1:]
    else:
        filenames = ["stella.out"]
    

    fig,ax = plt.subplots()


    text1 = ax.text(0.5, 0.995,  "", transform=None)
    

    tss = []
    phi2ss = []
    for filename in filenames:
        if os.path.isdir(filename):
            filename = filename + "/stella.out"
        ts,phi2s  =parse(filename)
        tss.append(ts)
        phi2ss.append(phi2s)
        ax.semilogy(ts,phi2s)
    ax.legend(filenames)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)


    plt.show()
