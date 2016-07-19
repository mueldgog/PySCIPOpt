from graphics import *
from random import randrange, seed

def displaySolution(data, edgemap):
    """
    Displays a solution of the semifluid packing problem.
    Input:
    data    - instance data
    edgemap - list containing all edge tuples (e,type,x,y)
    """
    D = data['D']
    H = data['H']
    l = data['l']
    h = data['h']

    seed(5)

    colors = ['#488101', '#ac6a77', '#1e99c2','#aa20d2', '#363e5a', '#fc6e9d','#879f0e', '#1679f9', '#d4ee4d','#557e23']
    for i in range(len(h)-10):
        colors.append(color_rgb(randrange(0,255), randrange(0,255), randrange(0,255)))

    scale = 500.0 / max(D,H)
    win = GraphWin('solution', D * scale + 100, H * scale + 100)

    container = Rectangle(Point(50,50), Point(D*scale + 50, H*scale + 50))
    container.draw(win)

    V    = [0]
    left = {0 : 0}
    up   = {0 : 0}

    for (e,t,xval,yval) in edgemap:
        # skip unused edges <=> y_e = 0
        if yval <= 0.5:
            continue

        assert e[0] in left and e[0] in up

        lt = l[t]
        ht = h[t] * xval
        rect = Rectangle(Point(scale*left[e[0]] + 50, H*scale + 50 - scale*up[e[0]]), Point(scale*(left[e[0]] + lt) + 50, H*scale + 50 - scale*(up[e[0]] + ht)))
        rect.setFill(colors[t])
        rect.draw(win)
        left[e[1]] = left[e[0]]
        up[e[1]] = up[e[0]] + ht
        left[e[0]] = left[e[0]] + lt

    win.getMouse()
    win.close()
