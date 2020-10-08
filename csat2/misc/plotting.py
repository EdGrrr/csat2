''' Plotting routines
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scipy.stats
import copy


def plt_bestfit(data_x, data_y, linec='r', stats=True, statspos='tl', *args, **kwargs):
    '''Plots a scatterplot and linear regression line, can cope with missing data'''
    plt.scatter(data_x, data_y, *args, **kwargs)
    ylim = np.array(plt.ylim())
    xlim = np.array(plt.xlim())
    mask = np.where(np.isfinite(data_x+data_y))
    a, b, r, p, t = scipy.stats.linregress(data_x[mask], data_y[mask])
    y_out = a*xlim+b
    x_out = (ylim-b)/a
    if y_out[0] < ylim[0]:
        if y_out[1] > ylim[1]:
            plt.plot(x_out, ylim, c=linec)
        else:
            plt.plot([x_out[0], xlim[1]], [ylim[0], y_out[1]], c=linec)

    else:
        if y_out[1] > ylim[1]:
            plt.plot([xlim[0], x_out[1]], [y_out[0], ylim[1]], c=linec)
        else:
            plt.plot(xlim, y_out, c=linec)
    if stats == True:
        xlim = np.array(plt.xlim())
        ylim = np.array(plt.ylim())
        if statspos == 'tl':
            x = xlim[0] + 0.05 * (xlim[1]-xlim[0])
            y = ylim[1] - 0.15 * (ylim[1]-ylim[0])
        elif statspos == 'bl':
            x = xlim[0] + 0.05 * (xlim[1]-xlim[0])
            y = ylim[0] + 0.05 * (ylim[1]-ylim[0])
        else:
            try:
                x = statspos[0]
                y = statspos[1]
            except:
                raise ValueError(
                    "Statpos should be 'tl', 'bl' or a two element tuple/list")
        text = 'Slope: %.4f\nIntercept: %.4f\nCorr: %.4f' % (a, b, r)
        plt.text(x, y, text)


def plt_sublabel(text, *args, **kwargs):
    '''Plots a sublabel in the top left corner of the plot, hopefully just outside the plotting area'''
    xlim, ylim = plt.xlim(), plt.ylim()
    plt.text(xlim[0], ylim[1]+0.02*(ylim[1]-ylim[0]), text, *args, **kwargs)


def plt_sublabel_index(i, *args, **kwargs):
    '''Will plot an alphabetic subplot label
    Starts with i=0 giving 'a)' '''
    alp = 'abcdefghijklmnopqrstuvwxyz'
    plt_sublabel(alp[i]+')', *args, **kwargs)


def plt_cbar(vmin, vmax, cmap, title, nticks=5, title_size=10,
             orientation='horizontal',
             plttype='heatmap', levels=None, aspect_ratio=0.03, norm=None, ticks=None,
             title_pos='above'):
    '''Plots a colorbar in an empty subplot
    Will do 'horizontal' or 'vertical' orientation
    Defaults to a heatmap colourbar, but if plttype='contourf' is passed,
    then a contourf-type colourbar will be plotted, using 'levels' to
    determine where the steps should be placed.   '''
    a = np.outer(np.arange(0, 1.005, 0.01), np.ones(int(aspect_ratio*100)))
    if orientation == 'horizontal':
        if plttype == 'heatmap':
            plt.imshow(a.transpose(), cmap=cmap, origin='upper',
                       interpolation='nearest', norm=norm)
        elif plttype == 'contourf':
            # Reset due to lack of interpolation
            a = np.outer(np.arange(0, 1.005, 0.01), np.ones(4))
            if levels == None:
                levels = np.linspace(0, 1, 10)
            else:
                levels = (levels-vmin)/float(vmax-vmin)
            plt.contourf(a.transpose(), levels=levels, cmap=cmap)
            ax = plt.gca()
            ax.set_aspect(1)
        plt.yticks([], [])
        if ticks:
            plt.xticks(np.linspace(0, 100, nticks), ticks, size=title_size)
        else:
            plt.xticks(np.linspace(0, 100, nticks), np.linspace(
                vmin, vmax, nticks), size=title_size)
        if title_pos == 'above':
            plt.title(title,size=title_size)
        elif title_pos == 'below':
            plt.xlabel(title,size=title_size)
    elif orientation == 'vertical':
        if plttype == 'heatmap':
            plt.imshow(a, cmap=cmap, interpolation='nearest',
                       origin='lower', norm=norm)
        elif plttype == 'contourf':
            a = np.outer(np.arange(0, 1.005, 0.01), np.ones(4))
            if levels == None:
                levels = np.linspace(0, 1, 10)
            else:
                levels = (levels-vmin)/float(vmax-vmin)
            plt.contourf(a, levels=levels, cmap=cmap)
            ax = plt.gca()
            ax.set_aspect(1.2)
        plt.xticks([], [])
        plt.tick_params(axis='y', which='both',
                        labelleft=False, labelright=True)
        if ticks:
            plt.yticks(np.linspace(0, 100, nticks), ticks, size=title_size)
        else:
            plt.yticks(np.linspace(0, 100, nticks), np.linspace(
                vmin, vmax, nticks), size=title_size)
        plt.tick_params(axis='y', which='both',
                        labelleft=False, labelright=True)
        plt.ylabel(title, size=title_size)
    return


####################
# Colorbars
####################
# NCL style WhiteBlueGreenYellowRed linear segmented colourmap
_wbgyr_cdict = {'red': ((0.0,  255./255, 255./255),
                        (0.125, 173./255, 173./255),
                        (0.25,  95./255,  95./255),
                        (0.375, 73./255,  73./255),
                        (0.5,  165./255, 164./255),
                        (0.625, 248./255, 248./255),
                        (0.75, 236./255, 236./255),
                        (0.875, 200./255, 200./255),
                        (1.0,  146./255, 146./255)),
                'green': ((0.0,  255./255, 255./255),
                          (0.125, 224./255, 224./255),
                          (0.25, 163./255, 163./255),
                          (0.375, 166./255, 166./255),
                          (0.5,  207./255, 207./255),
                          (0.625, 184./255, 184./255),
                          (0.75,  86./255,  86./255),
                          (0.875, 29./255,  29./255),
                          (1.0,   21./255,  21./255)),
                'blue': ((0.0,  255./255, 255./255),
                         (0.125, 248./255, 248./255),
                         (0.25, 214./255, 214./255),
                         (0.375, 120./255, 120./255),
                         (0.5,   81./255,  81./255),
                         (0.625, 73./255,  73./255),
                         (0.75,  41./255,  41./255),
                         (0.875, 38./255,  38./255),
                         (1.0,   25./255,  25./255))}
# And reverse
_wbgyr_cdict_r = copy.deepcopy(_wbgyr_cdict)
for i in _wbgyr_cdict_r.keys():
    _wbgyr_cdict_r[i] = [(1-j[0], j[1], j[2]) for j in _wbgyr_cdict_r[i]]
    _wbgyr_cdict_r[i].reverse()

cmap = LinearSegmentedColormap('WBGYR', _wbgyr_cdict)
cmap.set_bad('#D2D2D2')
plt.register_cmap(cmap=cmap)

cmap = LinearSegmentedColormap('WBGYR_r', _wbgyr_cdict_r)
cmap.set_bad('#D2D2D2')
plt.register_cmap(cmap=cmap)
