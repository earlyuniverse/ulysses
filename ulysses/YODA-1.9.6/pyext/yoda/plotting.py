# -*- python -*-

"""
Plotting utilities, particularly for interaction with matplotlib and Rivet make-plots
"""

import yoda
import sys
import numpy as np
import matplotlib as mpl


# TODO: Move to core objects
#     def same_binning_as(self, other):
#         if self.dim != other.dim:
#             return False
#         if not (other.x == self.x).all() and \
#                (other.exminus == self.exminus).all() and \
#                (other.explus == self.explus).all():
#             return False
#         if self.dim == 2:
#             return True
#         return (other.y == self.y).all() and \
#                (other.eyminus == self.eyminus).all() and \
#                (other.eyplus == self.eyplus).all()


def read_plot_keys(datfile):
    import re
    re_begin = re.compile("#*\s*BEGIN\s+PLOT\s*(\w*)")
    re_comment = re.compile("#.*")
    re_attr = re.compile("(\w+)\s*=\s*(.*)")
    re_end = re.compile("#*\s*END\s+PLOT\s+\w*")
    plotkeys = {}
    with open(datfile) as f:
        inplot = False
        name = None
        for line in f:
            l = line.strip()
            if re_begin.match(l):
                inplot = True
                name = re_begin.match(l).group(1)
            elif re_end.match(l):
                inplot = False
                name = None
            elif re_comment.match(l):
                continue
            elif inplot:
                m = re_attr.match(l)
                if m is None: continue
                plotkeys.setdefault(name, {})[m.group(1)] = m.group(2)
    return plotkeys


def mplinit(engine="MPL", font="TeX Gyre Pagella", fontsize=17, mfont=None, textfigs=True):
    """One-liner matplotlib (mpl) setup.

    By default mpl will be configured with its native MathText rendering
    backend, and a Palatino-like font for both text and math contexts, using
    'lower-case numerals' if supported. Setting the engine to 'TEX' will use
    standard mpl rendering, with calls to LaTeX for axis labels and other text;
    setting it to 'PGF' will use the TeX PGF renderer: both these modes are much
    slower than MPL mode, but the latter only supports a limited set of LaTeX
    macros and does not render as nicely as the TeX backends.

    The font and mfont optional arguments can be used to choose a different text
    font and math font respectively; if mfont is None, it defaults to the same
    as the text font. The textfigs boolean argument can be set false to disable
    the lower-case/text/old-style numerals and use 'upper-case' numerals
    everywhere. These options do not currently apply to the MPL rendering engine.
    """
    mpl.rcParams.update({
        "text.usetex" : (engine != "MPL"),
        "font.size"   : int(fontsize),
        "font.family" : "serif", #< TODO: make configurable? auto-detect?
        })

    texpreamble = [r"\usepackage{amsmath,amssymb}", r"\usepackage{mathspec}"]
    mfont = mfont if mfont else font
    fontopts = "[Numbers=OldStyle]" if textfigs else ""
    mfontopts = fontopts.replace("]", ",") + "Scale=MatchUppercase" + "]"
    texpreamble.append( r"\setmainfont{fopts}{{{font}}}".format(fopts=fontopts, font=font) )
    texpreamble.append( r"\setmathsfont(Digits,Latin){fopts}{{{font}}}".format(fopts=mfontopts, font=mfont) )

    if engine.upper() == "PGF":
        mpl.use("pgf")
        mpl.rcParams["pgf.preamble"] = texpreamble
    # TODO: Fix?
    # elif engine.upper() == "TEX":
    #     mpl.rcParams["text.latex.preamble"] = texpreamble

    return mpl

## Alias
initmpl = mplinit
setup_mpl = mplinit


def show():
    """
    Convenience call to matplotlib.pyplot.show()

    NOTE: done this way to avoid import of pyplot before mplinit()
    or mpl.use() has been (optionally) called.
    """
    import matplotlib.pyplot as plt
    plt.show()


def mk_figaxes_1d(ratio=True, title=None, figsize=(8,6)):
    "Make a standard main+ratio plot figure and subplot layout"

    ## We need to use pyplot here to set up the backend-specific canvas
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=figsize)
    #fig = mpl.figure.Figure(figsize=figsize, tight_layout=True)

    if title:
        fig.suptitle(title, horizontalalignment="left", x=0.13)

    ## Make axes. GridSpec may not be available, in which case fall back ~gracefully
    axmain, axratio = None, None
    if ratio:
        try:
            gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3,1], hspace=0)
            axmain = fig.add_subplot(gs[0])
            #axmain.hold(True)
            axratio = fig.add_subplot(gs[1], sharex=axmain)
            #axratio.hold(True)
            axratio.axhline(1.0, color="gray") #< Ratio = 1 marker line
        except:
            sys.stderr.write("matplotlib.gridspec not available: falling back to plotting without a ratio\n")
            ratio = False
    if not ratio:
        axmain = fig.add_subplot(1,1,1)
        #axmain.hold(True)

    return fig, (axmain, axratio)


def set_axis_labels_1d(axmain, axratio, xlabel=None, ylabel=None, ratioylabel=None):
    axmain.set_ylabel(ylabel, y=1, ha="right", labelpad=None)
    if axratio:
        axmain.xaxis.set_major_locator(mpl.ticker.NullLocator())
        axratio.set_xlabel(xlabel, x=1, ha="right", labelpad=None)
        axratio.set_ylabel(ratioylabel)
    else:
        axmain.set_xlabel(xlabel, x=1, ha="right", labelpad=None)


def mk_lowcase_dict(d):
    "Convert the keys of a str->obj dict to lower-case"
    return dict((k.lower(), v) for (k,v) in d.items())


# TODO: Needs generalisation for 2D marginal axes)
def setup_axes_1d(axmain, axratio, **plotkeys):

    ## Case-insensitize the plotkeys dict
    plotkeys = mk_lowcase_dict(plotkeys)

    ## Axis labels first
    xlabel = plotkeys.get("xlabel", "")
    ylabel = plotkeys.get("ylabel", "")
    ratioylabel = plotkeys.get("ratioylabel", "ratio")
    set_axis_labels_1d(axmain, axratio, xlabel, ylabel, ratioylabel)

    ## log/lin measures
    # TODO: Dynamic default based on data ranges?
    # TODO: take log axes and preference for round numbers into account in setting default axis limits
    xmeasure = "log" if yoda.util.as_bool(plotkeys.get("logx", False)) else "linear"
    ymeasure = "log" if yoda.util.as_bool(plotkeys.get("logy", False)) else "linear"
    ratioymeasure = "log" if yoda.util.as_bool(plotkeys.get("ratiology", False)) else "linear"
    axmain.set_xscale(xmeasure)
    axmain.set_yscale(ymeasure)
    if axratio:
        axratio.set_xscale(xmeasure)
        axratio.set_yscale(ratioymeasure)

    ## Plot range limits
    if "ymin" in plotkeys:
        axmain.set_ylim(bottom=float(plotkeys.get("ymin")))
    if "ymax" in plotkeys:
        axmain.set_ylim(top=float(plotkeys.get("ymax")))
    #
    if "xmin" in plotkeys:
        axmain.set_xlim(left=float(plotkeys.get("xmin")))
    if "xmax" in plotkeys:
        axmain.set_xlim(right=float(plotkeys.get("xmax")))
    #
    if axratio:
        # TODO: RatioSymmRange option
        # axratio.set_xlim([xmin-0.001*xdiff, xmax+0.001*xdiff]) # <- TODO: bad on a log scale!
        if "xmin" in plotkeys:
            axratio.set_xlim(left=float(plotkeys.get("xmin")))
        if "xmax" in plotkeys:
            axratio.set_xlim(right=float(plotkeys.get("xmax")))
        if "ratioymin" in plotkeys:
            axratio.set_ylim(bottom=float(plotkeys.get("ratioymin")))
        if "ratioymax" in plotkeys:
            axratio.set_ylim(top=float(plotkeys.get("ratioymax")))

    # TODO: Ratio plot manual ticks


def plot_hist_on_axes_1d(axmain, axratio, h, href=None, default_color="black", default_linestyle="-", **plotkeys):

    ## Case-insensitize the plotkeys dict
    hkeys = mk_lowcase_dict(h.annotationsDict())
    hkeys.update(plotkeys)
    plotkeys = hkeys

    # TODO: Split into different plot styles: line/filled/range, step/diag/smooth, ...?

    ## Styles
    default_color = plotkeys.get("color", default_color)
    marker = plotkeys.get("marker", plotkeys.get("polymarker", None)) # <- make-plots translation
    marker = {"*":"o"}.get(marker, marker) # <- make-plots translation
    mcolor = plotkeys.get("linecolor", default_color)
    errbar = plotkeys.get("errorbars", None)
    ecolor = plotkeys.get("errorbarscolor", default_color)
    line = plotkeys.get("line", None)
    lcolor = plotkeys.get("linecolor", default_color)
    lstyle = plotkeys.get("linestyle", default_linestyle)
    lstyle = {"solid":"-", "dashed":"--", "dotdashed":"-.", "dashdotted":"-.", "dotted":":"}.get(lstyle, lstyle) # <- make-plots translation
    lwidth = 1.4
    msize = 7

    ## If no drawing is enabled, default to a step line
    if not any([marker, line, errbar]):
        line = "step"

    ## MPL plotting
    # TODO: Split this into different functions for each kind of data preparation (and smoothing as an extra function?)
    # TODO: First convert h to scatter
    artists = None
    if errbar:
        artists = axmain.errorbar(h.xVals(), h.yVals(), xerr=h.xErrs(), yerr=h.yErrs(), color=ecolor, linestyle="none", linewidth=lwidth, capthick=lwidth) # linestyle="-", marker="o",
    if line == "step":
        artists = axmain.step(np.append(h.xMins(), h.xMax()), np.append(h.yVals(), h.yVals()[-1]), where="post", color=lcolor, linestyle=lstyle, linewidth=lwidth)
    elif line == "diag":
        artists = axmain.plot(h.xVals(), h.yVals(), color=lcolor, linestyle=lstyle, linewidth=lwidth)
    elif line == "smooth":
        from scipy.interpolate import spline
        xnew = np.linspace(min(h.xVals()), max(h.xVals()), 3*h.numBins)
        ynew = spline(h.xVals(), h.yVals(), xnew)
        artists = axmain.plot(xnew, ynew, color=lcolor, linestyle=lstyle, linewidth=lwidth)
    if marker:
        artists = axmain.plot(h.xVals(), h.yVals(), marker=marker, markersize=msize, linestyle="none", color=mcolor, markeredgecolor=mcolor)

    ## Legend entry
    if h.annotation("Title") and artists:
        artists[0].set_label(h.annotation("Title"))

    ## Ratio
    ratioartists = None
    if href and h is not href:
        # TODO: exclude and specify order via RatioIndex
        # assert h.same_binning_as(href)
        # TODO: log ratio or #sigma deviation
        yratios = np.array(h.yVals())/np.array(href.yVals())
        # TODO: Same styling control as for main plot (with Ratio prefix, default to main plot style)
        ## Stepped plot
        ratioartists = axratio.step(np.append(href.xMins(), href.xMax()), np.append(yratios, yratios[-1]), where="post", color=lcolor, linestyle=lstyle, linewidth=lwidth)
        # TODO: Diag plot
        # axratio.plot(href["x"], yratios, color="r", linestyle="--")
        # TODO: Smoothed plot

    return artists


def plot(hs, outfile=None, ratio=True, show=False, axmain=None, axratio=None, **plotkeys):
    """
    Plot the given histograms on a single figure, returning (fig, (main_axis,
    ratio_axis)). Show to screen if the second arg is True, and saving to outfile
    if it is otherwise non-null.
    """

    ## Case-insensitize the plotkeys dict
    plotkeys = mk_lowcase_dict(plotkeys)

    ## Handle single histo args
    if isinstance(hs, yoda.AnalysisObject):
        hs = [hs,]
        ratio = False

    ## Get data ranges (calculated or forced)
    xmin = float(plotkeys.get("xmin", min(h.xMin() for h in hs)))
    xmax = float(plotkeys.get("xmax", max(h.xMax() for h in hs)))
    xdiff = xmax - xmin
    # print xmin, xmax, xdiff
    # TODO: Tweak max-padding for top tick label... sensitive to log/lin measure
    ymin = plotkeys.get("ymin", min(min(h.yVals()) for h in hs))
    #print( max(max(h.yVals()) for h in hs) )
    ymax = plotkeys.get("ymax", 1.1*max(max(h.yVals()) for h in hs))
    ymin = float(ymin)
    ymax = float(ymax)
    ydiff = ymax - ymin
    # print ymin, ymax, ydiff

    ## Identify reference histo by annotation (unless explicitly disabled)
    href = None
    # TODO: Use ratio to setdefault RatioPlot in plotkeys, then use that to decide whether to look for href
    if ratio:
        for h in hs:
            hkeys = mk_lowcase_dict(h.annotationsDict())
            if yoda.util.as_bool(hkeys.get("ratioref", False)):
                if href is None:
                    href = h
                else:
                    #print "Multiple ratio references set: using first value = {}".format(href.path)
                    break
    if href is None: #< no ref found -- maybe all were disabled?
        ratio = False

    ## Make figure and subplot grid layout
    title = plotkeys.get("title", "")
    if not axmain:
        fig, (axmain, axratio) = mk_figaxes_1d(ratio and not axratio, title)
    else:
        fig = axmain.get_figure()

    ## Setup axes appearances
    axmain.set_xlim([xmin, xmax])
    axmain.set_ylim([ymin, ymax])
    if axratio:
        axratio.set_xlim([xmin, xmax])
        axratio.set_ylim(auto=True)
    setup_axes_1d(axmain, axratio, **plotkeys)

    # TODO: specify ratio display in log/lin, abs, or #sigma, and as x/r or (x-r)/r

    ## Draw ratio error band (do this before looping over cmp lines)
    # TODO: Actually we can call this when we hit the href, and force the zorder into groups: bands, lines, dots, legend, text, frame
    if axratio:
        ref_ymax_ratios = np.array(href.yMaxs())/np.array(href.yVals())
        ref_ymin_ratios = np.array(href.yMins())/np.array(href.yVals())
        # TODO: Diag: (needs -> limit handling at ends)
        # axratio.fill_between(href.x, ref_ymin_ratios, ref_ymax_ratios, edgecolor="none", facecolor=ratioerrcolor, interpolate=False)
        # Stepped:
        def xedges_dbl(h):
            edges = np.empty((2*len(h.xVals()),))
            edges[0::2] = h.xMins()
            edges[1::2] = h.xMaxs()
            return edges
        def dbl_array(arr):
            return sum(([x,x] for x in arr), [])
        ratioerrcolor = plotkeys.get("ratioerrcolor", "yellow")
        axratio.fill_between(xedges_dbl(href), dbl_array(ref_ymin_ratios), dbl_array(ref_ymax_ratios),
                             edgecolor="none", facecolor=ratioerrcolor)
        # TODO: Smoothed: (needs -> limit handling at ends)
        # Redraw ratio = 1 marker line:
        axratio.axhline(1.0, color="gray")

    COLORS = ["red", "blue", "magenta", "orange", "green"]
    LSTYLES = ["-", "--", "-.", ":"]

    ## Dataset plotting
    some_valid_label = False
    for ih, h in enumerate(hs):
        #print ih, h.path
        aa = plot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[ih % len(COLORS)], LSTYLES[ih % len(LSTYLES)])
        if aa and not aa[0].get_label().startswith("_"):
            # print "@@@", aa[0].get_label()
            some_valid_label = True

    ## Legend
    # TODO: allow excluding and specify order via LegendIndex
    if some_valid_label: #< No point in writing a legend if there are no labels
        pass #axmain.legend(loc=plotkeys.get("LegendPos", "best"), fontsize=plotkeys.get("LegendFontSize", "x-small"), frameon=False)

    ## Tweak layout now that everything is in place
    # TODO: merge tight_layout() into the Figure constructor, and maybe the ratio ticker when retrospectively drawing the zorder'ed err band
    if axratio:
        axratio.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4, prune="upper"))
    fig.tight_layout()

    ## Save to an image file if we were asked to
    if outfile:
        #print "Saving to " + outfile
        fig.savefig(outfile)

    ## Show to screen if requested
    if show:
        import matplotlib.pyplot as plt
        plt.show()

    ## Return the figure objects
    return fig, (axmain, axratio)


## Aliases
plot_hists_1d = plot
plot_hist_1d = plot


def _plot1arg(args):
    "Helper function for mplot, until Py >= 3.3 multiprocessing.pool.starmap() is available"
    return plot(*args)


def nplot(hs, outfiles=None, ratio=True, show=False, nproc=1, **plotkeys):
    """
    Plot the given list of histogram(s), cf. many calls to plot().

    hs must be an iterable, each entry of which will be the content of a single
    plot: the entries can either be single histograms or lists of histograms,
    i.e. either kind of valid first argument to plot().

    Outfiles must be an iterable corresponding to hs, and ratio may either be a
    bool or such an iterable.

    The return value is a list of the return tuples from each call to plot(), of
    the same length as the hs arg.


    MULTIPROCESSING -- *WARNING* CURRENTLY BROKEN

    The main point of this function, other than convenience, is that the Python
    multiprocessing module can be used to distribute the work on to multiple
    parallel processes.

    The nproc argument should be the integer number of parallel processes on
    which to distribute the plotting. nproc = None (the default value) will use
    Ncpu-1 or 1 process, whichever is larger. If nproc = 1, multiprocessing will
    not be used -- this avoids overhead and eases debugging.
    """

    argslist = []
    for i, hs_arg in enumerate(hs):
        outfile_arg = outfiles[i] if outfiles else None
        ratio_arg = ratio[i] if hasattr(ratio, "__iter__") else ratio
        show_arg = False #< we just do this once, at the end
        plotkeys_arg = plotkeys if type(plotkeys) is dict else plotkeys[i]
        argslist.append( (hs_arg, outfile_arg, ratio_arg, show_arg, None, None, plotkeys_arg) )
    #print argslist

    # TODO: make the multiprocessing work
    import multiprocessing
    nproc = nproc or multiprocessing.cpu_count() or 1
    if nproc > 1:
        pool = multiprocessing.Pool(processes=nproc)
        res = pool.map_async(_plot1arg, argslist)
        rtn = res.get()
    else:
        ## Run this way in the 1 proc case for easier debugging
        rtn = [_plot1arg(args) for args in argslist]

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return rtn
