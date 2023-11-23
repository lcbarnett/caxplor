# caxplor
Highly-efficient simulation, filtering and visualisation of 1D binary-state [cellular automata](https://en.wikipedia.org/wiki/Cellular_automaton).

### Building
This code requires a 64-bit little-endian architecture, and uses [X11/Xlib](https://www.x.org/releases/current/doc/libX11/libX11/libX11.html) for graphics. As yet, it has only been built and tested on Linux x86-64, but is in principle portable to MacOS with an X server, e.g.,  [XQuartz](https://www.xquartz.org/), or Windows with [WSL](https://learn.microsoft.com/en-us/windows/wsl/), [Cygwin](https://www.cygwin.com/) or an X server like [XMing](http://www.straightrunning.com/XmingNotes/) [^1]. It also reguires the [GD graphics library](https://libgd.github.io/pages/about.html); if you are on Linux, install the appropriate development package through your software manager.

To build, you will need the [Make](https://www.gnu.org/software/make/) build tool. In a terminal, navigate to the caxplor root directory and type 'make' to build. There is no installation; the executable is called 'caxplor'.

### Usage
To run the main 'CA explorer' routine with default parameters, type
```
./caxplor xplor
```
at the terminal prompt. This will display a random 1D CA with rule of size 5 [^2] in full(ish) screen, and display a list of available simulation parameters and a list of key commands (to be typed in the focused CA window, not the terminal!) The 'f' key pages the CA one screen forward, and 'i' re-initialises the CA. The simulation maintains a list of CA rules, and for each CA rule, a list of filter rules [^3]. There are two operational/display modes: 'exploring' and 'filtering'; the ESC (or 'm') key toggles the mode. In filtering mode, the CA screen displays in dark blue, and operations apply to filters on the current CA. Diagnostic information appears on the terminal, in particular the current CA/filter ID (as a hex string representation of the rule table) and Langton's &lambda;-parameter.

New random CA/filter rules are added to the corresponding list with SPACE (or 'n'); 'N' prompts at the terminal for a user-supplied CA/filter ID. You may scroll through the CA/filter lists with the left- and right-arrow (or 'j', 'k') keys, while 'J' and 'K' take you to the first and last CA/filter respectively. The current CA/filter may be deleted with the DEL (or 'd') key. The 'c' key prompts to reset the size of subsequent new CA/filter rules. The current CA/filter rule ID(s) may be saved by pressing 's', by default in a text file called 'saved.rt'. Screen images may be written to file (by default in PNG format) with the 'w' key. Press 'h' to redisplay the key command help, and 'q' to quit the simulation.

Simulation parameters may be set using command-line switches; e.g.,
```
./caxplor -rsiz 7 -fsiz 3 -rlam 0.3
```
will run CAs of size 7, filters of size 3 and an average Langton's &lambda; = 0.3. To see the available switches and their values, just run, e.g., `./caxplor xplor -i <other switches>`.

Saved CA rules (and their filter rules, if present) may also be loaded from a file:
```
./caxplor -irtfile my_saved_rules.rt
```
The entropy and 1-lag [transfer entropy](https://link.springer.com/book/10.1007/978-3-319-43222-9) aka [dynamical dependence](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.014304) for the current CA/filter may be calculated with the 'E' and 'D' keys respectively. This (experimental and undocumented) feature requires the [Gnuplot](http://www.gnuplot.info/) scientific graphing utility to be installed on your system.

There are currently a few (probably buggy/undocumented) routines for analysis and benchmarking and batch dynamical independence calculation, as well as a template for your own test routines, which may be run as `./caxplor ana`, `./caxplor bmark`, `./caxplor ddr` and `./caxplor test` respectively; you may edit these to taste.

Have fun!

### Remarks
For efficiency, CA simulation is performed at the bit level in 64-bit words. This means that the width of the CA (the number of cells in a row) may only be multiples of 64. The switch `-nwords` of the `xplor` simulation specifies the number of 64-bit words (if set to 0, the number of words is set to maximise screen width of the display). CAs wrap at the beginning/end of rows.

### Author
Lionel Barnett: lionelb@sussex.ac.uk

[^1]: Do let the author know (or fork and contribute) if you have any joy porting to other platforms.
[^2]: CA rules are specified by a rule table, which maps each short binary sequence of specified length to an associated binary value. The "size" of a rule is the length of the binary sequence operated on. A CA rule operates on a long binary sequence of specified length (the "width" of a CA "row") by sliding along the sequence, thus generating the cells of the output row. This process is applied recursively, with the output row becoming the input row for the next iteration.
[^3]: A filter rule is itself a CA rule, except that instead of being applied recursively, it is applied independently to each row of an existing CA. Filter rules need not be the same size as the CA rule which generated the CA.
