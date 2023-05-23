# caxplor
Highly-efficient simulation, filtering and visualisation of 1D binary-state cellular automata.

This code requires a 64-bit little-endian architecture, and uses [X11](https://en.wikipedia.org/wiki/X_Window_System) graphics. As yet, it has only been built and tested on Linux x86-64, but is in principle portable to MacOS with an X server, e.g.,  [XQuartz](https://www.xquartz.org/), or Windows with [WSL](https://learn.microsoft.com/en-us/windows/wsl/), [Cygwin](https://www.cygwin.com/) or an X server like [XMing](http://www.straightrunning.com/XmingNotes/) [^1].

### Building
To build, you will need the [Make](https://www.gnu.org/software/make/) build tool. In a terminal, navigate to the caxplor root directory and type 'make' to build. There is no installation; the executable is called 'caxplor'.

### Usage
To run the main 'CA explorer' routine with default parameters, type
```
./caxplor
```
This will display a random 5th-order 1D CA in full(ish) screen, and display a list of available parameters and a list of key commands - to be typed in the CA window, not the terminal! 'f' pages the CA one screen forward, and 'i' re-initialises the CA. The routine maintains a list of CA rules, and for each CA rule, a list of filter rules. There are two operational/display modes: 'exploring' and 'filtering'. In filtering mode, the CA screen displays in dark blue, and operations apply to filters on the current CA. Diagnostic information appears on the terminal, in particular the current CA/filter ID (as a hex string representation of the rule table) and Langton's &lambda;-parameter. New random CAs/filters are added to the corresponding list with 'n'; 'N' prompts at the terminal for a user-supplied CA/filter ID. You may scroll through the CA/filter lists with the left- and right-arrow (or 'j', 'k') keys, and the current CA/filter may be deleted with the DEL (or 'd') key. The current CA/filter rule ID(s) may be saved by pressing 's', by default in a text file called 'saved.rt'. Screen images may be written to file (by default in PNG format) with the 'w' key. Press 'h' to redisplay the key command help.

Simulation parameters may be set using command-line switches; e.g.,
```
./caxplor -B 7 -F 3 -rlam 0.3
```
will run CAs of order 7, filters of order 3 and an average CA &lambda; equal to 0.2.

There are a currently a couple of (probably buggy/undocumented) routines for analysis and benchmarking, as well as a template for test routines, which may be run as `./caxplor ana`, `.caxplor bmark` and `.caxplor bmark` respectively; you may edit these to taste.

Have fun!

### Author
Lionel Barnett: lionelb@sussex.ac.uk

[^1]: Do let the author know (or fork and contribute) if you have any joy porting to other platforms.
