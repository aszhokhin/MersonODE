#
# $Id: using.demo 3.38.2.6 1992/11/14 02:25:21 woo Exp $
#
# Requires data file "using.dat" from this directory,
# so change current working directory to this directory before running.
#


#set term post eps enhanced "Times-Roman" 18
#set output "d4eqg.ps"
#set term gif
#set output "d4eqg.gif"

set size 1.,.8
set xlabel '{/Times-Italic=20 Time}'
set ylabel '{/Times-Italic=20 Concentrations}'

plot 'D4EQ.DAT' using 1:2 title '{/Times-Italic 1}' with lines,\
     'D4EQ.DAT' using 1:3 title '{/Times-Italic 2}' with lines,\
     'D4EQ.DAT' using 1:4 title '{/Times-Italic 3}' with lines,\
     'D4EQ.DAT' using 1:5 title '{/Times-Italic 4}' with lines
pause -1 "Hit return to continue"
set size .8,.9

plot 'D4EQ.DAT' using 2:3 title '{/Times-Italic 1}' with lines
pause -1 "Hit return to continue"

