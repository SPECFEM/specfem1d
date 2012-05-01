#!/bin/csh

# run this script using " ./create_gnuplot_script_to_display_seismograms.csh > plot_all_snapshots.gnu "

echo "set term x11"
echo "#set term wxt"

echo set xrange \[0:3000\]
echo "#set yrange [-1.5:+1.5]"

foreach file ( snapshot* )

  echo plot \"$file\" w l lc 1
  echo "pause -1 'hit any key...'"

end

