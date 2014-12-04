#!/bin/csh

rm -f _________plot_all_snapshots.gnu

#echo "set term eps" >> _________plot_all_snapshots.gnu

echo set xrange \[0:50000\] >> _________plot_all_snapshots.gnu
echo "#set yrange [-1.5:+1.5]" >> _________plot_all_snapshots.gnu

foreach file ( ./OUTPUT_FILES/snapshot* )

  echo plot \"$file\" w l lc 1 >> _________plot_all_snapshots.gnu
  echo "pause -1 'hit any key...'" >> _________plot_all_snapshots.gnu

end

gnuplot _________plot_all_snapshots.gnu
rm -f _________plot_all_snapshots.gnu

