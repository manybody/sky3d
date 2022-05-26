#!/usr/bin/env bash

# automatically make animations in directories with *.tdd files
# - takes one argument, which is the name of a directory.  By default
# the final movie file has the directory name as its own base name, so
# it is better not to use "." as the directory (but it's not fatal).

# Things to-do:  parse for005 file to e.g. set geometry of plot up better

# USER NEEDS TO DEFINE LOCATION OF CUTS PROGRAM
CUTS=$HOME/Sky3D/Utils/Cuts/Cuts  

# Check gnuplot version
GP_MAJ=`gnuplot -V |awk '{print $2}' |awk -F. '{print $1}'`
GP_MIN=`gnuplot -V |awk '{print $2}' |awk -F. '{print $2}'`

if [ $GP_MAJ -lt 4 ] 
then
    echo At least v4.6 of gnuplot needed.  You have ${GP_MAJ}.$GP_MIN
    exit
fi
if [ $GP_MAJ -eq 4 ]
then
    if [ $GP_MIN -lt 6 ]
    then
	echo At least v4.6 of gnuplot needed.  You have v${GP_MAJ}.$GP_MIN
	exit
    fi
fi

if [ "x$1" = "x" ] 
then
    echo "$0 needs directory containing tdd files as argument"
    exit
fi

cd $1
step=`ls -1 *tdd |head -2|tail -1|cut -c1-6|sed s/^0*//`
end=`ls -1 *tdd |tail -1|cut -c1-6|sed s/^0*//`
oname=`echo $PWD|awk -F/ '{print $NF}'`
echo "Converting $(( $end / $step )) tdd files"
# run Cuts program on all files
ls -1 *tdd >$$-tmp
$CUTS <$$-tmp >/dev/null
rm $$-tmp 
echo "Running through gnuplot -> png"
cat >gnu$$.gnu <<EOF
#!/usr/bin/env gnuplot
# gnuplot script to create animation from bunch of *rxy.txt files as output 
# from Sky3d's Cuts program

set zrange [0:0.2]
set contour base
set cntrparam levels incr 0.02,0.02,0.14
unset surface
set size ratio -1
set view map
set key rmargin

do for [it=0:$end:$step] {
  fname=sprintf('%06irxz.txt',it)	
  set term png size 640,640
  set out sprintf('%06irxz.png',it)
  set title sprintf('$2 Iteration # %06i ',it)
  splot fname w l t 'density'
}
EOF

gnuplot gnu$$.gnu

 rm gnu$$.gnu *r??.txt 
echo "now making .avi"
# with the png files created, use mencoder to make a movie
[ -f $oname.avi ] && rm -f $oname.avi
mencoder mf://*r??.png -mf fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -o $oname.avi >/dev/null 2>&1
rm *r??.png

# make some other useful formats
echo "now making .mp4"
[ -f $oname.mp4 ] && rm -f $oname.mp4
ffmpeg -i $oname.avi $oname.mp4 >/dev/null 2>&1

if [ -f $oname.mp4 ] 
then
