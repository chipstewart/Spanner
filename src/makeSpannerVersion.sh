REV=`svnversion -n`
# info | grep Revision | awk ''{print $2}''`
#version -n`
echo "#define SpannerSvnVersion \"5.$REV\"" > SpannerVersion.h
