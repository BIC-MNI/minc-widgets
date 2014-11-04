# CMakeFiles.txt for the minc-widgets
#
# 
#
# Andrew Janke - a.janke@gmail.com
# Vladimir S. FONOV - vladimir.fonov@gmail.com

PROJECT(minc-widgets)

INSTALL(PROGRAMS ecat2hdf5/ecat2hdf5
    gennlxfm/gennlxfm
    img2mnc/img2mnc
    mincbigaverage/mincbigaverage
    mincnorm/mincnorm
    nlpfit/nlpfit
    volalign/volalign
    volcentre/volcentre
    volextents/volextents
    volflip/volflip
    voliso/voliso
    volmash/volmash
    volpad/volpad
    volrot/volrot
    volsymm/volsymm
    xfmavg/xfmavg
    xfmdecomp/xfmdecomp.pl
    xfmflip/xfmflip
  PERMISSIONS
    OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ  WORLD_EXECUTE WORLD_READ
  DESTINATION
    bin )

