########################################
#                                      #
#  Script to segment control mice      #
#                                      #
########################################

ITKdir="/home/xabiarta/ITK"
FastMarchingDir="$ITKdir/AirwaySegmentationFM_v2/bin"
imagedir="/home/xabiarta/Desktop/Images/Mice/LungAtlas"
cd $FastMarchingDir

./Tree_FM_v2 "$imagedir/w04_C501/w04_C501_HU.mhd" 230 221 0 "$imagedir/w04_C501/w04_C501_airways_v2.tiff"

