#UNIX shell executable file to to segment airways lungs with rule-based wavefront propagation

cd "/home/xabiarta/ITK/AirwaySegmentation_v4/bin"
imagedir="/home/xabiarta/Desktop/Images/Mice/LungAtlas"


./Tree_v4 "$imagedir/Atlas_011/Atlas_011_crop_closing.tif" 162 74 401 "$imagedir/Atlas_011/Atlas_011_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_012/Atlas_012_crop_closing.tif" 147 107 410 "$imagedir/Atlas_012/Atlas_012_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_013/Atlas_013_crop_closing.tif" 177 149 446 "$imagedir/Atlas_013/Atlas_013_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_014/Atlas_014_crop_closing.tif" 201 122 442 "$imagedir/Atlas_014/Atlas_014_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_015/Atlas_015_crop_closing.tif" 225 111 364 "$imagedir/Atlas_015/Atlas_015_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_016/Atlas_016_crop_closing.tif" 164 140 476 "$imagedir/Atlas_016/Atlas_016_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_017/Atlas_017_crop_closing.tif" 188 138 429 "$imagedir/Atlas_017/Atlas_017_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_018/Atlas_018_crop_closing.tif" 177 152 375 "$imagedir/Atlas_018/Atlas_018_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_019/Atlas_019_crop_closing.tif" 165 164 437 "$imagedir/Atlas_019/Atlas_019_airways_segmented_v4.tif"
