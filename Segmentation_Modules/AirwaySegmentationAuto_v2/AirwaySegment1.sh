#UNIX shell executable file to to segment airways lungs with rule-based wavefront propagation
imagedir="/home/xabiarta/Desktop/Images/Mice/LungAtlas"

cd "/home/xabiarta/ITK/AirwaySegmentation_v4/bin"

./Tree_FM "$imagedir/Atlas_001/Atlas_001_crop_closing.tif" 147 117 452 "$imagedir/Atlas_001/Atlas_001_airways_segmented_FM.tif"

./Tree_v4 "$imagedir/Atlas_002/Atlas_002_crop_closing.tif" 174 145 445 "$imagedir/Atlas_002/Atlas_002_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_003/Atlas_003_crop_closing.tif" 181 131 403 "$imagedir/Atlas_003/Atlas_003_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_004/Atlas_004_crop_closing.tif" 160 110 397 "$imagedir/Atlas_004/Atlas_004_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_005/Atlas_005_crop_closing.tif" 178 135 387 "$imagedir/Atlas_005/Atlas_005_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_006/Atlas_006_crop_closing.tif" 152 115 379 "$imagedir/Atlas_006/Atlas_006_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_007/Atlas_007_crop_closing.tif" 216 104 383 "$imagedir/Atlas_007/Atlas_007_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_008/Atlas_008_crop_closing.tif" 222 148 421 "$imagedir/Atlas_008/Atlas_008_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_009/Atlas_009_crop_closing.tif" 149 125 409 "$imagedir/Atlas_009/Atlas_009_airways_segmented_v4.tif"

./Tree_v4 "$imagedir/Atlas_010/Atlas_010_crop_closing.tif" 140 120 430 "$imagedir/Atlas_010/Atlas_010_airways_segmented_v4.tif"







