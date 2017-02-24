#!/bin/bash
SERIE_DIR=$1

BINS="/usr/lib64/ralph/filters"

echo "0 - Copying DICOM serie to a single file"
if [ ! -d source.dcm ]; then 
    $BINS/dicom_serie2image $SERIE_DIR source.dcm

    # http://doc1.cima.es/bugzilla/show_bug.cgi?id=13
    echo "0bis - Inverting the image to let the algorithms work (feel guilty!)"
    $BINS/inv_image source.dcm 0-input.dcm

    echo "1 - Applying median filter to smooth out small objects"
    $BINS/apply_median 0-input.dcm 1 1-median.dcm

    echo "2 - Thresholding with Hu"
    $BINS/hu_threshold 1-median.dcm 2-hu.dcm

    echo "3 - Converting to a LabelMap keeping only the 3 biggest objects (with luck background, lungs and bed )"
    $BINS/labelize 2-hu.dcm 3-map.dcm

    echo "4 - Choosing lungs from candidates objects"
    $BINS/choose_lungs 3-map.dcm 4-lungs_mask.dcm

    echo "5 - Extracting lungs from original image:"
    $BINS/apply_mask 0-input.dcm 4-lungs_mask.dcm 5-lungs.dcm

fi

#HU_air_real=$($BINS/read_air_calibration 0-input.dcm 2-hu.dcm|awk '{print $3}')
#echo "air_calibration is: "$HU_air_real

#echo "6 - Checking blood calibration:"
#HU_blood_real=$($BINS/read_aorta_calibration 0-input.dcm 2-hu.dcm|awk '{print $3}')

#echo "7 - Recalibrate if needed:"
#$BINS/calibrate 0-input.dcm $HU_air_real $HU_blood_real
