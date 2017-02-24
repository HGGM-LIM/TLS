bin_dir="/usr/lib64/ralph/filters"
converter="ImageReadDicomSeriesWrite"

KERTYPE="NOT-KNOWN"
SL="1" #We should read it from the original image header
IDS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30"

scan_dir="Online/scans"

read -ra TYPES <<< "Fixed Moving"

for PID in $IDS; do

    echo "Processing "$PID
    for TYPE in "${TYPES[@]}"; do
	image=$scan_dir/$PID"_"$TYPE".mhd"
	echo $image

	mkdir -p dicoms/$PID/$TYPE

	$bin_dir/$converter $image dicoms/$PID/$TYPE

	#Set Class ID to Enhanced CT Media Storage
        # We keep have problem with spacing
        #dcmodify -nb -v -m "SOPClassUID=1.2.840.10008.5.1.4.1.1.2.1" dicoms/$new_image

        #Set parameters
	# (0018,1210) - ReconstructionKernel
        # (0018,0050) - SliceThickness
        # (0010,0020) - PatientID
	dcmodify -nb -i "(0018,1210)=$KERTYPE" dicoms/$PID/$TYPE/*
        dcmodify -nb -i "(0018,0050)=$SL" dicoms/$PID/$TYPE/*
        dcmodify -nb -i "(0010,0020)=$PID" dicoms/$PID/$TYPE/*

	#gdcmdump dicoms/$new_image > $PID_final.txt

        # Then you can check differences with:
        #diff -y 01_original.txt 01_final.txt
    done
    
done