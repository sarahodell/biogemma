#!/bin/bash

#iget -K -P -T --retries 3 --lfrestart checkpoint-file /iplant/home/shared/panzea/genotypes/GBS/v27/ZeaGBSv27_publicSamples_rawGenos_AGPv3_20170206.h5
#iget -K -P -T --retries 3 --lfrestart checkpoint-file /iplant/home/shared/panzea/genotypes/GBS/v27/AllZeaGBSv2.7impV5_AnonDonors4k.tar.gz
#iget -K -P -T --retries 3 --lfrestart checkpoint-file /iplant/home/shared/panzea/genotypes/GBS/v27/ZeaGBSv27_publicSamples_imputedV5_AGPv3_20170206.h5

iget -K -P -T -r --retries 3 -X checkpoint-file2 /iplant/home/shared/panzea/hapmap3/hmp321/imputed/uplifted_APGv4
