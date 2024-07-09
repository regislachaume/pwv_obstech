#! /bin/sh

for ubx in *.ubx; 
do
    rnx=`echo $ubx | sed "s/.ubx/.rnx/"`
    convbin -ht GEODETIC -ha "UNKNWOWN/TWIVSP6037L" -ho "El Sauce/Obstech" -hr "UNKNOWN/UBLOX-ZED-F9P" -od -os -v 3 -hm SAUC00CHL -hp "1813185/-5195981/-3216438" -o $rnx $ubx; 
done
