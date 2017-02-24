#!/bin/bash

if [ $# -eq 0 ];
then 
   echo "./generatePersistence.sh TestNumber: (1, 2, or 3)";
else 
#   odb -d mysql --generate-query --generate-schema Measure.hxx;
#   odb -d mysql --generate-query --generate-schema ImageEntity.hxx;
   g++ -c Measure-odb.cxx;
   g++ -c ImageEntity-odb.cxx;
#   g++ -c testPersistence.cxx;
#   g++ -o testPersistence testPersistence.o Measure-odb.o ImageEntity-odb.o -lodb-mysql -lodb;
#   ./testPersistence --user odb_test --database odb_test $1 ;
fi
