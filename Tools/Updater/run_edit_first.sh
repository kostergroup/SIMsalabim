#!/bin/bash


#SimSS and ZimT:
for parfile in ../../*/simulation_setup.txt     # list directories in the form "/tmp/dirname/"
do
	./Updater -f $parfile	
done


#now the tests:
for parfile in ../../*/*/simulation_setup.txt     # list directories in the form "/tmp/dirname/"
do
	./Updater -f $parfile	
done

