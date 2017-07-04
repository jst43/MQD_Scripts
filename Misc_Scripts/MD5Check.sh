#!/bin/bash

if [ $# -eq 0 ]; then
	echo "No arguments provided"
	exit 1
fi

while getopts ":f:h" opt; do
	case $opt in
		f)
			Filepath=$OPTARG >&2
			echo "Filepath is $OPTARG" >&2
			;;
		h)
			echo "Usage $0 [ -f FILEPATH ]" >&2
			echo
			echo " -f	filepath to directory containing files and their"
			echo "		corresponding md5 files"
			echo " -h	display this help message"
			exit 1
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument" >&2
			exit 1
			;;

	esac
done

cd $Filepath

ls *.md5 > md5_list.txt

sed -i 's|.md5||' md5_list.txt


